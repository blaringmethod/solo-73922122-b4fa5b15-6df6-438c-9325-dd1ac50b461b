"""Snapshot sampling: deterministic paired-end read-pair draws.

Picks up to ``target`` read pairs of one sample from the frozen lane list of
a snapshot, for trial alignment. Allocation avoids letting high-throughput
lanes dominate:

* when the target covers every non-empty lane, each such lane contributes at
  least one pair; otherwise the earliest lanes in snapshot order are covered
  (one pair each) until the target is spent;
* the remaining quota is distributed in proportion to each lane's published
  read-pair count; rounding remainders and shares a lane cannot absorb
  (capacity shortfall) are filled in snapshot order;
* the grand total is min(target, pairs actually available).

Within a lane the draw is an index sample driven by a per-lane key derived
from the request seed, so identical requests return identical results and
duplicate read names in one lane each occupy their own slot. Selected pairs
are emitted in original file order, lanes in snapshot order.

Every lane that contributes pairs is fully re-scanned: read names, paired
format, per-line limits and the published pair count are all verified. Any
missing file, corruption or count mismatch fails the whole request. All I/O
is streamed record by record, so a 512 MiB lane fits the 128 MiB container.
"""
from __future__ import annotations

import hashlib
import random
from pathlib import Path
from typing import Any

from . import config, db, processing


class SampleDataError(Exception):
    """Published data failed integrity checks; the whole request fails."""


def allocate(counts: list[int], target: int) -> list[int]:
    """Split ``target`` pairs across lanes (snapshot order) per the rules."""
    alloc = [0] * len(counts)
    nonempty = [i for i, c in enumerate(counts) if c > 0]
    if target <= 0 or not nonempty:
        return alloc
    goal = min(target, sum(counts))
    if target < len(nonempty):
        # Not enough slots: cover as many lanes as possible, snapshot order.
        for i in nonempty[:target]:
            alloc[i] = 1
        return alloc
    # Every non-empty lane gets its guaranteed pair first.
    for i in nonempty:
        alloc[i] = 1
    remaining = goal - len(nonempty)
    if remaining <= 0:
        return alloc
    # Proportional to published read-pair counts, capped by spare capacity.
    total = sum(counts)
    for i in nonempty:
        share = remaining * counts[i] // total
        alloc[i] += min(share, counts[i] - alloc[i])
    leftover = goal - sum(alloc)
    # Rounding remainders and capacity shortfalls: fill in snapshot order.
    while leftover > 0:
        progressed = False
        for i in nonempty:
            if leftover == 0:
                break
            if alloc[i] < counts[i]:
                alloc[i] += 1
                leftover -= 1
                progressed = True
        if not progressed:  # everything saturated (goal == total available)
            break
    return alloc


def _lane_key(seed: str, lane: str, sample_id: str) -> int:
    """Deterministic per-lane RNG key derived from the request seed."""
    digest = hashlib.sha256(
        f"{seed}\x00{lane}\x00{sample_id}".encode("utf-8")).digest()
    return int.from_bytes(digest[:8], "big")


def select_indices(count: int, k: int, seed: str, lane: str,
                   sample_id: str) -> list[int]:
    """Sorted pair indices drawn within one lane; duplicates keep their slots."""
    rng = random.Random(_lane_key(seed, lane, sample_id))
    return sorted(rng.sample(range(count), k))


def _record_text(rec: processing.Record) -> bytes:
    """Rebuild the original four FASTQ lines of one mate."""
    return rec.name + b"\n" + rec.seq + b"\n" + rec.plus + b"\n" \
        + rec.qual + b"\n"


def _scan_lane(r1_path: Path, r2_path: Path, keep: set[int],
               lane: str) -> tuple[list[tuple[processing.Record,
                                              processing.Record]], int]:
    """Stream one lane, verify it, and collect the selected pairs.

    Returns (selected pairs in file order, total pairs seen). Raises
    SampleDataError on any format, limit or integrity violation.
    """
    selected: list[tuple[processing.Record, processing.Record]] = []
    total = 0
    try:
        pairs = processing._paired_records(
            r1_path, r2_path,
            size_budget=config.MAX_TOTAL_UNCOMPRESSED,
            line_limit=config.SAMPLE_MAX_LINE)
        for index, (r1, r2) in enumerate(pairs):
            for mate, rec in (("R1", r1), ("R2", r2)):
                token = rec.name[1:].split()[0] if len(rec.name) > 1 else b""
                if len(token) > config.SAMPLE_MAX_NAME:
                    raise SampleDataError(
                        f"lane {lane}: record {index} {mate} read name "
                        f"exceeds {config.SAMPLE_MAX_NAME} bytes")
            if index in keep:
                selected.append((r1, r2))
            total += 1
    except (processing.FastqFormatError,
            processing.SizeExceededError) as exc:
        raise SampleDataError(f"lane {lane}: {exc}") from exc
    return selected, total


def _published_file(path_str: str, lane: str) -> Path:
    """Resolve a published lane file, failing the request when absent."""
    root = config.PUBLISHED_DIR.resolve()
    path = Path(path_str).resolve()
    if root not in path.parents or not path.is_file():
        raise SampleDataError(f"lane {lane}: published file missing: "
                              f"{Path(path_str).name}")
    return path


def sample_snapshot(snapshot: dict[str, Any], sample_id: str, seed: str,
                    target: int) -> dict[str, Any] | None:
    """Draw up to ``target`` pairs of ``sample_id`` from the snapshot.

    Returns None when the sample has no pairs in any snapshot lane.
    Raises SampleDataError on missing/corrupt/mismatched published data.
    """
    conn = db.get_conn()
    lanes: list[dict[str, Any]] = []
    for lane in snapshot["lanes"]:
        row = conn.execute(
            "SELECT * FROM lane_counts WHERE lane=? AND sample_id=?",
            (lane, sample_id)).fetchone()
        count = row["read_pairs"] if row is not None else 0
        lanes.append({"lane": lane, "row": row, "count": count})

    counts = [info["count"] for info in lanes]
    available = sum(counts)
    if available == 0:
        return None
    plan = allocate(counts, target)

    per_lane: list[dict[str, Any]] = []
    r1_out: list[bytes] = []
    r2_out: list[bytes] = []
    returned = 0
    for info, want in zip(lanes, plan):
        if info["count"] == 0:
            continue  # empty lane: not a candidate, contributes nothing
        entry = {"lane": info["lane"], "available_pairs": info["count"],
                 "selected_pairs": 0}
        if want > 0:
            row = info["row"]
            r1_path = _published_file(row["r1_file"], info["lane"])
            r2_path = _published_file(row["r2_file"], info["lane"])
            keep = set(select_indices(info["count"], want, seed,
                                      info["lane"], sample_id))
            selected, actual = _scan_lane(r1_path, r2_path, keep, info["lane"])
            if actual != info["count"]:
                raise SampleDataError(
                    f"lane {info['lane']}: published count {info['count']} "
                    f"but files hold {actual} pairs")
            for r1, r2 in selected:
                r1_out.append(_record_text(r1))
                r2_out.append(_record_text(r2))
            entry["selected_pairs"] = len(selected)
            returned += len(selected)
        per_lane.append(entry)

    # latin-1 maps bytes 1:1 into JSON-safe text; ASCII content is untouched.
    return {
        "snapshot_id": snapshot["snapshot_id"],
        "sample_id": sample_id,
        "seed": seed,
        "requested_pairs": target,
        "returned_pairs": returned,
        "available_pairs": available,
        "lanes": per_lane,
        "r1_fastq": b"".join(r1_out).decode("latin-1"),
        "r2_fastq": b"".join(r2_out).decode("latin-1"),
    }
