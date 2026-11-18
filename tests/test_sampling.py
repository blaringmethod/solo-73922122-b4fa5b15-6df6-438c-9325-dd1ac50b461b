"""Snapshot sampling API: allocation rules, determinism, integrity checks."""
import gzip
import io
import json
import time
from pathlib import Path

import pytest

from app import sampling

BARCODES = [
    {"sample_id": "s1", "barcode_r1": "AAAAAAAA", "barcode_r2": "CCCCCCCC"},
    {"sample_id": "s2", "barcode_r1": "GGGGGGGG", "barcode_r2": "TTTTTTTT"},
]

R1_BAR = b"AAAAAAAA"
R2_BAR = b"CCCCCCCC"


def _records(names, insert1=b"ACGTACGT", insert2=b"TGCATGCA"):
    """Paired records whose barcodes hit sample s1 exactly."""
    r1 = [(n, R1_BAR + insert1, b"I" * (8 + len(insert1))) for n in names]
    r2 = [(n, R2_BAR + insert2, b"J" * (8 + len(insert2))) for n in names]
    return r1, r2


def _post(client, lane, barcodes, r1_path, r2_path):
    with open(r1_path, "rb") as f1, open(r2_path, "rb") as f2:
        return client.post(
            "/jobs",
            data={"lane": lane, "barcodes": json.dumps(barcodes)},
            files={"r1": ("R1.fastq.gz", f1, "application/gzip"),
                   "r2": ("R2.fastq.gz", f2, "application/gzip")},
        )


def _wait(client, job_id, timeout=15):
    deadline = time.time() + timeout
    while time.time() < deadline:
        row = client.get(f"/jobs/{job_id}").json()
        if row["status"] in ("published", "failed"):
            return row
        time.sleep(0.05)
    raise AssertionError("timeout waiting for job")


def publish_lane(client, helpers, tmp_path, lane, names, barcodes=BARCODES,
                 insert1=b"ACGTACGT", insert2=b"TGCATGCA"):
    """Publish one lane whose s1 bucket holds the given read names."""
    d = tmp_path / lane
    d.mkdir()
    r1, r2 = _records(names, insert1, insert2)
    p1, p2 = helpers.make_pair(r1, r2, d)
    resp = _post(client, lane, barcodes, p1, p2)
    assert resp.status_code == 201, resp.text
    done = _wait(client, resp.json()["job"]["job_id"])
    assert done["status"] == "published", done
    return lane


def make_snapshot(client, helpers, tmp_path, lane_counts):
    """Publish lanes (snapshot order) with n s1-pairs each; return snapshot."""
    for i, n in enumerate(lane_counts):
        publish_lane(client, helpers, tmp_path, f"lane{i}",
                     [f"r{k}".encode() for k in range(n)])
    return client.post("/snapshots").json()


def sample(client, snap_id, sample_id="s1", seed="seed-1", pairs=10):
    return client.get(
        f"/snapshots/{snap_id}/samples/{sample_id}/sample",
        params={"seed": seed, "pairs": pairs})


def parse_fastq(text):
    lines = text.split("\n")
    assert lines[-1] == ""  # trailing newline after the last record
    lines = lines[:-1]
    assert len(lines) % 4 == 0
    out = []
    for i in range(0, len(lines), 4):
        out.append({
            "name": lines[i],
            "seq": lines[i + 1],
            "plus": lines[i + 2],
            "qual": lines[i + 3],
        })
    return out


# ---------------- allocation unit tests ----------------

def test_alloc_every_nonempty_lane_gets_one_when_slots_suffice():
    assert sampling.allocate([4, 6, 10], 3) == [1, 1, 1]
    assert sampling.allocate([4, 6, 10], 4) == [2, 1, 1]  # leftover in order


def test_alloc_covers_earliest_lanes_when_slots_short():
    assert sampling.allocate([5, 5, 5], 2) == [1, 1, 0]
    assert sampling.allocate([0, 5, 0, 5, 5], 2) == [0, 1, 0, 1, 0]
    assert sampling.allocate([5, 5, 5], 1) == [1, 0, 0]


def test_alloc_proportional_remainder_and_capacity_fill():
    # 7 spare over counts 4/6/10: floors 1/2/3 -> 2/3/4, leftover 1 to lane0.
    assert sampling.allocate([4, 6, 10], 10) == [3, 3, 4]
    # Lane0/1 saturate at 2 each; their shortfall flows to lane2 in order.
    assert sampling.allocate([2, 2, 100], 60) == [2, 2, 56]
    # Heavy shortfall: small lanes fill completely, the rest goes downstream.
    assert sampling.allocate([1, 1, 8], 10) == [1, 1, 8]


def test_alloc_total_is_min_of_target_and_available():
    assert sampling.allocate([3, 4], 200) == [3, 4]
    assert sampling.allocate([0, 0], 5) == [0, 0]
    assert sampling.allocate([3, 4], 0) == [0, 0]
    got = sampling.allocate([9, 1, 9], 11)
    assert sum(got) == 11 and got[1] == 1


def test_alloc_never_exceeds_lane_capacity():
    for counts, target in [([1, 2, 3], 6), ([2, 2, 2], 5), ([7, 0, 3], 9),
                           ([1, 1, 1, 1], 4), ([5, 6], 11)]:
        got = sampling.allocate(counts, target)
        assert all(g <= c for g, c in zip(got, counts))
        assert sum(got) == min(target, sum(counts))


# ---------------- API behaviour ----------------

def test_sampling_deterministic_and_proportional(client, helpers, tmp_path):
    snap = make_snapshot(client, helpers, tmp_path, [4, 6, 10])
    r1 = sample(client, snap["snapshot_id"], pairs=10)
    assert r1.status_code == 200, r1.text
    body = r1.json()
    # Same request again -> byte-identical response.
    r2 = sample(client, snap["snapshot_id"], pairs=10)
    assert r2.text == r1.text

    assert body["requested_pairs"] == 10
    assert body["returned_pairs"] == 10
    assert body["available_pairs"] == 20
    got = {l["lane"]: l["selected_pairs"] for l in body["lanes"]}
    assert got == {"lane0": 3, "lane1": 3, "lane2": 4}  # allocate([4,6,10],10)
    avail = {l["lane"]: l["available_pairs"] for l in body["lanes"]}
    assert avail == {"lane0": 4, "lane1": 6, "lane2": 10}

    # R1/R2 texts hold the same pairs in the same order; names match.
    r1_recs = parse_fastq(body["r1_fastq"])
    r2_recs = parse_fastq(body["r2_fastq"])
    assert len(r1_recs) == len(r2_recs) == 10
    for a, b in zip(r1_recs, r2_recs):
        assert a["name"] == b["name"]
        assert a["plus"] == "+" and b["plus"] == "+"
        assert len(a["seq"]) == len(a["qual"])
    # Original order within each lane: names are r<index>, ascending per lane.
    names = [r["name"][1:] for r in r1_recs]
    assert names[:3] == sorted(names[:3], key=lambda s: int(s[1:]))
    assert names[3:6] == sorted(names[3:6], key=lambda s: int(s[1:]))
    assert names[6:] == sorted(names[6:], key=lambda s: int(s[1:]))
    # Trimmed inserts (barcode removed at demux time).
    assert r1_recs[0]["seq"] == "ACGTACGT"
    assert r2_recs[0]["seq"] == "TGCATGCA"


def test_sampling_covers_first_lanes_when_target_small(
        client, helpers, tmp_path):
    snap = make_snapshot(client, helpers, tmp_path, [5, 5, 5])
    body = sample(client, snap["snapshot_id"], pairs=2).json()
    got = {l["lane"]: l["selected_pairs"] for l in body["lanes"]}
    assert got == {"lane0": 1, "lane1": 1, "lane2": 0}
    assert body["returned_pairs"] == 2


def test_sampling_total_capped_by_available(client, helpers, tmp_path):
    snap = make_snapshot(client, helpers, tmp_path, [3, 4])
    body = sample(client, snap["snapshot_id"], pairs=200).json()
    assert body["returned_pairs"] == 7
    got = {l["lane"]: l["selected_pairs"] for l in body["lanes"]}
    assert got == {"lane0": 3, "lane1": 4}


def test_sampling_seed_changes_draw_but_stays_deterministic(
        client, helpers, tmp_path):
    snap = make_snapshot(client, helpers, tmp_path, [60])
    a = sample(client, snap["snapshot_id"], seed="alpha", pairs=10).json()
    b = sample(client, snap["snapshot_id"], seed="beta", pairs=10).json()
    a_again = sample(client, snap["snapshot_id"], seed="alpha", pairs=10).json()
    assert a == a_again
    names_a = [r["name"] for r in parse_fastq(a["r1_fastq"])]
    names_b = [r["name"] for r in parse_fastq(b["r1_fastq"])]
    assert names_a != names_b  # 60 choose 10: collision is negligible


def test_sampling_duplicate_names_each_take_a_slot(client, helpers, tmp_path):
    names = [b"dup", b"dup", b"dup", b"solo"]
    publish_lane(client, helpers, tmp_path, "dupLane", names)
    snap = client.post("/snapshots").json()
    body = sample(client, snap["snapshot_id"], pairs=4).json()
    assert body["available_pairs"] == 4
    assert body["returned_pairs"] == 4
    names_out = [r["name"] for r in parse_fastq(body["r1_fastq"])]
    assert names_out == ["@dup", "@dup", "@dup", "@solo"]  # original order


def test_sampling_excludes_lanes_published_after_snapshot(
        client, helpers, tmp_path):
    publish_lane(client, helpers, tmp_path, "early", [b"a", b"b"])
    snap = client.post("/snapshots").json()
    publish_lane(client, helpers, tmp_path, "late", [b"c", b"d", b"e"])
    body = sample(client, snap["snapshot_id"], pairs=10).json()
    assert [l["lane"] for l in body["lanes"]] == ["early"]
    assert body["available_pairs"] == 2
    assert body["returned_pairs"] == 2


def test_sampling_empty_sample_404_and_unknown_snapshot_404(
        client, helpers, tmp_path):
    publish_lane(client, helpers, tmp_path, "laneX", [b"a"])
    snap = client.post("/snapshots").json()
    assert sample(client, snap["snapshot_id"], sample_id="ghost"
                  ).status_code == 404
    assert sample(client, "no-such-snapshot").status_code == 404


def test_sampling_pairs_range_validated(client, helpers, tmp_path):
    publish_lane(client, helpers, tmp_path, "laneY", [b"a"])
    snap = client.post("/snapshots").json()
    assert sample(client, snap["snapshot_id"], pairs=0).status_code == 422
    assert sample(client, snap["snapshot_id"], pairs=201).status_code == 422
    resp = client.get(
        f"/snapshots/{snap['snapshot_id']}/samples/s1/sample",
        params={"pairs": 1})
    assert resp.status_code == 422  # seed is required
    assert sample(client, snap["snapshot_id"], pairs=1).status_code == 200
    assert sample(client, snap["snapshot_id"], pairs=200).status_code == 200


def test_sampling_skips_lanes_where_sample_absent(client, helpers, tmp_path):
    # laneA: s1 + s2; laneB: only s2 in its barcode table.
    publish_lane(client, helpers, tmp_path, "laneA", [b"a", b"b"])
    only_s2 = [BARCODES[1]]
    d = tmp_path / "laneB"
    d.mkdir()
    r1 = [(b"x", b"GGGGGGGGACGT", b"I" * 12)]
    r2 = [(b"x", b"TTTTTTTTTGCA", b"J" * 12)]
    p1, p2 = helpers.make_pair(r1, r2, d)
    resp = _post(client, "laneB", only_s2, p1, p2)
    assert _wait(client, resp.json()["job"]["job_id"])["status"] == "published"
    snap = client.post("/snapshots").json()
    body = sample(client, snap["snapshot_id"], pairs=10).json()
    assert [l["lane"] for l in body["lanes"]] == ["laneA"]
    assert body["available_pairs"] == 2


# ---------------- integrity failures: the whole request fails ----------------

def _lane_file(data_dir, lane, sample_id, mate):
    return data_dir / "published" / "lanes" / lane / \
        f"{sample_id}.{mate}.fastq.gz"


def test_sampling_missing_file_fails(client, helpers, tmp_path, data_dir):
    snap = make_snapshot(client, helpers, tmp_path, [3])
    _lane_file(data_dir, "lane0", "s1", "r1").unlink()
    resp = sample(client, snap["snapshot_id"], pairs=2)
    assert resp.status_code == 409
    assert "missing" in resp.json()["detail"].lower()


def test_sampling_corrupt_gzip_fails(client, helpers, tmp_path, data_dir):
    snap = make_snapshot(client, helpers, tmp_path, [3])
    victim = _lane_file(data_dir, "lane0", "s1", "r2")
    victim.write_bytes(victim.read_bytes()[:-8])  # truncate gzip stream
    resp = sample(client, snap["snapshot_id"], pairs=2)
    assert resp.status_code == 409


def test_sampling_count_mismatch_fails(client, helpers, tmp_path, data_dir):
    snap = make_snapshot(client, helpers, tmp_path, [3])
    from app import db as dbmod
    dbmod.get_conn().execute(
        "UPDATE lane_counts SET read_pairs=99 WHERE lane='lane0' "
        "AND sample_id='s1'")
    resp = sample(client, snap["snapshot_id"], pairs=2)
    assert resp.status_code == 409
    assert "count" in resp.json()["detail"].lower()


def test_sampling_oversize_sequence_line_fails(client, helpers, tmp_path):
    # 8-base barcode + 1001-base insert -> published seq line is 1001 bytes.
    publish_lane(client, helpers, tmp_path, "longSeq", [b"a"],
                 insert1=b"A" * 1001)
    snap = client.post("/snapshots").json()
    resp = sample(client, snap["snapshot_id"], pairs=1)
    assert resp.status_code == 409
    assert "exceeds" in resp.json()["detail"]


def test_sampling_oversize_read_name_fails(client, helpers, tmp_path):
    publish_lane(client, helpers, tmp_path, "longName", [b"n" * 300])
    snap = client.post("/snapshots").json()
    resp = sample(client, snap["snapshot_id"], pairs=1)
    assert resp.status_code == 409
    assert "read name" in resp.json()["detail"]


def test_sampling_limits_accept_exact_boundary(client, helpers, tmp_path):
    # Insert of 1000 bases -> published seq line is exactly 1000 bytes.
    publish_lane(client, helpers, tmp_path, "edgeLane",
                 [b"n" * 256], insert1=b"A" * 1000, insert2=b"T" * 1000)
    snap = client.post("/snapshots").json()
    resp = sample(client, snap["snapshot_id"], pairs=1)
    assert resp.status_code == 200, resp.text
    rec = parse_fastq(resp.json()["r1_fastq"])[0]
    assert len(rec["seq"]) == 1000
    assert rec["name"] == "@" + "n" * 256


def test_sampling_mate_name_mismatch_fails(client, helpers, tmp_path,
                                           data_dir):
    snap = make_snapshot(client, helpers, tmp_path, [2])
    # Rewrite the published R1 with a renamed second record (desync).
    victim = _lane_file(data_dir, "lane0", "s1", "r1")
    with gzip.open(victim, "rb") as g:
        content = g.read()
    content = content.replace(b"@r1\n", b"@zzz\n")
    with gzip.open(victim, "wb") as g:
        g.write(content)
    resp = sample(client, snap["snapshot_id"], pairs=2)
    assert resp.status_code == 409


# ---------------- line-limit reader unit tests ----------------

def test_bounded_readline_rejects_long_lines(tmp_path):
    from app import processing
    r1 = tmp_path / "R1.fastq.gz"
    r2 = tmp_path / "R2.fastq.gz"
    with gzip.open(r1, "wb") as g:
        g.write(b"@x\n" + b"A" * 5000 + b"\n+\n" + b"I" * 5000 + b"\n")
    with gzip.open(r2, "wb") as g:
        g.write(b"@x\nCCCC\n+\nIIII\n")
    with pytest.raises(processing.FastqFormatError, match="exceeds 1000"):
        list(processing._paired_records(r1, r2, size_budget=1 << 20,
                                        line_limit=1000))
    # Without the limit the same file parses fine (demux path unchanged).
    pairs = list(processing._paired_records(r1, r2, size_budget=1 << 30))
    assert len(pairs) == 1 and len(pairs[0][0].seq) == 5000
