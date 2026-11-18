# Sequencing Lane Demultiplexing API

FastAPI + SQLite service that splits mixed paired-end gzip FASTQ uploads into
per-sample files. Designed to run as **two container instances on one shared
volume** with a **128 MiB memory limit**.

## What it does

* Accepts a lane number, a dual-index barcode table (1–8 samples, two
  8-base barcodes each) and paired `R1`/`R2` gzip FASTQ files.
  Decompressed input across both mates is capped at **512 MiB**.
* Validates both mates: valid gzip, well-formed FASTQ, equal record counts,
  identical read names in the same order, matching sequence/quality lengths.
  Any corruption fails the **whole lane**.
* Assigns a pair to a sample when the combined R1+R2 barcode matches with
  **at most one substitution** and the candidate sample is **unique**.
  Everything else (no hit, >1 substitution, or a tie) goes to quarantine.
* Assigned reads have both 8-base barcodes and their corresponding quality
  characters trimmed; quarantined reads are kept intact.
* A lane becomes visible only when **all** per-sample paired files, the
  quarantine paired files and the count table exist (atomic publish:
  files staged, renamed, then a single DB transaction).
* Same lane + identical input (hash of barcode table and both files) reuses
  the existing task; same lane with different input is rejected with `409`.
* Failed tasks can be retried (`POST /jobs/{id}/retry`); uploads are kept
  until a lane is published.
* Snapshots freeze the published lane list at creation time. Sample/
  quarantine downloads through a snapshot only ever concatenate files from
  that fixed list — lanes published later never mix in.
* A **sampling endpoint** draws a deterministic subset of a sample's read
  pairs from a snapshot for trial alignment, spreading the quota across
  lanes so high-throughput lanes cannot hog it (see below).

## Snapshot sampling (trial alignment)

`GET /snapshots/{id}/samples/{sample_id}/sample?seed=...&pairs=1..200`

Draws up to `pairs` read pairs of one sample from the snapshot's frozen
lane list and returns per-lane hit counts plus the paired-end FASTQ text
(`r1_fastq` / `r2_fastq`, positionally paired, lanes in snapshot order,
pairs in original file order within each lane).

Allocation of the pair quota:

* If the target covers every non-empty lane, each contributes **at least
  one pair**; otherwise the earliest lanes in snapshot order are covered
  with one pair each until the target is spent.
* The remaining quota is split **in proportion to each lane's published
  read-pair count**; rounding remainders and shares a lane cannot absorb
  (capacity shortfall) are filled in snapshot order.
* The grand total is `min(target, pairs actually available)`.

Within a lane the draw is an index sample keyed by a per-lane hash of the
request `seed`, so **identical requests return identical results** (verified
byte-identical across both container instances) and duplicate read names in
one lane each occupy their own slot.

Integrity: every contributing lane is fully re-scanned — read names, paired
format and the published pair count are verified, and each mate's sequence/
quality lines must be ≤ 1000 bytes and read names ≤ 256 bytes. (Record
lines, headers included, are read with a 1000-byte bound so a pathological
single-line file cannot be buffered whole.) A missing file, corruption,
count mismatch or over-limit line fails the **whole request** with `409`. Per-lane decompressed data may reach 512 MiB and is
streamed record-by-record (bounded line reads), so the 128 MiB container
limit is respected; the scan runs in the threadpool, so other lanes keep
being served meanwhile.

## Running (Docker Compose, two instances on a shared volume)

```bash
docker compose up --build
# instance A: http://localhost:8001   instance B: http://localhost:8002
# OpenAPI demo UI:  http://localhost:8001/docs  (also on 8002)
```

Both containers mount the named volume `demux-shared` at `/data`, are capped
at 128 MiB (`mem_limit`), and automatically recover on restart.

### Quick demo

```bash
python scripts/make_demo_fastq.py /tmp/demo          # makes R1/R2.fastq.gz
# In http://localhost:8001/docs -> POST /jobs:
#   lane: demoLane
#   barcodes: [{"sample_id":"sampleA","barcode_r1":"AAAAAAAA","barcode_r2":"CCCCCCCC"},
#              {"sample_id":"sampleB","barcode_r1":"GGGGGGGG","barcode_r2":"TTTTTTTT"}]
#   r1/r2: the generated gzip files
curl http://localhost:8001/lanes/demoLane             # counts
curl -OJ http://localhost:8001/lanes/demoLane/samples/sampleA/r1
curl -OJ http://localhost:8001/lanes/demoLane/quarantine/r1

# Deterministic sampling for trial alignment (also try on :8002 — identical):
SNAP=$(curl -s -X POST http://localhost:8001/snapshots | python -c 'import json,sys; print(json.load(sys.stdin)["snapshot_id"])')
curl -s "http://localhost:8001/snapshots/$SNAP/samples/sampleA/sample?seed=run1&pairs=5"
```

## API surface

| Method | Path | Purpose |
|---|---|---|
| POST | `/jobs` | Submit a lane (multipart: `lane`, `barcodes` JSON, `r1`, `r2`) |
| GET | `/jobs`, `/jobs/{id}` | Task status (`queued`/`processing`/`published`/`failed`) |
| POST | `/jobs/{id}/retry` | Requeue a failed task |
| GET | `/lanes`, `/lanes/{lane}` | Published lanes and per-bucket counts |
| GET | `/lanes/{lane}/samples/{sample_id}/{r1\|r2}` | Download a sample mate |
| GET | `/lanes/{lane}/quarantine/{r1\|r2}` | Download a quarantine mate |
| POST | `/snapshots` | Freeze the current published lane list |
| GET | `/snapshots`, `/snapshots/{id}` | List/inspect snapshots |
| GET | `/snapshots/{id}/samples/{sample_id}/{mate}` | Snapshot-scoped concatenated download |
| GET | `/snapshots/{id}/samples/{sample_id}/sample` | Deterministic read-pair sampling (`seed`, `pairs` 1–200) |
| GET | `/snapshots/{id}/quarantine/{mate}` | Snapshot-scoped quarantine download |
| GET | `/health` | Liveness, instance id, current job, queue depth |

## Concurrency and crash safety

* **SQLite WAL** with `BEGIN IMMEDIATE` job claiming: exactly one instance
  ever picks up a given job; the other keeps serving requests.
* Each instance processes one lane at a time in a background thread, so a
  slow lane never blocks uploads, queries or downloads on the event loop.
* A job runs under an exclusive `flock` on its work directory. Startup
  recovery removes scratch/staging debris only when no lock is held, and
  never deletes directories that have a published-lane DB row — so another
  instance's in-flight work and all published data are preserved.
* Publish order: write/verify outputs in `work/`, move to
  `published/lanes/.tmp.<job>`, atomic `rename` to the final lane directory,
  then insert metadata in one DB transaction. Inputs are deleted only after
  the commit, so a crash at any point leaves no half-visible lane.

## Memory

All I/O is streamed record-by-record (1 MiB network chunks, line-based gzip
reads, one open gzip writer per bucket). No file is ever buffered whole.
The automated streaming test processes hundreds of MiB of decompressed
data with RSS growth below 64 MiB.

## Layout

```
app/            FastAPI app, SQLite layer, streaming splitter, worker,
                snapshot sampling (allocation + seeded draw)
scripts/        demo FASTQ generator, two-instance integration check
tests/          pytest suite (happy path, corruption, retry, recovery,
                barcode rules, size budget, streaming memory, sampling)
Dockerfile      python:3.11-slim, non-root user, /data volume, healthcheck
docker-compose.yml  two 128 MiB instances sharing one named volume
```

## Tests

```bash
python -m pytest tests/ -q
python scripts/integration_two_instances.py   # real two-process check
```
