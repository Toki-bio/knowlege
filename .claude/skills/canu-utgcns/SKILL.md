---
name: canu-utgcns
description: Internals of canu's utgcns (unitig consensus, the 5-consensus stage of canu/HiCanu). Use when utgcns hangs or stalls single-threaded, when debugging consensus parameters (-e, -maxcoverage, -threads), or when explaining why consensus is single-threaded while overlaps finish fast.
---

# canu utgcns internals

> Type: reference (facts verified against the pinned commit) with marked `notes`.

## Pinned source

| Field | Value |
|---|---|
| Repo | https://github.com/marbl/canu |
| Version | master r10522 "canu master +34 changes" (what bioconda `canu 2.3` installs as of 2026-10) |
| Commit | `86d7d3b8816a91372d68eceefb9a59ecfbea3eaa` |
| Digested on | 2026-11-03 |
| Language / build | C/C++ (kmer/merylutils), per-directory Makefile.boilermake |
| License | Public domain (US Gov't work), with GNU-licensed edlib in src/ext/ |

## What it does in three sentences

utgcns is canu's unitig consensus: it takes a partitioned
tig + its reads and produces a FASTA-like `.cns` output (and, on success, feeds
canu's contigs). For each tig it first greedily stitches a linear "template"
from the reads' layout positions, aligning each read to the template tail with
edlib and retrying/backtracking on failure. It then aligns all reads to the
template, builds a partial-order multialignment, and votes per-column.

## Where to look (read on demand)

| Question | File |
|---|---|
| How the stage is driven, what files exist, what `-threads` really does | `references/architecture.md` |
| Step-by-step trace of ONE tig through generatePBDAG, incl. all goto-retry loops | `references/control-flow.md` |
| The known stall: why it runs "forever", single-threaded, on hard unitigs; workarounds | `references/pitfalls.md` |

## Ten facts worth knowing before touching the code

1. `processTigs()` (src/utgcns/utgcns-processTigs.C:123) loops tigs **strictly
   sequentially**; there is no threading inside the consensus compute at all
   (grep for threads in unitigConsensus.C: zero hits). `-threads` is accepted
   (utgcns.C:183) and used only by I/O-side helpers (the 7 extra threads sleep
   while 1 spins - the observed signature of a stuck tig).
2. `cnsParameters::skipTig()` prints the FIRST half of each progress-table row
   when a tig is loaded; processTigs.C:135 prints the second half after
   `generate()` returns. A half-written last line in consensus.NNNNNN.out =
   the tig named at its start never finished (how to name a stalled tig).
3. Template stitching (`generateTemplateStitch`, unitigConsensus.C:386) is a
   backtracking search with `goto retryCandidate`/`retryTig`/`alignAgain`,
   NOT a bounded linear pass. Cost explodes when read placements disagree.
4. edlib calls use EDLIB_MODE_HW + EDLIB_TASK_PATH with allowed distance
   k = olapLen * bandErrRate (unitigConsensus.C:610-612); bandErrRate escalates
   in ERROR_RATE_FACTOR steps up to _errorRate (unitigConsensus.C:702-714).
   edlib cost is ~O(N*k), so every escalation multiplies per-call cost.
5. Whole-tig restart: if reads run out before the tig end, errorRate *= 1.25
   and the ENTIRE template is rebuilt (`retryTig`, unitigConsensus.C:522-528).
6. `badToAdd` (failed-read memo) is CLEARED when `allowLargerShift` flips
   (unitigConsensus.C:748-754), so failed reads are retried under looser
   error rates - backtracking has no monotone progress measure.
7. Reads that never align are appended verbatim to the template
   (unitigConsensus.C:795-824), so even "success" can build a garbage template.
8. `-maxcoverage 40` maps to `params.maxCov` used by
   `tig->filterContains(params.maxCov, false)` (processTigs.C:132) BEFORE
   consensus - it cap-selects reads per window; it does not cap read count
   fed to edlib inside the stitch loop.
9. There is NO per-tig timeout or iteration cap; canu perl (Consensus.pm:302-
   355) only detects a failed partition AFTER the child process exits.
10. Partitioning (Consensus.pm:114-151) balances partitions by summed
    read-area, an estimate that does not predict stitch difficulty; one bad
    tig stalls its whole partition (observed: they stall, neighbors finish).

## Open questions / unverified

- Where exactly the 8th+ threads are spawned by the utgcns binary (suspected
  readBuffer/sqStore loader helpers); verified only that consensus compute is
  single-threaded and that 7 of 8 threads sleep during a stall.
- Whether an input exists that makes the backtracking genuinely infinite
  (livelock) vs "merely" exponential; both match the observed 38h zero-progress.

## Change log

- 2026-11-03: initial digest at 86d7d3b88, after the trizona HiCanu run
  stalled 38h in two partitions (see pitfalls.md, "case study").
