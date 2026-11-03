# utgcns pitfalls, the known stall, and workarounds (commit 86d7d3b88)

## The stall (verified case 2026-11-03, S. trizona HiFi, canu from env 2.3)

Symptom signature (all verified by live /proc + logs):
- One utgcns per stuck partition: exactly 1 of 8 threads R at ~107% CPU,
  the other 7 S; /proc/PID/io read_bytes/write_bytes frozen; RSS ~17 GB
  fluctuating +-1 GB; no CPU-time anomaly (ticks accrue normally).
- consensus.NNNN.out progress table stops mid-line (unflushed second half);
  nothing else logged since. Diagnose the tig by its half-line tigID.
- Sibling partitions with LARGER estimated work (Consensus.pm `Area`)
  finished in 25-90 min; stuck partitions had ~half the Area. Area does not
  predict stitch difficulty.

Mechanism (from code, control-flow.md): generateTemplateStitch is a
goto-backtracking search (retryCandidate clears the badToAdd memo and loosens
allowLargerShift; retryTig rebuilds the template at errorRate*1.25). On a unitig
whose read placements are mutually inconsistent - heterozygous SV unitigs,
divergent haplotypes, expanded repeats - reads keep failing to extend, and the
search re-runs a growing space of (read order, template trim, error rate)
combinations, each step a full edlibAlign. There is no timeout, no iteration
cap, no progress guarantee: runtime is unbounded in practice ("glitch",
indefinite). Single-threaded by design, so -threads does not help.

Trigger data type: the SAME reads that break hifiasm's cleanup (het SVs,
GC-depressed coverage) are exactly what makes canu's unitig layout
self-inconsistent -> stitch backtracking.

## What NOT to conclude / do

- Do not wait it out; there is no evidence it converges (case above: 38 h
  with zero tig progress on ordinary-size tigs 6469 [897 reads, 0.94 Mb] and
  5109 [710 reads, 0.63 Mb]).
- Do not kill just the utgcns and walk away: Consensus.pm relaunches the
  partition. Expect it: after killing 3881814/3881919 the driver spawned fresh
  utgcns within seconds. To stop, kill the canu perl driver tree by explicit
  PIDs as well (children first). Outputs of the 20 good partitions persist
  (consensus.sh exits 0 if jobid.cns exists) - rerunning canu resumes.
- gdb attach to see the stuck frame needs ptrace_scope relaxed:
  `sudo sysctl -w kernel.yama.ptrace_scope=0` (Ubuntu 26.04 default 1 blocks
  same-uid attach). Do this BEFORE killing if a live backtrace is wanted.

## Workarounds (ordered)

1. Skip the offending tig(s): identify stalled tigIDs from half-lines, then
   rerun canu (it resumes missing partitions only) with higher -e/-em, or use
   utgcns's import/skip options on a repackaged partition. A tig that still
   fails after `canuIterationMax` attempts is reported failed and SKIPPED,
   not fatal ("tig N failed" -> WARNING) - production continues without it.
2 Loosen: bandErrRate k grows with -e; but note -e also multiplies edlib cost.
   Small -e (0.02) + `corMaxEvidenceCoverage`/`-maxcoverage` cap 20 for repeats
   are the community-reported mitigations (Biostars/SEQanswers; canu issues
   #324/#543 family).
3. Give up on canu consensus for the bad unitigs: use the layout (utgcns can
   dump layouts) + racon external, or accept canu contigs built from the good
   partitions' cns and treat the failed tigs with a different tool.
4. For diploid/het genomes request primary-only or use hifiasm-style phasing
   BEFORE consensus; HiCanu consensus on unpurged het data is the failure
   regime (observed here; matches hifiasm-side findings that het SV unitigs
   are the hard structures in this genome).

## Parameter map (utgcns CLI)

| Flag | Meaning | Code point |
|---|---|---|
| -e / -em | errorRate / max errorRate; k=olapLen*bandErrRate, bandErrRate ends at _errorRate | utgcns.C args; unitigConsensus.C:564,702 |
| -maxcoverage | reads per window kept by filterContains BEFORE stitching | processTigs.C:132 |
| -threads | not used in consensus compute | utgcns.C:183 |
| -pbdagcon -edlib | the (only maintained) algorithm+aligner pair | generate() dispatch |
| -min-coverage | minimum column coverage | |

## Verified incident log (therioserver, /home/toki/gof/asm2)

- Launched 2026-11-01 13:58 (asm2_hicanu.sh, canu -p tri -d hicanu
  genomeSize=2.8g useGrid=false maxThreads=56 -pacbio-hifi ERR16751478).
- 22 partitions; 20 done by Nov 2 ~04:00; 20/21 stuck on tigs 6469/5109.
- CPU-advance /60s test separated "slow" from "stuck" (ticks +60s, io static).
- Killed 2026-11-03 17:17: driver relaunched immediately (as Consensus.pm
  dictates); driver tree killed by explicit PIDs; waiter self-guarded exit
  ("no HiCanu contigs"). 20 .cns preserved across kill.
