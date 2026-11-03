# canu utgcns architecture (commit 86d7d3b88)

## Pipeline position

Canu stages: read loading (seqStore) -> overlap (ovlStore) -> 4-unitigger
(ctgStore) -> **5-consensus (this)** -> outputs `tri.contigs.fasta`.
The perl driver (canu/Consensus.pm) partitions tigs, then runs per-partition
`utgcns` via generated `consensus.sh` scripts. With `useGrid=false` each
`sh -c ./consensus.sh N` is a direct child of the canu perl (ps --forest shows
the tree; a defunct sh there = an already-exited job).

## Files (src/utgcns/)

| File | Role |
|---|---|
| utgcns.C | CLI parsing (-P, -R, -T, -O, -e, -em, -maxcoverage, -threads, -pbdagcon, -edlib), calls processTigs |
| utgcns-parameters.C | cnsParameters: loads partition reads (sqStore blob format, magic `sqFile__`), skipTig() prints progress-row first half |
| utgcns-processTigs.C | the sequential per-tig loop; prints progress-row second half (line 135) |
| unitigConsensus.C | the algorithm: generate() dispatch, generateTemplateStitch (386), generatePBDAG (1098), PDAG voting |
| unitigPartition.C / layoutToPackage.C | partitioning support (utgcns -partition mode) |
| libpbutgcns/ | bundled pbdagcon-derived multialignment (abMultiAlign) |

## Job/rerun protocol (verified against server run + Consensus.pm)

- consensus.sh exits 0 immediately if `ctgcns/NNNN.cns` already exists
  (resumable: killing and rerunning canu only redoes missing partitions).
- canu detects failure by ABSENCE of the .cns file after the child exits
  (Consensus.pm:302-327), retries the partition up to `canuIterationMax`
  (Consensus.pm:330-355: "Consensus jobs failed, tried N times, giving up").
- A killed/crashed child without .cns => driver relaunches it - killing a
  single utgcns does NOT stop the run; kill the canu perl to stop retries
  (outputs already written are preserved; the stage resumes from them).
- Watcher-style wrappers elsewhere must check for the FINAL fasta
  (asm2_after_hicanu.sh pattern: exit if tri.contigs.fasta missing).

## Partition blob format (what you can read without running anything)

`tri.ctgStore/partition.NNNN` = sqStore blob, one per partition, ~equal bytes
(all partition files similar size; size does NOT predict difficulty). The
companion `partitioning` **text file** is the key diagnostic: per-tig table
`Tig Reads Length Area MemoryGB Partition`, sortable; `Area` = sum of
read-length*depth work used by Consensus.pm to balance partitions.
A stalled partition's stuck tig is found by matching the half-printed last
line of consensus.NNNNNN.out against the tig ID column (see SKILL.md fact 2).

## Threading

`-threads` is parsed (utgcns.C:183-184, setNumThreads) but the consensus
compute in unitigConsensus.C contains no thread primitives; its comment at
line 1156 ("Construct the graph from the alignments. This is not thread
safe.") confirms single-threaded construction. Extra threads seen in ps
during a stall (7 of 8 sleeping) come from I/O-side helpers only.
