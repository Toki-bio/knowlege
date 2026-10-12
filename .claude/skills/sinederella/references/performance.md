# Performance: where the wall-clock goes (Q4)

Line numbers at `90d8ac4`. All timings are the author's, quoted from docs
(servers: therioserver, KIT; not reproduced here).

## Measured run times

| Genome / bank | Total | Breakdown | Source |
|---|---|---|---|
| hg38 chr21, 46.7 Mb, 2 Alu consensuses | 20 min full | `--add` AluY 10 min; `--exclude` 2 min; `--resume` 1 min; publish 2.5 min | `docs/RELEASE_CHECK_2026-10-05.md:22-26` |
| mm39 chr19, 61.4 Mb, B1 + B2 | 15 min | | `:55` |
| rsi *Rhinolophus sinicus*, 2.1 Gb, 21 consensuses, 66,569 loci | **84 min** | step 1 incl. satellite stage **61 min**, step 2 **12 min**, step 3 1 min, step 4 4 min, length versions 1 min, audit 3 min | `:75` |
| tbr *Taphozous*, 2.3 Gb, 2 consensuses (VES 620k copies) | 54 min + 2 min publish; `--add` 4 MEG consensuses 9 min | | `:101,116` |
| rle *Rousettus*, 2 consensuses | 4 min; `--add` 8 min | | `:101,116` |
| single `sear` query vs 2.3 Gb tbr | **51-76 s per consensus** (unmasked 69 s Rhin-1, 64 s MEG-RL) | THREADS=12 | `docs/MASKING.md:121-132` |
| planted 8 Mb, 120 bp family with 1,500 copies | 37-46 s per query; ~4 s with the copies N-masked; ~3 s with them deleted | `-T 1` | `docs/MASKING.md:19-24` |
| satellite stage, rsi, 2,739 kind-B runs | **13 h 24 min** serial all-vs-all with `-z 11` -> **9.5 min** lag-1/2 pairs, `-z -1`, 32 processes; whole stage 718 s; in a run 11 min 33 s | | `docs/SATELLITES.md:282-324` |
| flankscan on rsi 67k copies / tbr 622k copies | 3.5 min / 18 min (TRF 8.5, fs3 ssearch 7.7) | | `flankscan/HANDOFF.md:163` |
| genome-wide TRF (ground truth, cobra 1.8 Gb) | ~9 min on 44 cores; windowed TRF two orders of magnitude less | | `docs/SATELLITES.md:175,192` |
| SubFam final alignment, 595 Timema consensuses | > 9 h with N columns as sequence (1.2.2) -> 2 s (1.2.4) | | `SubFam:29`, `tools/vendor/SubFam.sh:501` |
| publish alignments | "hours of MAFFT on a big run" | | `SINEderella:438`, audit D26 |

## Where the time goes, by mechanism

1. **Step 1 ssearch36 = bank size x genome.** One `sear` call per consensus,
   serial (`step1:109-123`); each is 20 serial ssearch36 runs with `-T THREADS`
   (`sear:389-465`). ~1 min per consensus per 2 Gb at 12 threads
   (`docs/MASKING.md:126-130`) -> 21 consensuses ~20-25 min of the rsi step 1;
   the rest of its 61 min is the satellite stage (11.5 min after the fix) and
   SubFam/MAFFT. `notes`: the author's own measurement says the Smith-Waterman
   scan "costs the same" whatever the query, and masking gains only when the
   query would align to many masked copies (`docs/MASKING.md:131-132`);
   conversely on the planted genome 2 % masked bases cut the time tenfold
   because "the time is the known copies being aligned": windows that score
   above the 8-bit limit are redone in 16 bits and then re-aligned for display
   (`docs/MASKING.md:44-46`). So cost = scan term + (number of copies) x
   (alignment term); hit-rich families (VES 620k) are alignment-bound.
2. **Step 2 = 10 x consensuses x copies.** `ssearch36` with all consensuses as
   query against each 20,000-copy part, 10 cycles (`step2:181-189`): 12 min for
   21 x 66,569 on rsi. `--add` avoids ~75 % of it by carrying firm copies
   forward (`MANUAL.md:195`, `step2:546-593`).
3. **Satellite kind-B check** was the pathological case (13 h) until it aligned
   only lag-1/2 pairs with `-z -1` in parallel single-threaded processes
   (`tools/satellite_kindB_verify.py:129-137`, `docs/SATELLITES.md:282-306`).
   A k-mer pre-screen was tried and rejected: pairs at 70-85 % identity share
   no 12-20-mers (`docs/SATELLITES.md:338-344`).
4. **Genome I/O:** the satellite tools stream the whole genome once per
   consensus (audit D10, `docs/AUDIT_2026-10-05.md:97-100`); `mask_genome.py`
   and `sanitize_fasta` are single passes.
5. **MAFFT L-INS-i** on 100-row plates with long flanks is the publish cost;
   threaded refinement is 3-7x faster but changes the alignment, so plates use
   1 thread each and run subfamilies in parallel (`step8a:53-57`,
   `tests/bench/run_bench.sh`).

## Parallelism knobs

| Knob | Reaches | Code |
|---|---|---|
| `THREADS` env | `sear -T`/minimap2 `-t`, step 2/3/4 `-T`, satellite `--threads` (process count), length variants, fs9 | `SINEderella:108-111`, `step1:140`, `step2:488` |
| `SINEderella_multi --parallel N --threads N` | N genomes at once | `SINEderella_multi:31-32,309` |
| `PAR_SF`, `MAFFT_THREADS` | step 8a plates | `step8a:58-59` |
| `CONSENSUS_AUDIT_JOBS` (8) | consensus audit rebuilds | `SINEderella:413` |
| not honoured: step 1 MAFFT `--thread nproc`, SubFam `nproc` | | `step1:300`, audit D8 |

No consensus-level or part-level parallelism exists in step 1; `sear_multi`
(one ssearch36 call per part with the whole bank as query) is the author's
untested answer to the per-query overhead (`sear_multi:4-7`,
`benchmark_sear.sh`, no results committed).

## Author's statements on speed and sensitivity

- "Keep ssearch36 ... Do not move the search to blastn, and do not use
  `ssearch36 -S` as the mask" (`docs/MASKING.md:5-15`); minimap2 "no
  SINE-scale sensitivity and no case" (`:60`).
- "ssearch36 loses hits when the library is huge and one family dominates";
  fixed in flankscan by 20,000-record chunks and `-Z 20000` (`HANDOFF.md:150-153`,
  `fs3_partners.sh:85-91`). `notes`: `sear`'s parts carry no `-Z`.
- `-z 11` is non-deterministic, ~1 % count noise (`docs/COMPOSITES.md:88`,
  `HANDOFF.md:48-49`); every family count in the release check is "within the
  `-z 11` noise (+-1 %)" (`RELEASE_CHECK:87`).
- Recall is bounded by the 0.8 length rule, not divergence (`RELEASE_CHECK:30`).
- Old copies at 61-68 % identity (rsi MEG-T2) are assigned by votes but cannot
  be placed by flank junctions; nhmmer profile HMMs rescued 9 of 145 "nomain"
  copies and 0 of 109 MEG-T2: "old copies, not a search-engine limit"
  (`HANDOFF.md:190-193`; audit perspectives item 5).
- `step2_asSINEment.sh:6`: "No alignment steps - fast assignment only".
- `--add` exists for efficiency (`MANUAL.md:182-195`, `:1356-1358`).
- The only explicit TODO about the search: use raw ssearch36 hits instead of
  `sear`'s merged hits for satellites (`docs/SATELLITES.md:170,193,206`), done
  via `rawhits.tsv`. `extract_alignments.sh:28` has a "TODO (future)" block
  (historical script).
