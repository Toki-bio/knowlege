# Is `bwa mem` usable as a sensitive repeat finder? (Type: notes)

> Interpretation built on the verified facts in `seeding-chaining.md` and `extension-scoring-output.md`. Timing figures are order-of-magnitude estimates, not measurements.

## The user's symptoms, traced to code

| Symptom | Mechanism | Where |
|---|---|---|
| "misses tons of hits" | (1) SMEM seeding reports only the *longest* exact match at each query position; high-copy shorter k-mers are hidden behind any longer match to a few copies | `bwt.c:289-351`, `bwamem.c:147-157` |
| | (2) Re-seeding fires once and only for SMEMs with <= 10 occurrences (`-s`) and length >= 1.5k (`-r`); the third pass only emits seeds with < 20 occurrences (`-y`) | `bwamem.c:163-164,176`, `bwt.c:370` |
| | (3) `-c 500` subsamples each seed to 500 copies | `bwamem.c:304-305` |
| | (4) No exact 19-mer exists in ~74% of 75%-identity copies (and ~44% at 80%) | math in `seeding-chaining.md` section 5 |
| | (5) Default `-A1 -B4` gives negative expected score at < 80% identity, so extensions die by Z-drop and fall under `-T 30` | `bwamem.c:57`, `ksw.c:494-500`, `bwamem.c:1048` |
| | (6) Everything overlapping the best hit is secondary and hidden without `-a`; XA is all-or-nothing above `-h 5` | `bwamem.c:528-543,1049`, `bwamem_extra.c:144` |
| "enormous index" | 1.75 B/base, fixed; ~5.4 GB RAM for human; it is the price of FM-index rank+locate, not a misconfiguration | `index-and-memory.md` |

Items 1-3 are the ones a user cannot see: they remove candidate copies before any alignment exists, so `-a`/`-T`/`-D` cannot bring them back.

## Closest configuration

Goal: make every pass emit high-copy seeds, stop subsampling, score for divergence, output everything.

```
bwa mem -t T \
  -k 12 -r 1 -s 1000000 -y 1000000 -c 1000000 \     # seeding: short seeds, re-seed everything, LAST-like tiling with no occurrence cap, keep all occurrences
  -D 0 -W 0 -N 1000000000 \                          # chaining: never drop a chain for weight; no chain count cap (default already)
  -A 1 -B 1 -O 1 -E 1 -L 0 -d 50 -w 100 \            # scoring for ~20-30% divergence (same as -x pacbio/ont2d), free clipping
  -T 25 -a -h 1000000 -z 0.1 -Y \                    # output: lower threshold, all secondaries, XA never truncated (irrelevant with -a), soft clips
  genome.fa sine_consensus.fa > hits.sam
```
Why each piece (all verified lines in the other files):
- `-k 12`: exact 12-mers exist in ~90% of 75%-identity copies over 300 bp; random background on a 3 Gbp genome ~370 occurrences/12-mer, so true-copy chains (2+ seeds, weight >= 24) still outrank most noise chains (weight 12), and `-T 25` removes noise after extension. `-k 11` doubles the noise per seed for little gain.
- `-r 1 -s 1000000`: re-seed every pass-1 SMEM of length >= k from its midpoint, regardless of occurrence (`bwamem.c:163`), recovering one level of shorter/more-frequent matches.
- `-y 1000000`: pass 3 then emits a `(k+1)`-mer at every tiling position with any occurrence count (`bwt.c:370`) — the only systematic source of high-copy seeds.
- `-c 1000000`: `step = 1`, all occurrences located (`bwamem.c:304-305`); `frac_rep` still -> 1 so MAPQ is 0 everywhere (ignore MAPQ, use `AS`).
- `-D 0`: disables weight-based chain dropping (`384`) **and** the secondary output cut at 50% of the primary (`1050`). `-W 0` is already the default.
- `-A1 -B1 -O1 -E1 -L0`: `+0.5/base` expected at 75% identity instead of `-0.25`; `-d 50` scaled to A=1 (default 100 is fine too). Keep `-w 100` unless copies have long indels.
- `-a`: print secondaries (flag 256, MAPQ 0, SEQ `*`, but CIGAR/`AS`/`NM`/`MD` present). `-T 25` ~ 50 aligned bases at 75% identity.
- `-Y`: supplementary pieces soft-clipped so their SEQ is complete.

Post-process: treat every SAM line (flags 0, 16, 256, 272, 2048, 2064) as a hit; filter on `AS`, alignment length, and `NM`; de-duplicate by coordinates.

## What that costs (`notes`, order of magnitude)

Per 300 bp query against a 3.1 Gbp genome with a 10,000-copy family:
- Seeding: ~25 tiled seeds x 10,000 occurrences x `bwt_sa` (~16 random LF steps each, `bwt.c:86-96`) ~ 4M cache-missing reads -> roughly 0.5-2 s.
- Chaining: 250k seed insertions into a B-tree, ~10k chains; `mem_chain_flt` is O(chains x kept) ~ 1e8 cheap comparisons.
- Extension: ~10k chains x 2 x `ksw_extend2` on ~300 x (band 201) cells ~ 1e9 scalar cells -> seconds.
- Primary marking `mem_mark_primary_se_core` is O(n^2) in alignments (~1e8 for 10k).
Expect seconds to a minute per high-copy consensus, i.e. fine for a library of hundreds of consensus sequences, not for millions of reads. Memory: index 5.4 GB + per-query region arrays (10k x ~100 B) — negligible. Random hits: with `-c` uncapped, a 12-mer hitting 370 random places costs 370 locates + 370 extensions each; `-T` removes them from output but not from run time.

Reducing SA interval to 8 via `bwt2sa -i 8` (`index-and-memory.md`) cuts locate cost ~4x at +1.5 GB of `.sa` for human.

## What still cannot be fixed by flags

- Copies with no exact `(k+1)`-mer in common with the consensus are never found (seed-and-extend is the architecture). At 70% identity and k=12 that is ~10-15% of copies; at 65% the majority.
- Pass 1 still hides high-copy k-mers behind longer rarer matches; only passes 2 (one level) and 3 (tiling) compensate, so coverage of copies is uneven across the query.
- `frac_rep`, `sub`, `sub_n` make MAPQ meaningless; `XS` is suppressed on secondaries.
- Chaining requires seeds within `-w` of the same diagonal; copies with large indels relative to the consensus fragment into separate chains/alignments (`bwamem.c:228`), later possibly merged by `mem_patch_reg` only when nearly co-linear (`447-449`).
- Output volume: `-a` with 10k copies prints 10k lines per query; fine for a few hundred queries.

## Honest bottom line

`bwa mem` is a read mapper: its seeding, occurrence cap, chain filter, secondary/MAPQ model and default scoring all exist to find *the* placement of a near-identical read quickly and to suppress repeats. It can be bent into an exhaustive seed-and-extend aligner with the flag set above, and for 10-20% diverged copies it will then recover most of them at seconds per query, but it is not a sensitive repeat annotator: it has no profile/HMM scoring, no sensitivity guarantee, and it hides high-copy seeds by design. `bwasw` and `bwa aln` are worse for this (see `bwasw-aln-fastmap.md`). Use `fastmap -l 12 -w 1000000` first to see what exact anchors exist; if a copy has none, bwa cannot find it. For the actual task (SINE consensus vs genome at 10-30% divergence), tools built for divergent repeat search (RepeatMasker/nhmmer-style profile search, BLASTN/LAST-style ungapped-seed + spaced-seed local aligners, or minimap2's `asm20`-class presets) are the appropriate class; `bwa mem` is at best a fast first pass for the closest copies. (This last sentence is an opinion, not derived from the bwa code.)
