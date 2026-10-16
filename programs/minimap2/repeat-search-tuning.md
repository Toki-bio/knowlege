---
title: "minimap2: Sensitive repeat / SINE search with minimap2"
type: reference
generated: true
---

> Part of the [minimap2](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/minimap2/references/repeat-search-tuning.md)

{% raw %}
# Sensitive repeat / SINE search with minimap2

> Type: reference for the code facts (every claim carries `file:line` at commit
> `3c28777`); the probability tables and the verdict are **notes** (my own
> model and simulation, 2026-10-11), marked as such.

Problem setting: SINE consensus sequences (~100-500 bp) as queries against a
whole genome (~3 Gbp) where the family has ~10^5 copies at 10-30 % divergence,
and minimap2 "misses SINEs easily".

## 1. Where a SINE copy is lost, in pipeline order

| Stage | Mechanism | Default effect on a 300 bp SINE with 10^5 copies | Code |
|---|---|---|---|
| Sketch | only ~2/(w+1) of k-mers become minimizers | w=10 keeps ~18 % of the query's k-mers | `sketch.c:117-138` |
| Seed lookup | query minimizer with more than `mid_occ` reference hits is **dropped**; only `round(run_bp / 500)` lowest-occurrence seeds per high-occurrence run are rescued, each <= 4095 occurrences | every family-wide seed is high-occurrence -> all dropped except **one** rescued seed for the whole 300 bp query (300/500 rounds to 1), and only if it has <= 4095 hits | `seed.c:106-111, 56-96` (`73`: `max_high_occ = (pe-ps)/dist + .499`), `options.c:37-38` |
| Chain | chain needs `>= -n 3` anchors **and** score `>= -m 40` | with one rescued seed no chain can form | `lchain.c:52,67` |
| Primary/secondary | all hits of one query overlap on the query -> one primary, the rest secondary; keep `<= -N 5` secondaries with score `>= 0.8 x primary` | at most 6 copies reported per query even if seeding were perfect | `hit.c:165, 265-267`, `map.c:210` |
| Alignment | only with `-c`/`-a`; filters `mlen < -m`, `dp_max < -s` (80 for map-ont) | 300 bp at 72 % identity scores ~ 300*(0.72*2 - 0.28*4) = 96 before gaps, so `-s 80` starts biting below ~72 % | `map.c:217`, `hit.c:307-311`, `options.c:52` |
| Output | secondaries printed unless `--secondary=no` | PAF `tp:A:S` lines | `map.c:612-613` |

The `rl:i` tag is the diagnostic: if `rl` is close to the query length, the
seeds were thrown away by the occurrence filter, not by chaining
(`seed.c:117-122`, `format.c:447`).

## 2. Which minimizers survive the occurrence filter

`mid_occ` is derived once per index part (`options.c:73-79`):

```
mid_occ = count at the (1 - f) quantile over DISTINCT minimizers + 1   (index.c:217; f = -f = 2e-4)
mid_occ = max(min_mid_occ=10, min(max_mid_occ=1e6, mid_occ))            (options.c:75-78)
```

Printed with `-v3` as `mid_occ = N` (`options.c:81-82`). For a mammalian genome
at k15/w10 this is typically a few hundred (unverified here; measure it). A
consensus k-mer is retained in roughly `C * p^k` copies of a family of `C`
copies at identity `p` (plus `~ 2 * 3e9 / 4^k` random background hits):

| k | p | copies retaining the k-mer (C = 1e5) | background per k-mer, 3 Gbp | exceeds a mid_occ of ~500? |
|---|---|---|---|---|
| 15 | 0.75 | 1,300 | ~6 | yes -> dropped |
| 15 | 0.85 | 8,700 | ~6 | yes |
| 15 | 0.95 | 46,000 | ~6 | yes |
| 11 | 0.75 | 4,200 | ~1,400 | yes (but see k=11 note) |
| 11 | 0.85 | 16,700 | ~1,400 | yes |

So with defaults **every minimizer that is diagnostic of the family is
dropped**, and the better conserved the k-mer, the more certainly it is
dropped. The one rescued seed per 500 bp (`seed.c:73`) must also have
<= `max_max_occ` = 4095 hits (`seed.c:89-91`, `options.c:37`), which a
conserved SINE k-mer violates. Net: a high-copy SINE consensus produces zero
or one anchor per copy and no chain. This is by design: README
`README.md:325-328` describes exactly this filter.

To keep the seeds: `-f INT` with INT >= 1 sets `mid_occ` directly and skips the
clamp (`main.c:330`), e.g. `-f 1000000`; or `-f 0` which ends up at
`max_mid_occ` = 1,000,000 via the clamp (`index.c:204`, `options.c:77-78`).
For `mid_occ` beyond 1e6 use `-f 10000000` or `-U 10,100000000`
(`main.c:322-324`). `-e` (occ_dist) and `max_max_occ` become irrelevant once
`mid_occ` exceeds the family's copy count because `n > max_occ` is then false
(`seed.c:65, 109-111`). At k=11 on a 3 Gbp genome the mean 11-mer already
occurs ~1,400 times, so the derived `mid_occ` at `-f 2e-4` will be large
(unverified value; check `-v3`), but the clamp to 1e6 and the per-run rescue
rules still apply.

## 3. Chaining math: why 300 bp at 75-85 % identity fails `-n 3` / `-m 40`

Chain score for colinear anchors without indels equals the number of query
bases covered by the anchors' k-mer spans (`comput_sc`: `sc = min(q_span,
dg)`, `lchain.c:124-126`; first anchor contributes `q_span`, `lchain.c:171`).
Indels subtract `0.8*0.01*k*|dd| + 0.5*log2(|dd|+1)` per diagonal shift
(`lchain.c:129-135`, `map.c:273`): about 1 point per 1-3 bp indel at k=15.
Chains are emitted only if `score >= -m` and `anchors >= -n` (`lchain.c:52,67`).

**Notes (model + simulation, not minimap2 code).** Independent substitutions at
rate `1-p`; expected conserved k-mers in a 300 bp query = `(301-k) p^k`;
minimizer density `~2/(w+1)` (`w=1` -> every k-mer, `sketch.c:117-138`).
Simulation: 400 random 300 bp consensus/copy pairs per cell, hashed canonical
minimizers as in `sketch.c`, counting shared minimizers at homologous
positions and the bases they cover (= chain score with no indels):

| k | w | identity | E[shared minimizers] | P(pass `-n3 -m40`) | P(pass `-n2 -m25`) | P(pass `-n1 -m20`) | E[chain score] |
|---|---|---|---|---|---|---|---|
| 15 | 10 | 0.75 | 0.6 | 0.01 | 0.06 | 0.09 | 8 |
| 15 | 10 | 0.80 | 1.5 | 0.11 | 0.29 | 0.32 | 18 |
| 15 | 10 | 0.85 | 3.8 | 0.49 | 0.74 | 0.78 | 41 |
| 15 | 10 | 0.90 | 8.9 | 0.92 | 0.98 | 0.99 | 85 |
| 15 | 1 | 0.75 | 3.7 | 0.09 | 0.27 | 0.35 | 17 |
| 15 | 1 | 0.85 | 24.5 | 0.90 | 0.95 | 0.97 | 77 |
| 13 | 3 | 0.75 | 3.3 | 0.15 | 0.42 | 0.45 | 22 |
| 13 | 3 | 0.80 | 7.1 | 0.54 | 0.83 | 0.84 | 44 |
| 13 | 1 | 0.75 | 6.5 | 0.23 | 0.53 | 0.56 | 27 |
| 13 | 1 | 0.80 | 15.8 | 0.70 | 0.89 | 0.90 | 54 |
| 11 | 10 | 0.75 | 1.7 | 0.06 | 0.21 | 0.41 | 16 |
| 11 | 10 | 0.85 | 7.2 | 0.83 | 0.93 | 0.97 | 60 |
| 11 | 3 | 0.75 | 5.6 | 0.36 | 0.65 | 0.76 | 33 |
| 11 | 3 | 0.80 | 12.3 | 0.82 | 0.96 | 0.97 | 62 |
| 11 | 3 | 0.85 | 21.6 | 0.98 | 1.00 | 1.00 | 99 |
| 11 | 1 | 0.75 | 12.4 | 0.55 | 0.78 | 0.84 | 43 |
| 11 | 1 | 0.80 | 25.2 | 0.91 | 0.99 | 0.99 | 77 |
| 11 | 1 | 0.85 | 49.2 | 1.00 | 1.00 | 1.00 | 123 |

Analytic expectations agree (k15 w10 at 85 %: `286 * 0.85^15 * 2/11 = 4.5`;
k11 w1 at 75 %: `290 * 0.75^11 = 12.2`). Reading: with the default sketch a
copy at 80 % identity has ~1-2 shared minimizers and passes `-n 3 -m 40` about
10 % of the time *even if nothing were filtered*; at 75 % essentially never.
Getting to >= 80 % per-copy sensitivity at 75 % identity needs `k <= 11` **and**
`w = 1` **and** `-n 1 -m 20`; indels (common in SINE copies) lower every row.

## 4. Alignment stage at high divergence

- Scoring break-even identity is `b / (a + b)`: default `-A2 -B4` -> 67 %;
  `asm20` (`-A1 -B4`) -> **80 %**, so `asm20` is worse than default for 20-30 %
  divergence despite its name (`options.c:147-148`). `-A2 -B2` -> 50 %;
  `-A1 -B1` plus small `-O/-E` is the "BLASTN-like" end of the range, subject
  to `(O1+E1)+(O2+E2) <= 127` (`options.c:251-255`).
- Z-drop `-z 400` barely matters for 300 bp (max attainable score 600,
  `ksw2.h:171-182`); `-z 100` or lower would start truncating.
- `-s` filters on `dp_max` = best-scoring segment of the final CIGAR
  (`hit.c:310`; `align.c:254`): default 80 (`options.c:52`). At `-A2 -B4` a 300
  bp copy at 75 % identity scores ~150 minus gaps, at 70 % ~60: set `-s 20-40`.
- `-m` is reused after alignment as a floor on `mlen` (exact matched bases
  from the CIGAR, `hit.c:309`); `-n` as a floor on anchor count (`hit.c:307`).
- DP windows are small for short queries: extension is bounded by the unaligned
  query length plus `(l*a - q)/e` reference bases (`align.c:723-733, 749-758`),
  so cost per hit is ~`300 * (300 + 2*1.5*bw)` cells; the band `1.5*bw+1` with
  `-r 500` is wider than the query (`align.c:660`), so `-r 50` is enough for a
  300 bp query and makes every DP cheaper and chaining stricter (`lchain.c:122`).
- `--end-bonus` default -1 (`options.c:51`) means extension never prefers
  reaching the query end (`ksw2_extz2_sse.c:304-306`); `sr` uses 10. For
  consensus-vs-copy full-length alignments `--end-bonus 10` reduces soft clips.
- Hits per query: without `-P`, `1 + min(-N, #passing -p)` chains are aligned
  and printed (secondaries as `tp:A:S`, SAM flag 0x100, `format.c:556`). `-p 0`
  turns `mm_select_sub` off entirely (`hit.c:257`), so `-p 0` alone keeps every
  chain regardless of `-N`; `-P` additionally skips `mm_set_parent` so every
  hit is `tp:A:P` with its own MAPQ (`map.c:208-212, 219-223`).

## 5. A sensitive-repeat parameter set, and what each flag costs

```sh
minimap2 -t 16 -k 11 -w 1 -f 1000000 -U 10,100000000 \
         -n 1 -m 20 -s 30 -r 50 -g 100 -z 200 \
         -A 2 -B 2 -O 4,24 -E 2,1 --end-bonus 10 \
         -P -c -K 1M --secondary=yes ref.fa sine_consensus.fa > hits.paf
```

| Flag | Why | What breaks | Code |
|---|---|---|---|
| `-k 11` (or 13) | 4-7x more conserved k-mers than k15 at 75-85 % | random 11-mer matches: each query 11-mer has ~1,400 random hits in 3 Gbp (both strands), so a 300 bp query yields ~4e5 random anchors **before** the true ones; k=9 gives ~2.3e4 per k-mer, ~7e6 anchors per query: impractical. `k <= 28` only hard limit (`sketch.c:84`) | `index.c:53` (bucket bits capped at 2k), `lchain.c:163-166` (36 B/anchor) |
| `-w 1` | every k-mer indexed; no sampling loss | index size: 3 Gbp -> ~3e9 minimizers -> ~48 GB during build (16 B each, `index.c:306-313`) + ~24 GB final positions (`index.c:248`) + 1.5 GB sequence; `-w 3` halves it. Build is single-pass but `worker_post` sorts every bucket (`index.c:236`) | `parameters-and-presets.md` memory formula |
| `-f 1000000` / `-f 0` | keep family-wide seeds | anchors per query ~ `sum over seeds of occurrence`: 1e5 copies x ~12-50 shared seeds = 1e6-5e6 anchors (`seed.c:124`) -> 40-180 MB per thread per query, chaining O(anchors x predecessors within 100 bp) is fine, but `mm_gen_regs`/`mm_set_parent` then see ~1e5 chains (`hit.c:52, 125`) | `main.c:329-331` |
| `-U 10,100000000` | lift the 1e6 clamp if any k-mer exceeds it (k=9 or satellite-like) | none | `main.c:322-324`, `options.c:77-78` |
| `-n 1 -m 20` | accept single-anchor and short chains | spurious chains: at k=11 w=1 a single anchor scores 11 < 20, so 2 random anchors within `-r`/`-g` of each other on one diagonal are needed; expected ~`4e5 * 4e5 * 100 / 6e9 ~ 10` per query at `-g 100`. `-m 25-30` removes most of them; `-n 2` too | `lchain.c:52,67` |
| `-s 30` | let 70-75 % identity alignments through | more noise alignments pass | `hit.c:310` |
| `-r 50 -g 100` | SINE copies have few long indels; narrow band = fewer spurious links and cheaper DP | copies with a > 50 bp indel split into two hits (`align.c:848-868`) | `lchain.c:118-122`, `align.c:660` |
| `-A2 -B2`, `-O4,24 -E2,1` | break-even 50 % identity | more false extension into flanks; Z-drop still trims | `align.c:659` |
| `-P` (or `-p 0 -N 1000000`) | report every copy | `-P` disables `--secondary=no` (error) and skips `mm_update_dp_max`; output lines = chains | `options.c:266-270`, `align.c:1114` |
| `-c` | needed for CIGAR/identity (`de:f`); without it output is chains only, fast | 1e5 ksw2 calls per query: each ~300 x 450 cells, ~1e10 cells per query; a few seconds to a minute per query per thread (unverified) | `map.c:217` |
| `-K 1M` | all hits of a mini-batch stay in RAM until written (`map.c:578-641`); 1e5 hits x (`mm_reg1_t` ~80 B + CIGAR) per query x queries per batch | smaller batches reduce threading efficiency (`-K` doc `minimap2.1:597-603`) | `options.c:56`, `bseq.c:99-100` |
| `--secondary=yes` | default; only matters if a preset set `NO_PRINT_2ND` (`sr`, `splice:sr`) | | `main.c:282-283` |

Not useful here:
- `--rmq`: RMQ chainer for long gaps/inversions in long queries (`map.c:275-277`,
  `lchain.c:250-368`); no sensitivity gain for 300 bp queries; refuses `--sr`
  (`options.c:209-213`).
- `--sr` / `-x sr`: k21/w11 (`options.c:152`) -> `0.75^21 * 280 = 0.7`
  conserved 21-mers per copy; `MM_F_SR` also fixes `mid_occ = 1000`
  (`options.c:166`), uses longest exact stretch + ungapped-first alignment
  (`align.c:664-669, 823-834`), rechains with `max_occ` 5000 only when no chain
  formed (`map.c:293-316`), hides secondaries (`NO_PRINT_2ND`, `options.c:153`).
  Worst possible choice for diverged repeats. `--sr` alone only sets the flag.
- `-H` (HPC): collapses homopolymers, so A-rich SINE tails shrink and HPC
  k-mers span more bases (`sketch.c:94-104`); it does not change mismatch
  tolerance. Incompatible with `--sr` (`align.c:655`) and `--qstrand`.
- `-2`: only adds a third pipeline thread for I/O (`main.c:201`, `map.c:677`).
- `-T` (SDUST): would remove low-complexity *query* minimizers; marked unused
  in the help (`main.c:378`).

## 6. Throughput and memory, 3 Gbp genome, `-k 11 -w 1` (derived, not measured)

- Index: ~3e9 minimizers. Build peak ~ 48 GB (bucket arrays, `index.c:306-313`)
  + up to 24 GB positions as buckets convert (`index.c:248`) + 1.5 GB sequence
  with `-c`; final ~26 GB. With `-w 3` roughly half; default k15 w10 ~ 5-9 GB.
  `-I` default 8 Gbp (`options.c:11`) keeps it in one part.
- Mapping: `-t` threads each hold anchors + chaining arrays (36 B/anchor,
  `lchain.c:163-166`); a query with 2e6 anchors needs ~70 MB; arena resets
  above 500 MB (`map.c:371-376`). Random background anchors alone are ~4e5 per
  300 bp query at k=11 (1.6e3 at k=15, 2.6e4 at k=13, 6.7e6 at k=9).
- Pipeline: reader / `kt_for` mapper / writer (`map.c:541-643`, `kthread.c`);
  with ~1e5 hits per query the writer and memory per batch, not CPU, dominate,
  hence small `-K`.

## 7. Verdict (notes)

minimap2 *can* be bent into a repeat finder for copies at >= ~85 % identity:
`-k 11..13 -w 1..3 -f 1000000 -n 2 -m 25 -s 30 -r 50 -p 0` (or `-P`) with
`-A2 -B2` recovers most copies in the simulation table and no hard-coded
limit stops it. Three things make it the wrong tool for *high-copy, 20-30 %
diverged* elements rather than merely a mis-tuned one:

1. **The occurrence filter is the opposite of the task.** `-f`/`mid_occ`,
   the one-seed-per-500 bp rescue and the 4095 cap (`seed.c:56-111`,
   `options.c:37-38`) exist to discard exactly the seeds a 10^5-copy family
   produces, and MAPQ is penalised by `rep_len` (`hit.c:445`). Disabling it
   means indexing and chaining millions of anchors per query.
2. **Minimizer sampling plus exact k-mer seeds has no mismatch tolerance.**
   At 75 % identity even `k=11, w=1` leaves ~45 % of copies without a
   qualifying chain (table above); every `w > 1` drops a further `1 - 2/(w+1)`
   of the few conserved k-mers at random. There are no spaced seeds, no
   mismatch-tolerant seeds, no profile/HMM scoring; `-k 9` is the only lever
   left and it generates ~7e6 random anchors per query on 3 Gbp.
3. **Cost scales with copy number, not with sensitivity.** Every output copy is
   a chain plus a ksw2 alignment, held in memory per `-K` batch; the default
   `-N 5`/`-p 0.8` exists because the design assumes one true locus.

For the stated problem (SINE consensus vs genome, 10-30 % divergence, 10^5
copies) a k-mer-exhaustive or profile-based search (RepeatMasker/nhmmer,
BLASTN/megablast with word 7-11, LAST, or a dedicated seed-and-extend with
mismatch-tolerant seeds) matches the problem's structure; minimap2 is a
mapper whose defaults encode the assumption that repeats are noise. Use it
with the settings above only as a fast first pass for the well-conserved
(>= 85 %) fraction, and read `rl:i` to see what it threw away.
{% endraw %}
