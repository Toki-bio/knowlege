---
title: "bwa: Index construction, on-disk size and memory"
type: reference
generated: true
---

> Part of the [bwa](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/bwa/references/index-and-memory.md)

{% raw %}
# Index construction, on-disk size and memory

Answers: "is the bwa index really enormous?" — it is ~1.75 bytes per genome base, linear, with fixed sampling rates baked into `bwa index`.

## What is indexed

- The FM-index text is the forward genome **and** its reverse complement, so `seq_len = N = 2*G` (`bntseq.c:306-312` adds the revcomp when `for_only == 0`; `bwtindex.c:271` calls it that way for the BWT step). This is what lets `bwa mem` search both strands with one bidirectional extension (`bwt_extend`, `bwt.c:262-275`) and never reverse-complements the read.
- The shipped `.pac` is forward-only (`bwtindex.c:305`, `for_only = 1`); reverse-strand coordinates `>= l_pac` are mapped back with `bns_depos` (`bntseq.h:87-90`).
- Every non-ACGT reference base is replaced by a pseudo-random base (srand48(11), `bntseq.c:266,295-296`), and the N-runs are logged in `.amb`. Spurious seeds can land in former N regions; `bns_cnt_ambi` counts them.

## BWT construction algorithms (`-a`)

| `-a` | Code | When auto-chosen | Memory / limits |
|---|---|---|---|
| `is` | `bwt_pac2bwt(.., use_is=1)` -> `is_bwt` (`bwtindex.c:92-93`, `is.c:208-223`) | `l_pac <= 50 Mbp` (`bwtindex.c:275`; note this `l_pac` is already 2G) | 1 byte/base text copy (`bwtindex.c:83`) + `4*(N+1)` bytes int SA (`is.c:211`); `int n` so N must be < 2^31, i.e. G < ~1 Gbp. Man page: "does not work for databases longer than 2GB" (`bwa.1:69-70`). For human it would need ~31 GB and overflows. |
| `bwtsw` | `bwt_bwtgen2(pac, bwt, block_size)` (`bwtindex.c:281`, `bwt_gen.c:1592-1599`), BWT-SW incremental construction, `-b` block (default 10,000,000) | `l_pac > 50 Mbp` | man page: ~5 GB RAM and ~3 h for human (`bwa.1:758,767`); "does not work for short genomes" (usage `bwtindex.c:241`). |
| `rb2` | `bwt_pac2bwt(.., use_is=0)`: rope/RLE incremental insertion (`bwtindex.c:95-119`) | never auto | not measured here (unverified) |

After construction, every path does `bwt_bwtupdate_core` (interleave Occ, `bwtindex.c:150-172`) and `bwt_cal_sa(bwt, 32)` (`bwtindex.c:316`), the latter a single sequential pass of N LF-mapping steps (`bwt.c:77-81`), each a random memory access into the BWT. `notes:` for N = 6.2e9 this pass alone is tens of minutes.

## File sizes (exact formulas, from the writers)

Let G = forward length, N = 2G.

| File | Bytes | Where |
|---|---|---|
| `.pac` | `G/4 + 1` (+1 pad if `G%4==0`) | `bntseq.c:316-323` |
| `.bwt` | `40 + 4*bwt_size`, `bwt_size = ceil(N/16) + 8*(ceil(N/128)+1)` words -> `N/4` text + `N/4` Occ = `N/2` = **G bytes** | `bwtindex.c:74` (text words), `155-156` (Occ words: `n_occ * sizeof(bwtint_t)` = 8 words = 32 B per 128 bases), `bwt.c:389-391` |
| `.sa` | `56 + 8*((N+32)/32 - 1)` -> `N/4` = **G/2 bytes** | `bwt.c:73` (`n_sa`), `bwt.c:400-404` |
| `.ann`, `.amb` | text, proportional to #sequences and #N-runs | `bntseq.c:65-95` |

Per forward base: `.bwt` 1.0 B + `.sa` 0.5 B + `.pac` 0.25 B = **1.75 B/base**.

| Genome | `.bwt` | `.sa` | `.pac` | Total |
|---|---|---|---|---|
| 3.1 Gbp (human) | 3.1 GB | 1.55 GB | 0.78 GB | ~5.4 GB |
| 1 Gbp | 1.0 GB | 0.5 GB | 0.25 GB | 1.75 GB |
| 100 Mbp | 100 MB | 50 MB | 25 MB | 175 MB |

The man page's "sampe uses ~5.4GB" (`bwa.1:761-763`) matches this; `bwa mem` loads all three (`bwa.c:299-312`), so its resident set is ~5.4 GB for human plus batch buffers. `fastmap` and `bwasw` skip `.pac` at load (`fastmap.c:444`, `bwtsw2_main.c:83`): ~4.65 GB. `bwa aln` loads only `.bwt`: ~3.1 GB (`bwtaln.c:171-175`; man page says 3.2 GB, `bwa.1:761`).

So "enormous" is relative: ~1.75x the genome in bytes, about the same order as other FM-index mappers; it is not tunable from `bwa index`. For comparison, 2-bit `.pac` alone is 0.25 B/base — the BWT text (0.25) + Occ (0.25) + SA (0.5) are the cost of O(1) rank and O(32) locate.

## Knobs that exist (and the ones that don't)

- **SA sampling** (`sa_intv`): hard-coded 32 in `bwa index` (`bwtindex.c:316`). Change only via the manual pipeline `bwa fa2pac; bwa pac2bwt(gen); bwa bwtupdate; bwa bwt2sa -i INT` (`bwtindex.c:188-207`); must be a power of 2 (`bwt.c:67-68`). Cost model: `bwt_sa` walks up to `intv-1` LF steps (`bwt.c:86-96`, average ~intv/2 random accesses) per located occurrence. With `-c` raised to tens of thousands, locate cost dominates seeding; `-i 8` quadruples `.sa` (G*2 bytes) but makes each locate ~4x cheaper. `notes:` useful when hammering high-copy seeds.
- **Occ interval**: `OCC_INTERVAL = 128`, compile-time, with a warning not to change it because `bwt_bwt()`/`bwt_occ_intv()` macros hard-code the layout (`bwt.h:36-38,73-75`).
- **Shared memory**: `bwa shm idx` stages the whole index (`bwa_idx2mem`, `bwa.c:362-401`; `bwashm.c:16-89`) so many `bwa mem` processes share one copy; `mem` probes shm first (`fastmap.c:362`).
- **`-6`**: writes `<prefix>.64.*` names (`bwtindex.c:248`); loaders probe `.64.bwt` first (`bwa.c:253`).

## Run-time memory beyond the index

- Batch of reads: `chunk_size(10 Mbp) * n_threads` bases (`fastmap.c:394`) plus per-read `mem_alnreg_v`; per-thread `smem_aux_t` scratch.
- Extension scratch: `ksw_extend2` allocates `qlen*5` profile + `8*(qlen+1)` (`ksw.c:423-424`); `ksw_global2` allocates a `min(qlen, 2w+1) * tlen` backtrack matrix (`ksw.c:548-549`).
- README FAQ 8: memory spikes are usually batches full of centromeric/repeat reads (coordinate-sorted input). A batch of repeat consensus queries with raised `-c` behaves the same way.
{% endraw %}
