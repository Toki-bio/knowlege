# Architecture

Full commit SHA: `3c28777e7e2dcc90f825de1b9f17a89cca7d4452` (tag 2.31, r1302,
2026-05-19, "Release minimap2-2.31 (r1302)"). Repo: https://github.com/lh3/minimap2

## Directory map

Flat C repo; everything is at the top level.

| Path | Role | Key symbols |
|---|---|---|
| `main.c` | CLI: option table, two-pass parsing (preset first), index loop, dispatch to mapping | `long_options` `main.c:24-100`, `main()` `main.c:132-526` |
| `options.c` | Defaults, presets, option validation, `mid_occ` derivation | `mm_idxopt_init` :5, `mm_mapopt_init` :14, `mm_mapopt_update` :69, `mm_set_opt` :91, `mm_check_opt` :202 |
| `minimap.h` | Public API and structs (`mm_idx_t`, `mm_reg1_t`, `mm_mapopt_t`, `MM_F_*` flags) | flags :10-50, structs :77-212 |
| `mmpriv.h` | Private structs (`mm_seed_t`), anchor flag bits, prototypes | `MM_SEED_*` :19-22, `mm_seed_t` :45-51 |
| `sketch.c` | (w,k)-minimizer sketching incl. HPC | `mm_sketch` :77-143 |
| `index.c` | Index build pipeline, bucketed khash, occurrence quantile, dump/load, BED/junction/alt readers | `mm_idx_get` :93, `mm_idx_cal_max_occ` :198, `worker_post` :226, `mm_idx_gen` :389, `mm_idx_dump` :475, `mm_idx_load` :516, `mm_idx_reader_read` :621 |
| `seed.c` | Query-side occurrence filter, seed lookup, high-occurrence seed selection, `rep_len` | `mm_seed_mz_flt` :5, `mm_seed_collect_all` :30, `mm_seed_select` :56, `mm_collect_matches` :98 |
| `map.c` | Per-query pipeline (`mm_map_frag_core`), anchor array construction, multi-threaded file pipeline, output step | `collect_seed_hits` :168, `chain_post` :206, `align_regs` :215, `mm_map_frag_core` :227, `worker_for` :425, `worker_pipeline` :541, `mm_map_file_frag` :663 |
| `lchain.c` | Chaining: O(n*iter) DP chainer and RMQ (minigraph) chainer, backtracking | `comput_sc` :113, `mg_lchain_dp` :148, `mg_lchain_rmq` :250, `mg_chain_backtrack` :27 |
| `hit.c` | Chains -> `mm_reg1_t`, primary/secondary assignment, secondary selection, filters, MAPQ | `mm_gen_regs` :52, `mm_set_parent` :125, `mm_select_sub` :255, `mm_filter_regs` :301, `mm_set_mapq2` :432 |
| `align.c` | Base-level alignment driver around ksw2; seed sanitation; Z-drop re-test; inversion rescue | `mm_align_pair` :336, `mm_align1` :645, `mm_align1_inv` :916, `mm_align_skeleton` :1048 |
| `ksw2.h`, `ksw2_extz2_sse.c`, `ksw2_extd2_sse.c`, `ksw2_exts2_sse.c`, `ksw2_ll_sse.c`, `ksw2_dispatch.c` | SIMD DP: single-affine extension (`extz2`), dual-affine (`extd2`), splice-aware (`exts2`); SSE2/SSE4.1 runtime dispatch | flags `ksw2.h:8-18`, `ksw_apply_zdrop` `ksw2.h:171-182` |
| `esterr.c` | Per-hit divergence estimate (`dv` tag) from minimizer survival | `mm_est_err` :30-64 |
| `format.c` | PAF and SAM writers, cs/ds/MD tags | `write_tags` :397, `mm_write_paf4` :425, `mm_write_sam3` :522 |
| `pe.c` | Paired-end pairing, multi-segment secondary selection | `mm_pair` :81 |
| `jump.c`, `splitidx.c`, `sdust.c`, `bseq.c` | Splice jumps; `--split-prefix` merge; SDUST low-complexity masking (`-T`, effectively unused); FASTA/Q reader with per-batch base budget | `mm_bseq_read3` `bseq.c:80` |
| `kthread.c`, `kalloc.c`, `khash.h`, `ksort.h`, `kvec.h`, `krmq.h` | klib: work-stealing `kt_for`, 3-stage `kt_pipeline`, per-thread arena allocator, hash, radix sort, RMQ tree | `kt_for` `kthread.c:54`, `kt_pipeline` `kthread.c:130` |
| `python/`, `setup.py`, `pyproject.toml` | `mappy` Cython binding | |
| `misc/paftools.js` | PAF utilities (js) | |
| `test/` | Tiny FASTA fixtures (MT-human.fa, MT-orang.fa, x3s-*, inversion toys); no automated test target | |

## Module dependency sketch (top-down)

```
main.c
 |- options.c          (defaults, -x presets, mm_check_opt, mm_mapopt_update -> index.c:mm_idx_cal_max_occ)
 |- index.c            (mm_idx_reader_* -> mm_idx_gen: bseq.c -> sketch.c -> buckets -> worker_post/khash)
 '- map.c  mm_map_file_frag -> kt_pipeline(read | kt_for(worker_for) | write)
       worker_for -> mm_map_frag -> mm_map_frag_core
            |- sketch.c   mm_sketch (query minimizers)        [+ sdust.c if -T]
            |- seed.c     mm_seed_mz_flt, mm_collect_matches -> index.c:mm_idx_get
            |- map.c      collect_seed_hits[_heap] -> sorted anchor array a[]
            |- lchain.c   mg_lchain_dp | mg_lchain_rmq  -> chains u[] + compacted a[]
            |- hit.c      mm_gen_regs -> mm_set_parent -> mm_select_sub      (chain_post)
            |- esterr.c   mm_est_err (dv)
            |- align.c    mm_align_skeleton -> mm_align1 -> mm_align_pair -> ksw2_*   (only with -c/-a)
            |- hit.c      mm_filter_regs, mm_set_parent, mm_select_sub, mm_set_sam_pri, mm_set_mapq2
            '- pe.c       mm_pair (paired-end only)
       write step -> format.c mm_write_paf4 | mm_write_sam3
```

## Entry points

| Command / mode | Entry | Notes |
|---|---|---|
| `minimap2 [-x preset] ref.fa qry.fq` | `main.c:132` -> `mm_idx_reader_read` `index.c:621` -> `mm_mapopt_update` `main.c:465` -> `mm_map_file` `map.c:688` | Outer loop over index parts (`-I`) `main.c:437-506`; query files re-read per part |
| `minimap2 -d ref.mmi ref.fa` | same loop; `mm_idx_dump` via `index.c:631`; no query -> `continue` `main.c:488-491` | Index stores `k,w,b,flag` and overrides later `-k/-w/-H` `index.c:626-627` |
| Paired/multi-segment (`-x sr` with 2 files, `--frag`) | `mm_map_file_frag` `map.c:663` with `MM_F_FRAG_MODE` `main.c:493-500` | Fragments grouped by `mm_qname_same` `map.c:566-571` |
| `--split-prefix` multi-part merge | `mm_split_merge` `map.c:693` after all parts | Hits of each part written to temp files then merged `map.c:476-539` |
| Library | `mm_map` `map.c:392`, `mm_map_frag` `map.c:380`, `mm_idx_str` `index.c:421` | `example.c`, `python/mappy` |

## Design decisions the authors made (and why, if stated)

- **Minimizers, not all k-mers.** Density ~2/(w+1) cuts index size and query
  work; the sketch is strand-symmetric (canonical k-mer by hash, `sketch.c:106-112`)
  and skips palindromic k-mers (`sketch.c:108`). k is capped at 28 (56 bits) and
  w < 256 (`sketch.c:84`).
- **Hash the k-mer before taking the minimum** (`hash64`, `sketch.c:28-38`) so
  poly-A-like low-complexity k-mers are not systematically chosen.
- **Drop frequent minimizers at query time, not index time** (`-f` ->
  `mid_occ`, `options.c:73-79`, `seed.c:106-112`). The index keeps every
  occurrence; filtering is a mapping option so one index serves all presets.
  Rationale in README algorithm overview `README.md:325-328`: frequent seeds are
  "repetitive" and cost time without disambiguating the read's origin.
- **Rescue a few high-occurrence seeds** (`-e occ_dist`, `seed.c:56-96`) so reads
  lying entirely in repeats still get a hit, but cap it hard (`max_max_occ`
  4095, `options.c:37`; 128 per streak, `seed.c:54`).
- **Chaining before alignment; alignment optional.** PAF without `-c` is just
  chains (`map.c:217`). Chaining uses a concave gap cost (linear + 0.5*log2) and
  heuristics `--max-chain-skip`/`--max-chain-iter` to avoid quadratic time
  (`lchain.c:173,182-184`; man page `minimap2.1:275-287`).
- **Primary/secondary by query overlap (`-M`)**, then `-p`/`-N` pruning
  *before* alignment (`map.c:206-213`) so only `1 + N` chains per query are
  ever aligned; this is where repeats are discarded for speed.
- **ksw2 extension with Z-drop**, two-pass (approximate then exact Z-drop,
  `align.c:840-844`) and inversion rescue via `zdrop_inv` (`align.c:864-867`,
  `mm_align1_inv` `align.c:916`).
- **Arena allocator per thread** (`kalloc`, `mm_tbuf_t`) reset when it grows past
  `--cap-kalloc` (`map.c:366-377`) to keep RSS bounded on pathological queries.
