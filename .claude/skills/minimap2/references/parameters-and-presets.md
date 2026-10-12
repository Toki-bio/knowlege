# Parameters and presets

All line numbers at commit `3c28777`. "Default" = `mm_idxopt_init` /
`mm_mapopt_init` (`options.c:5-67`), which is also `-x map-ont` / `lr`
(`options.c:96`).

## Presets (`mm_set_opt`, `options.c:91-193`)

Only the fields a preset changes are listed; everything else stays at default.
`min_dp_max` is `-s`; `mid_occ`/`max_occ` are the `-f` thresholds.

| Preset | k | w | HPC | Scoring a,b / O1,O2 / E1,E2 | zdrop | -n / -m / -s | bw / bw_long | max_gap | -p / -N | mid_occ bounds (-U) | Other | Lines |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| default, `map-ont`, `lr` | 15 | 10 | no | 2,4 / 4,24 / 2,1 | 400,200 | 3 / 40 / 80 | 500 / 20000 | 5000 | 0.8 / 5 | [10, 1e6] | `end_bonus -1`, `min_ksw_len 200`, `occ_dist 500`, `max_max_occ 4095` | `options.c:5-67, 96` |
| `ava-ont` | 15 | 5 | no | default | default | 3 / 100 / 80 | 2000 / 2000 | 5000 | 0.0 / 5 | default | `-P -D --dual=no --no-long-join`, `occ_dist 0` | `97-102` |
| `map-pb`, `map10k` | 19 | 10 | **yes** | default | default | default | default | default | default | default | | `103-104` |
| `ava-pb` | 19 | 5 | yes | default | default | 3 / 100 / 80 | 500 / 500 | 5000 | 0.0 / 5 | default | as ava-ont | `105-110` |
| `lr:hq` | 19 | 19 | no | default | default | default | default | 10000 | default | [50, 500] | | `111-114` |
| `map-hifi`, `map-ccs` | 19 | 19 | no | 1,4 / 6,26 / 2,1 | default | 3 / 40 / **200** | default | 10000 | default | [50, 500] | | `111-118` |
| `lr:hqae` | 25 | 51 | no | default | default | default | default | default | 0.8 / 100 | [50, 500] | `--rmq`, `rmq_inner 5000`, `occ_dist 200`, `chain_gap_scale 5` | `119-126` |
| `map-iclr-prerender` | 15 | 10 | no | 2,6 / 10,50 / 2,1, transition 1 | default | default | default | default | default | default | | `127-130` |
| `map-iclr` | 19 | 10 | no | 2,6 / 10,50 / 2,1, transition 4 | default | default | default | default | default | default | | `131-134` |
| `asm5` | 19 | 19 | no | 1,19 / 39,81 / 3,1 | 200,200 | 3 / 40 / 200 | 1000 / 100000 | 10000 | 0.8 / 50 | [50, 500] | `--rmq` | `135-144` |
| `asm10` | 19 | 19 | no | 1,9 / 16,41 / 2,1 | 200,200 | 3 / 40 / 200 | 1000 / 100000 | 10000 | 0.8 / 50 | [50, 500] | `--rmq` | `135-146` |
| `asm20` | 19 | **10** | no | 1,4 / 6,26 / 2,1 | 200,200 | 3 / 40 / 200 | 1000 / 100000 | 10000 | 0.8 / 50 | [50, 500] | `--rmq` | `135-149` |
| `sr`, `short` | 21 | 11 | no | 2,8 / 12,24 / 2,1 | 100,100 | **2 / 25 / 40** | 100 / 100 | 100 | 0.5 / 20 | `mid_occ = 1000` fixed, `max_occ 5000` | `MM_F_SR\|FRAG_MODE\|NO_PRINT_2ND\|2_IO_THREADS\|HEAP_SORT`, `end_bonus 10`, `max_frag_len 800`, `-K 50M`, FR orientation | `151-168` |
| `splice`, `cdna` | 15 | 5 | no | 1,2 / 2,32 / 1,0 | 200,100 | default | 200000 / 200000 | 2000 (ref 200000) | default | default | splice flags, `noncan 9`, `junc_bonus 9`, `junc_pen 5`, `max_sw_mat 0` | `169-178` |
| `splice:hq` | 15 | 5 | no | 1,4 / 6,24 / 1,0 | 200,100 | default | as splice | as splice | default | default | `noncan 5` | `179-180` |
| `splice:sr` | 15 | 5 | no | 1,4 / 6,24 / 1,0 | 200,100 | 3 / 25 / 40 | as splice | as splice | 0.8 / 10 | default | `min_ksw_len 20`, FR, weak pairing, `-K 100M` | `181-190` |

Notes:
- `-s` default is `min_chain_score * a` = 80 (`options.c:52`), *not* 40. The
  help text prints the live value (`main.c:388`).
- `-x` is applied in the first argv pass (`main.c:150-163`); later flags win.
- `sr` fixes `mid_occ = 1000` (`options.c:166`), so `-f` with a fraction is
  ignored for `sr` (because `mid_occ > 0` short-circuits `mm_mapopt_update`,
  `options.c:73`); an integer `-f N` after `-x sr` still overrides.

## Seed-frequency parameters

| Flag / field | Default | Effect | Code |
|---|---|---|---|
| `-f FLOAT` (<1) -> `mid_occ_frac`, `mid_occ=0` | 2e-4 | `mid_occ` = occurrence at the `(1-f)` quantile of *distinct* minimizers + 1, clamped to `[min_mid_occ, max_mid_occ]` | `main.c:329`, `options.c:73-79`, `index.c:198-220` |
| `-f INT1[,INT2]` (>=1) | | `mid_occ = INT1` directly (no clamp); `max_occ = INT2` enables the second-pass rechain | `main.c:330-331`, `map.c:293-316` |
| `-f 0` | | `mid_occ_frac = 0` -> `cal_max_occ` returns `INT32_MAX` -> clamped to `max_mid_occ` = **1,000,000** | `index.c:204`, `options.c:77-78` |
| `-U lo[,hi]` / `--min-occ-floor lo` | 10, 1000000 | bounds on the derived `mid_occ`; presets `hifi/asm/lr:hq*` use [50,500] | `main.c:322-324, 234`, `options.c:19-20` |
| `-e INT` -> `occ_dist` | 500 | in each run of too-frequent seeds keep `round(run_bp/occ_dist)` lowest-occurrence seeds (<= 128) | `main.c:198`, `seed.c:56-96` |
| `max_max_occ` (no CLI flag) | 4095 | rescued seeds must still have <= this many occurrences | `options.c:37`, `seed.c:89-91` |
| `--q-occ-frac FLOAT` | 0.01 | drop query minimizers that repeat > `mid_occ` times *within the query* and > this fraction; 0 disables | `main.c:254`, `seed.c:5-28` |
| `-T INT` (SDUST) | 0 | mask low-complexity query minimizers; help text hidden, "never used; might be buggy" | `main.c:188, 378`, `map.c:33-57` |
| `-H` | off | HPC minimizers (index-time property) | `main.c:169`, `sketch.c:94-104` |

Diagnostics: `-v3` prints `mid_occ = N` (`options.c:81-82`); `rl:i` tag = query
bp covered by dropped seeds (`seed.c:117-122`, `format.c:447`); `--dbg-seed-occ`
prints each seed's occurrence and filter bit (`main.c:264`, `seed.c:115-116`).

## Chaining parameters

| Flag | Default | Effect | Code |
|---|---|---|---|
| `-g NUM` max_gap | 5000 | max query gap (and ref gap unless `max_gap_ref`/`-F`) between consecutive anchors in a chain; also bounds alignment window extension | `main.c:173`, `map.c:262-271`, `lchain.c:118-120,172`, `align.c:724-728` |
| `-r NUM1[,NUM2]` bw, bw_long | 500, 20000 | max `\|dr - dq\|` between consecutive anchors; `max_drop` in backtracking; alignment band `1.5*bw+1`; `bw_long` for RMQ re-chain and long-join gap filling | `main.c:319-321`, `lchain.c:122,151`, `align.c:660-661,811-813`, `map.c:283-292` |
| `-n INT` min_cnt | 3 | min anchors per chain; also post-alignment filter and split threshold | `main.c:189`, `lchain.c:52,67`, `hit.c:307`, `align.c:864` |
| `-m INT` min_chain_score | 40 | min chaining score; also post-alignment `mlen` filter; `2*m` used by `mm_fix_bad_ends` | `main.c:190`, `lchain.c:36,52,67`, `hit.c:309`, `align.c:675` |
| `--chain-gap-scale` | 0.8 | `chn_pen_gap = scale * 0.01 * k` per bp of diagonal shift | `main.c:247`, `map.c:273` |
| `--chain-skip-scale` | 0.0 | linear penalty per bp of gap length | `main.c:248`, `map.c:274` |
| `--max-chain-skip` | 25 | early exit from predecessor scan | `main.c:221`, `lchain.c:182-184` |
| `--max-chain-iter` | 5000 | max predecessors scanned per anchor | `main.c:222`, `lchain.c:173` |
| `--rmq[=yes\|no]`, `--rmq-inner` | off, 1000 | use `mg_lchain_rmq` (minigraph) chainer; `rmq_size_cap` 100000, rescue size/ratio 1000/0.1 | `main.c:301-303, 259`, `options.c:31-34`, `map.c:275-277` |
| `-M FLOAT` mask_level, `--mask-len`, `--hard-mask-level` | 0.5, inf | secondary if query overlap with a better chain > `M` of the shorter | `main.c:178, 251, 241`, `hit.c:165` |
| `-p FLOAT` pri_ratio | 0.8 | keep secondary if `score >= p * primary` (or within `2k`); `0` disables the filter | `main.c:177`, `hit.c:257,265` |
| `-N INT` best_n | 5 | max secondaries kept per primary (before alignment!) | `main.c:176`, `hit.c:265-267`, `map.c:210` |
| `-P` / `--all-chain` | off | `MM_F_ALL_CHAINS`: no primary/secondary assignment at all; `-p/-N` ignored; incompatible with `--secondary=no` | `main.c:181`, `map.c:208,219`, `options.c:266-270` |
| `-X` | off | `-P -D --dual=no --no-long-join` (all-vs-all) | `main.c:182` |
| `--secondary=yes\|no` | yes | `no` -> `MM_F_NO_PRINT_2ND`: secondaries still computed, just not printed | `main.c:282-283`, `map.c:612-613` |

## Alignment parameters

| Flag | Default | Effect | Code |
|---|---|---|---|
| `-c` / `-a` | off | enable base alignment (`MM_F_CIGAR`); without it no DP, no `-s`/`-z` effects | `main.c:179,183`, `map.c:217` |
| `-A`, `-B`, `-b` (transition), `--score-N` | 2, 4, 0, 1 | match, mismatch, transition mismatch (0 = same as -B), N score; `0 <= score-N < B` | `main.c:191-193,236`, `align.c:11-36`, `options.c:256-260` |
| `-O INT1[,INT2]`, `-E INT1[,INT2]` | 4,24 / 2,1 | dual affine: gap of length g costs `min(O1+g*E1, O2+g*E2)`; need `E1>E2`, `O1+E1<O2+E2`, sum <= 127; `extz2` used when O1==O2 && E1==E2 | `main.c:344-349`, `options.c:246-255`, `align.c:356-359` |
| `-z INT1[,INT2]` zdrop, zdrop_inv | 400, 200 | stop extension when score falls `> zdrop + gap_len*e` below the running max; re-test with `zdrop_inv` for inversions | `main.c:341-343`, `ksw2.h:171-182`, `align.c:843-868` |
| `--end-bonus INT` | -1 (sr: 10) | bonus for reaching the query end in extension (`mqe + end_bonus > max` -> `reach_end`) | `main.c:227`, `ksw2_extz2_sse.c:304-306` |
| `-s INT` min_dp_max | 80 (sr 40, asm/hifi 200) | drop aligned hits whose peak segment score `dp_max` is below | `main.c:194`, `options.c:52`, `hit.c:310` |
| `--min-dp-len` min_ksw_len | 200 | gap filling DP only when both gaps since last DP >= this | `main.c:223`, `align.c:810` |
| `--max-clip-ratio` | 1.0 | drop if both ends clipped above this fraction (1.0 = never) | `main.c:233`, `hit.c:311` |
| `--cap-sw-mem` max_sw_mat | 1e8 cells | skip (mark zdropped) any single DP larger than this | `main.c:242`, `align.c:349-351` |
| `--no-end-flt` | off | skip `mm_fix_bad_ends` | `main.c:240`, `align.c:671` |
| `rank_min_len`, `rank_frac` (no CLI) | 500, 0.9 | divergence-aware `dp_max` rescaling for queries >= 500 bp | `options.c:60-61`, `align.c:1114-1117` |

## Index and resource parameters

| Flag | Default | Effect | Code |
|---|---|---|---|
| `-k`, `-w` | 15, 10 | `k <= 28`, `w < 256` asserted | `main.c:167-168`, `sketch.c:84` |
| `-I NUM` batch_size | 8G | reference bases per index part; multi-part = query files re-read per part, no SAM `@SQ` unless `--split-prefix` | `main.c:196`, `options.c:11`, `index.c:321`, `main.c:445-455` |
| `--bucket-bits` | 14 | `1<<b` buckets, `b <= 2k` | `main.c:216`, `index.c:53` |
| `-d FILE` | | dump index | `main.c:170`, `index.c:631` |
| `--idx-no-seq` | | do not store reference bases (set automatically without `-d`/`-c`/`-a`) | `main.c:229, 352-353` |
| `-t INT` | 3 | mapping threads (`kt_for`) and index post-processing threads | `main.c:138,171` |
| `-K NUM` mini_batch_size | 500M (sr 50M) | query bases per mini-batch; all hits of a batch stay in memory until written | `main.c:197`, `options.c:56`, `map.c:551-552` |
| `-2` | off | 3 pipeline threads instead of 2 (overlap read/write with mapping) | `main.c:201`, `map.c:677` |
| `--cap-kalloc` | 500M | reset a thread's arena when its capacity exceeds this | `main.c:253`, `map.c:371-376` |
| `--max-qlen` | 0 | skip queries longer than this | `main.c:243`, `map.c:244` |
| `--split-prefix` | | map against each index part, merge hits afterwards | `main.c:239`, `map.c:693-736` |

## Memory formula (derived from the code, not measured)

Let `L` = total reference bases, `d` = minimizer density ~ `2/(w+1)` (`w=1`:
1.0, i.e. every non-palindromic k-mer), `N = d*L` minimizers, `D` = distinct
minimizers, `N1` = singletons.

- Build phase: `16 B * N` for bucket arrays (`index.c:306-313`; `kv_push`
  doubling can transiently exceed this), plus `0.5 B * L` for `S` when
  sequence is stored (`index.c:356-361`).
- Final index: `8 B * (N - N1)` positions (`index.c:248`) + khash with
  `n_buckets = next_pow2(D / 0.77)` slots of 16 B (`index.c:246-247`,
  `khash.h:180`) + `0.5 B * L` sequence + per-bucket overhead (`1<<14`
  buckets).
- Peak during `worker_post` is between the two since buckets are converted one
  at a time (`index.c:226-278`).
- Per mapping thread: `36 B * n_anchors` (16 B anchors + 20 B DP arrays,
  `lchain.c:163-166`) plus `n_chains * sizeof(mm_reg1_t)` plus ksw2 buffers
  (`O(band * qlen)`), plus all `mm_reg1_t`/CIGARs of the current `-K` batch.

Worked estimates for a 3 Gbp genome (see `repeat-search-tuning.md` for the
k=11, w=1 case): default k15 w10 -> N ~ 5.5e8 -> ~8.7 GB build, ~4.4 GB
positions + hash (D ~ 4e8 distinct, mostly singletons stored inline so the
`p[]` array is smaller) + 1.5 GB sequence if `-c/-a`. k11 w1 -> N ~ 3e9 ->
~48 GB build, ~24 GB positions (almost nothing is a singleton at k=11),
hash negligible (D <= 2.1e6).
