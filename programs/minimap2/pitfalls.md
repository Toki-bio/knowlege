---
title: "minimap2: Pitfalls and gotchas"
type: reference
generated: true
---

> Part of the [minimap2](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/minimap2/references/pitfalls.md)

{% raw %}
# Pitfalls and gotchas

All line numbers at commit `3c28777`.

- **`-f 0` is not "no filter".** `-f FLOAT < 1` sets `mid_occ_frac` and zeroes
  `mid_occ` (`main.c:329`); `mm_idx_cal_max_occ` returns `INT32_MAX` for
  `f <= 0` (`index.c:204`) and `mm_mapopt_update` clamps to `max_mid_occ` =
  1,000,000 (`options.c:77-78`). Use `-f INT` (>= 1) for an exact cutoff
  (`main.c:330`) or raise `-U lo,hi` (`main.c:322-324`).
- **`-s` default is 80, not 40.** `min_dp_max = min_chain_score * a`
  (`options.c:52`); `sr` sets 40, `asm*`/`map-hifi` set 200. The help prints
  the live value (`main.c:388`).
- **The one-seed rescue rounds.** `max_high_occ = (int)((pe - ps) / dist +
  .499)` (`seed.c:73`): a 300 bp query of all-repetitive seeds rescues exactly
  1 seed at `-e 500`; a 240 bp query rescues 0. Rescued seeds must have
  `<= max_max_occ` 4095 hits (`seed.c:89-91`), a constant with no CLI flag
  (`options.c:37`). Rescue is skipped entirely when `occ_dist == 0` (ava
  presets, `options.c:102,110`) (`seed.c:106-111`).
- **`-N` is applied before alignment** (`chain_post`, `map.c:206-213`) and
  again after (`map.c:219-223`); `-p` and `-N` do nothing with `-P`
  (`map.c:208`) and `-p 0` disables `mm_select_sub` completely (`hit.c:257`),
  not just the ratio test. `-P` + `--secondary=no` is an error
  (`options.c:266-270`). `-N 0` is rewritten to the previous `-N` plus
  `--secondary=no` (`main.c:356-359`).
- **`-x` must come first in spirit**: presets are applied in a separate first
  pass (`main.c:150-163`) so any explicit flag wins regardless of order; but
  `-x sr` fixes `mid_occ = 1000` (`options.c:166`), which makes a later `-f
  0.0002` a no-op because `mm_mapopt_update` only recomputes when
  `mid_occ <= 0` (`options.c:73`).
- **A `.mmi` overrides `-k/-w/-H`** with only a warning (`index.c:626-627`);
  README `README.md:126-131`.
- **No sequence in the index without `-d`, `-c` or `-a`** (`main.c:352-353`),
  so a PAF-only run cannot be upgraded to CIGAR by re-using an index built
  that way (`main.c:439-444` errors "doesn't contain sequences").
- **`-m` is also a post-alignment filter on `mlen`** (exact matched bases from
  CIGAR, `hit.c:309`) and `-n` on anchor count (`hit.c:307`); `2 * -m` is the
  target of end trimming (`align.c:675`). Lowering `-m` therefore changes
  alignment too.
- **Gap filling between seeds is skipped for gaps < `min_ksw_len` 200**
  (`align.c:810`) unless it is the last seed: for a 300 bp query there is
  effectively one global DP from the first to the last anchor plus two
  extensions, so interior seed quality matters little but the first/last
  anchor positions matter a lot (`mm_fix_bad_ends`, `align.c:527-561`).
- **Scoring matrices are int8 and gap costs are checked:** `(O1+E1)+(O2+E2) <=
  127`, `E1 > E2`, `O1+E1 < O2+E2` when dual (`options.c:246-255`),
  `0 <= --score-N < -B` (`options.c:256-260`), `zdrop >= zdrop_inv`
  (`options.c:261-265`). `asm20` has break-even identity 80 % (`a=1,b=4`,
  `options.c:148`) despite the name suggesting 20 % divergence.
- **Chaining distance parameters are not symmetric in the code.** In
  `comput_sc` the query gap is tested against `max_dist_x` *and* `max_dist_y`
  (`lchain.c:118,120`) while the reference window is bounded via `st` by
  `max_dist_x` (`lchain.c:172`); both default to `-g`, but `-F`/`max_gap_ref`
  only change the reference side (`map.c:266-271`).
- **`mm_update_dp_max` silently rescales scores** for queries >= 500 bp when the
  second-best is within 0.9 of the best (`align.c:1114-1117, 1022-1046`),
  changing which hits survive `-s`; disabled by `-P`, `--sr`, `--split-prefix`.
- **All hits of a `-K` batch live in memory until the batch is written**
  (`map.c:578-641`): with many hits per query (repeats, `-P`) lower `-K`
  (default 500 Mbp of query, `options.c:56`).
- **`--print-seeds` / `--print-aln-seq` force `-t 1`** (`main.c:220,224`).
- **Reverse-strand anchor coordinates are flipped** onto the query's reverse
  complement (`map.c:190`) and the query strand mode `--qstrand` flips the
  reference instead (`map.c:191-194`); `--qstrand` is rejected with `-a`,
  `-H`, `--frag`, `--splice` (`options.c:271-275`).
- **Symmetric (palindromic) k-mers are never minimizers** (`sketch.c:108`);
  odd `k` avoids them entirely.
- **`-T` (SDUST) is undocumented and flagged "never used; might be buggy"**
  (`main.c:378`).
- **Index file is not portable** across endianness/architectures: raw structs
  (`index.c:475-514`).
- **MAPQ is depressed by repetitive seeds** through `uniq_ratio = sum_sc /
  (sum_sc + rep_len)` (`hit.c:445, 452`), so a correct unique hit of a query
  that *also* contains a repeat gets a low MAPQ.

## Defaults that surprise

| Parameter | Default | Why it matters |
|---|---|---|
| `-s` (`min_dp_max`) | 80 (= `-m 40` x `-A 2`), `options.c:52` | docs/folklore say 40; 300 bp alignments below ~72 % identity are dropped |
| `-f` -> `mid_occ` clamp | `[10, 1e6]`, `options.c:19-20` | `-f 0` = 1e6, not infinite |
| `occ_dist` `-e` | 500, `options.c:38` | one rescued repeat seed per 500 bp of query |
| `max_max_occ` | 4095, `options.c:37` | rescued seeds must be rarer than this; no flag |
| `-I` | 8 Gbp, `options.c:11` | README/FAQ still say 4G (`README.md:319`, `FAQ.md:16`) |
| `-K` | 500 Mbp, `options.c:56` | README says 200M (`README.md:322`) |
| `-g` | 5000, `options.c:27` | man page says 10k (`minimap2.1:170-173`) |
| `-r` | 500,20000, `options.c:26` | alignment band is `1.5 * bw + 1` (`align.c:660`), i.e. 751 |
| `--chain-gap-scale` | 0.8, `options.c:35` | man page says 1.0 (`minimap2.1:289-290`) |
| `end_bonus` | -1, `options.c:51` | extension never favours reaching the query end unless `--end-bonus > 0` |
| `-t` | 3, `main.c:138` | |
| `min_ksw_len` | 200, `options.c:53` | interior gap filling only across >= 200 bp |
| `rank_min_len`/`rank_frac` | 500 / 0.9, `options.c:60-61` | hidden `dp_max` rescaling for queries >= 500 bp |

## Known bugs / open upstream issues (as of pinned commit)

- `NEWS.md` for 2.31 (lines 1-30): fixed long-standing bugs in supplementary/
  secondary flagging and an out-of-bounds access in inversion Smith-Waterman;
  output may differ from 2.30 in rare cases.
- `mm_idx_reader_eof` comment: "in extremely rare cases, mm_bseq_eof() might
  not work" (`index.c:637`).
- `mm_align1` TODO: Z-drop sometimes triggers because the next seed placement
  is wrong (`align.c:848`).
- `mg_lchain_dp` TODO: not verified for more than 2^31 anchors (`lchain.c:150`).
- `-T` SDUST path untested (`main.c:378`).
{% endraw %}
