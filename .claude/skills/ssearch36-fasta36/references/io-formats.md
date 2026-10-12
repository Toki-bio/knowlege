# Input and output formats

## Inputs accepted

| Format | Reader | Notes / quirks |
|---|---|---|
| FASTA (query and library) | `agetlib()` `src/nmgetlib.c:574-705`; mmap variant `src/mmgetaa.c` | lines starting `#` ignored; `@` = stdin query, `file:50-150` subranges (`fasta36.1:432-447`); `>`-less first line fails |
| `@file.list` / FASTLIBS | `src/lib_sel.c` | indirect lists of library files; `-l` |
| BLAST2 formatted DB, GCG, EMBL, PIR | `src/ncbl2_mlib.c`, other `*getlib` in `nmgetlib.c` | library type digit after file name or in FASTLIBS |
| PSSM | `-P file` (`initfa.c:1238-1246`) | ssearch/gg/glsearch only |
| lower-case masking | `-S` (`initfa.c:1288-1291`) | lower case = X in scan only |

Query type is guessed from composition unless `-n`/`-p` (`fasta36.1:106-110`);
stdin queries need `-n`.

## Output format selection (`-m`, parsed in `src/doinit.c:875-955`)

| `-m` | markx bits (`src/defs.h:145-163`) | What you get | Alignments printed? |
|---|---|---|---|
| 0,1,2,3 | `MX_ATYPE` variants | classic hit list + alignments | yes (`-d` limits) |
| 4,5 | `MX_AMAP` | alignment "map" | yes |
| 6 | `MX_HTML` | HTML | yes |
| 8 | `MX_M8OUT`; sets `ashow=0`, `std_output=0` (`doinit.c:923-926`) | BLAST tabular, one row per HSP | no |
| 8C | + `MX_M8COMMENT` (`927`) | tabular with `# Fields:` header and `# N hits found` (`mshowbest.c:250-298`) | no |
| 8Cc / 8Cd / 8CC / 8CD / 8CB | `show_code` ALIGN / ALIGN+EXT / CIGAR / CIGAR+EXT / BTOP (`928-932`) | extra `aln_code`/`BTOP` column | no |
| 8CB**l** / **L** / **s** / **r** / **f** / **Z** | `MX_M8_BTAB_LEN`, `+DOMINFO`, `_SIM`, `_RAW`, `_ADD_LFILE`, `_DB_SIZE` (`934-955`) | add q/s lengths, similarity, raw score, library file name, db size | no |
| 9 / 9c / 9C / 9d / 9D / 9B / 9i / 9I | `MX_M9SUMM` + `show_code` (`910-918`) | normal output; hit list gains coordinate columns (+ encoded alignment) | yes (`-d 0` to suppress) |
| 10 | `MX_M10FORM` | `>>name` blocks with `; key: value` lines per HSP | yes (as `;`-tagged) |
| 11 | `MX_M11OUT` | lav (lalign36) | |
| B / BB | `MX_MBLAST` / `+MX_MBLAST2` (`878-889`) | BLAST-like alignments / full BLAST mimic | yes |
| `"F8C out.txt"` | separate file (`doinit.c:~760`) | several formats from one run (`fasta36.1:284-287`) | |

Formats combine: `-m 1 -m 9c -m 6` (`fasta36.1:270-271`).

## `-m 8` / `-m 8C` columns (`src/mshowbest.c:580-586,657-676`)

`query id, subject id, % identity, alignment length, mismatches, gap opens,
q. start, q. end, s. start, s. end, evalue, bit score` [+ `raw score` with `r`]
[+ `eval2` when `-z 2x`] [+ `aln_code`/`BTOP` with `c/C/B`] [+ annotations].
With `l`: `query id, query length, subject id, subject length, ...`
(`582-583`). Subject id is the first word of the description plus the `:k`
chunk suffix and `_r` strand suffix (`mshowbest.c:482-496`). Start/end are
whole-sequence coordinates (`d_start1/d_stop1` from `calc_coord`,
`633-636`); for reverse-strand hits the *query* coordinates are reported on the
reverse-complemented query (`frame == 1`, `comp_lib9.c:1133-1141`). One row per
HSP: the `do {...} while (cur_ares_p = cur_ares_p->next)` loop (`504-741`).
Gaps and `%identity` include gap positions by default (`-XG` for ungapped,
`doc/readme.v36:57-62`).

## `-m 9c` columns appended to the hit line (`mshowbest.c:641-650`)

`%_id %_sim sw alen an0 ax0 pn0 px0 an1 ax1 pn1 px1 gapq gapl fs [code]` =
identity, similarity, SW score, aligned length, query align start/stop, query
sequence start/end (full length ends), library align start/stop, library
sequence start/end, gaps in query, gaps in library, frameshifts, then the
encoded alignment (`c`: `=12+3-2=...` style, `C`: CIGAR with `F`/`R` for
frameshifts, `B`: BTOP; `fasta36.1:263-269`, `fasta_guide.tex:957-989`).
Secondary HSPs are separate lines starting `+-` (`mshowbest.c:548`).

## `-m 10` (`src/mshowalign2.c:588-608,761-786`; `src/c_dispn.c:155,176`)

Per hit `>>name description` then per HSP `; sw_s-w opt: N`, `; sw_z-score:`,
`; sw_bits:`, `; sw_expect:`, `; sw_ident:`, `; sw_sim:`, `; sw_overlap:`,
followed by `>name` and `; sq_len`, `; sq_offset`, `; sq_type`,
`; al_start: ...`, `; al_stop: ...`, `; al_display_start:` for query then
library (`c_dispn.c:155,176`), then the aligned strings. Guide: "no longer
reliably tested; `-m 8CBL` is easier to parse" (`fasta_guide.tex:991`).

## Getting all suboptimal alignments with coordinates on both sides

1. Keep `do_rep` on (default) and raise `e_cut_r`: `-E "0.01 1"` makes both
   thresholds 0.01 (`initfa.c:1138`); with the DNA default `e_cut_r` is 0.001
   (`initfa.c:1586`), so secondary copies need E() <= 1e-3.
2. Use `-m 8CBl` (q/s lengths + BTOP) or `-m 9c -d 0` for one line per HSP with
   both coordinate pairs; `-m 8CBr` adds the raw score so you can re-threshold.
3. Use `-b =N` / `-b $` to override the E-based hit count when E-values are
   unreliable (`initfa.c:1058-1080`), `-F` to drop the strongest.
4. Remember the `:k` suffix: HSPs in different chunks of the same chromosome
   appear under different subject ids; merge on the name before the `:` and
   de-duplicate hits whose coordinates fall in the 150-residue overlap
   (`fasta_guide.tex:2126-2136`).
5. `-Xg` keeps `gi|` prefixes; `-C N` sets the name width in alignments
   (`fasta36.1:170-172`). `-w`, `-W` control alignment line width/context.

## Header lines worth parsing

Statistics summary: `Expectation_n fit: rho(ln(x))= ... mu= ... mean_var= ...
Lambda= ...` (`scaleswn.c:531-534`); algorithm line `Smith-Waterman (SSE2,
Michael Farrar 2006) (36.3.8j1 ...)` and `+5/-4 matrix (5:-4), open/ext: -12/-4`
(`dropgsw2.c:662-691`); for fasta36 `ktup: 6, E-join: 0.25 (frac), E-opt: 0.05
(frac), width: 16` (`dropnfa.c:454-480`).
