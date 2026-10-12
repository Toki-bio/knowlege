# Nucleotide lookup tables and seeding

Files: `src/algo/blast/core/blast_nalookup.c`, `include/algo/blast/core/blast_nalookup.h`,
`lookup_wrap.c`, `blast_nascan.c`, `na_ungapped.c`, `blast_options.c`.

## 1. Which table is built: `BlastChooseNaLookupTable` (`blast_nalookup.c:46-196`)

Inputs: `word_size` (W, `-word_size`), `approx_table_entries` (number of
unmasked query positions over both strands, from `EstimateNumTableEntries`,
`lookup_wrap.c:127-128`), `max_q_off`. Output: table type and `lut_width`
(`lut_word_length`, the k-mer actually hashed, <= W).

| W | approx_table_entries | type / lut_width | lines |
|---|---|---|---|
| any with template | - | MB, width = W (discontiguous forces MB) | :65-68 |
| 4,5,6 | any | Small, width W | :76-80 |
| 7 | <250 / else | Small 6 / Small 7 | :82-87 |
| 8 | <8500 / else | Small 7 / Small 8 | :89-94 |
| 9 | <1250 / <21000 / else | Small 7 / Small 8 / MB 9 | :96-107 |
| 10 | <1250 / <8500 / <18000 / else | Small 7 / Small 8 / MB 9 / MB 10 | :109-123 |
| 11 | <12000 / <180000 / else | Small 8 / MB 10 / MB 11 | :125-135 |
| 12 | <8500 / <18000 / <60000 / <900000 / else | Small 8 / MB 9 / MB 10 / MB 11 / MB 12 | :137-153 |
| >=13 (e.g. 28) | <8500 / <300000 / else | Small 8 / MB 11 / MB 12 | :155-164 |

Small -> standard `BlastNaLookupTable` if entries >= 32767 or max offset >= 32768 (:171-176),
or if small-table packing overflows (`lookup_wrap.c:140-149`, `blast_nalookup.c:226-230`).
Note: a 300 bp query (~600 entries) with `-task blastn` (W=11) gets a **Small table of 8-mers**;
with megablast (W=28) also a Small table of 8-mers. Queries are concatenated per
batch (`GetQueryBatchSize`, blastn 100 kb, megablast 5 Mb;
`blast_input_aux.cpp:104-113`), so many queries -> larger entries -> MB table.

Scan step (stride over the subject): `scan_step = W - lut_word_length + 1`
(small `:394`, standard `:566`, MB contiguous `:1346-1348`); discontiguous MB
uses `scan_step = 1` (:1339-1343). Any exact W-mer contains a lut-word at one
of `scan_step` consecutive offsets, so striding is lossless for exact W-mers.
Table sizes: `backbone_size = 4^lut_width` (small `:392`, MB `hashsize` `:1270`,
`hashtable` is `Int4[hashsize]` -> 64 MB for width 12, 16 KB for width 7).

## 2. Building and masking the table

- Words are added only from `lookup_segments` = complement of query masks
  (`BLAST_MainSetUp`, `blast_setup.c:634-637`). Masks = dust (default on, level
  20 window 64 linker 1, `blast_options.c:46-48`; `-dust`), `-lcase_masking`
  (FASTA lowercase, `blast_fasta_input.cpp:397-412`; off by default
  `cmdline_flags.cpp:116`), `-filtering_db`, `-window_masker_taxid/-window_masker_db`
  (`blast_args.cpp:348-367`, invoked in `blast_objmgr_tools.cpp:191-226`).
  `blastn-short` clears filtering (`blast_options_handle.cpp:360`).
- `-soft_masking` (`mask_at_hash`, default **true** for nucleotide,
  `cmdline_flags.cpp:168`): masked bases are kept for extension; only words
  are dropped. With `mask_at_hash` and W > lut_width the table records
  `masked_locations` and every seed is re-checked by `s_IsSeedMasked`
  (`na_ungapped.c:459-491`) inside `s_TypeOfWord` (:508-629) so that the full
  W-mer does not touch a masked word. If `-soft_masking false`, the query
  residues are replaced by N (`BlastSetUp_MaskQuery`, `blast_setup.c:621-623`).
- Words containing ambiguity codes (`val & BLAST2NA_MASK`) are never indexed
  (`blast_nalookup.c:1046-1050`, :757-762).
- Unmasked regions shorter than W add nothing (:1036-1037).
- **Word-frequency / repeat filter:** `LookupTableOptions::db_filter` +
  `max_db_word_count` (`blast_options.h:217-223`). When set, `s_FillPV` marks
  query words, `s_ScanSubjectForWordCounts` scans the *whole database* once
  counting 16-mers with 4-bit saturating counters (`:1138-1240`), and
  `s_FillContigMBTable` skips words whose count >= `max_db_word_count`
  (:1053-1064). Only reachable through `-limit_lookup` of the mapping
  (magicblast) arguments (`blast_args.cpp:1495-1502,1575-1576`), default
  count 30 (`MAX_DB_WORD_COUNT_MAPPER`, `blast_options.h:174`), and
  validation rejects it for word size < 16 (`blast_options.c:1424-1429`).
  **blastn/megablast/dc-megablast have no such cutoff**: a high-copy
  repeat word produces a hit at every genomic copy. Poly-A/poly-T word
  removal (`s_RemovePolyAWords` :919-962) is also mapping-only (:1118-1120).
- Lookup table options validation (`LookupTableOptionsValidate`
  `blast_options.c:1284-1432`): nucleotide word size >= 4 (:1318-1324), <= 100
  (`DBSEQ_CHUNK_OVERLAP`, :1325-1333).

## 3. Discontiguous megablast templates

- Allowed: word size 11 or 12; template length 16, 18 or 21; type 0 coding,
  1 optimal, 2 both (`s_DiscWordOptionsValidate` `blast_options.c:1248-1282`;
  CLI constraints `blast_args.cpp:709-736`, `-template_type coding|optimal|coding_and_optimal`,
  `-template_length 16|18|21`).
- Template -> enum: `s_GetDiscTemplateType` (`blast_nalookup.c:606-651`);
  patterns drawn in `blast_nalookup.h:191-215`, e.g. 11-of-18 coding
  `10,110,110,010,110,110,1`, 11-of-18 optimal `111,010,010,110,010,111`.
- `s_FillDiscMBTable` (:659-843): `ComputeDiscontiguousIndex` (:783) masks the
  template's 1-positions out of an 18/21/16-base window into a 22- or 24-bit
  key; `coding_and_optimal` fills `hashtable2` with the next enum (coding+1 =
  optimal) (:690-696). Scanners: `s_MB_DiscWordScanSubject_11_18_1`,
  `_11_21_1`, generic `_1`, two-template `_TwoTemplates_1`
  (`blast_nascan.c:2619-2628`), stride 1.
- For discontiguous tables `word_length == lut_word_length == template_length`
  in the word finder (`na_ungapped.c:1633-1636`) and extension is
  `s_BlastNaExtendDirect` (:1807-1808): no exact-extension step; the hit goes
  straight to the two-hit/ungapped stage. dc-megablast default window 40.
- Database index cannot be combined with templates (`setup_factory.cpp:444-447`).

## 4. Scanning the subject (`blast_nascan.c`)

`BlastChooseNucleotideScanSubject` (:2982-2992) picks a specialised scanner by
table type, width and stride: small table (:1298-1375, e.g. `_8_4` for width 8
stride 4 = `-task blastn` on a short query, `_7_1` for `blastn-short`), MB
(:2602-2691, `_11_1Mod4`, `_12`/`16` generic). All read the 2-bit packed
subject; `s_MBScanSubject_Any` (:1482-1626) shows the generic pattern
(shift/mask a 32-bit window, `PV_TEST`, walk the `next_pos` chain). Hits are
written to `offset_pairs` (max `OFFSET_ARRAY_SIZE` 4096 + longest chain,
`lookup_wrap.c:232-269`).

## 5. Exact extension to the full word (`na_ungapped.c`)

`BlastChooseNaExtend` (:1799-1845): if `lut_word_length == word_length` or
discontiguous -> `s_BlastNaExtendDirect` (:972-1037); else aligned variants or
`s_BlastNaExtend` (:1060-1202) / `s_BlastSmallNaExtend` (:1494). `s_BlastNaExtend`
extends the lut hit left then right base-by-base and drops it unless
`ext_left + ext_right >= word_length - lut_word_length` (:1109-1146): the
**exact** W-mer requirement. Masked-seed and double-word checks happen in
`s_TypeOfWord` (:508-629).

## 6. Word finder loop

`BlastNaWordFinder` (`na_ungapped.c:1605-1703`): for each scan range
(`s_DetermineScanningOffsets`, honouring subject soft-mask ranges
`seq_ranges`), `scansub(...)` then `extend(...)`; afterwards
`Blast_InitHitListSortByScore`. With database soft masking
(`-db_soft_mask`) the generic scanner/extender are used (:1659-1676) so hits
never start in masked subject ranges.
`MB_IndexedWordFinder` (:1705-1797) replaces scanning with index lookups for
indexed OIDs and falls back to `BlastNaWordFinder` otherwise (:1737-1741).
