# Pitfalls and gotchas

- **`-g -3` is probably not what the manual says.** In FASTA36 `-g` takes a
  numeric gap-extension penalty; `-g -3` sets it to -3 and `-3` is never seen
  as the "forward strand only" switch. SINEderella `MANUAL.md:1218-1230`
  documents `-g` as "Global/gapped alignment" and `-3` as "forward + reverse
  complement". Both strands are searched regardless (that is the default), so
  results are unaffected, but anyone re-implementing "the same parameters"
  elsewhere should test which reading holds (`sear:396`). Unverified: ssearch36
  is not installed in this container.
- **`-E 2` has no `-Z`.** E-values are computed against each part file
  (~1/20 of the windows), so the effective threshold tightens with genome size.
  SINEderella's flankscan uses `-Z 20000` for exactly this reason
  (`flankscan/fs3_partners.sh:85-91`); `sear` never did.
- **Length filter column.** The ssearch36 branch keeps a window when the
  **longest merged hit length** (`-o max` of col 4) >= 0.8 x query
  (`sear:405-407`); the minimap2 branch uses the merged span `$3-$2`
  (`sear:337`). Commit `85421e7` ("filter on merged span ... was checking wrong
  column") suggests the intent is the span; check `plus_to_add_*` on a real
  run before quoting either.
- **Output prefix from the path.** `taxname=$(echo "$WORKING_BANK" | cut -c1-3)`
  (`sear:251`): with an absolute bank path in minimap2 mode outputs are named
  `/ho-<q>.bed`. The ssearch36 branch overwrites it from the part file name
  (`sear:387`). Fixed in SINEderella's copy with `basename` (`sear:254`,
  `docs/AUDIT_2026-10-05.md` D5).
- **Cascade text vs code.** Usage and log say `asm5 -> asm10 -> asm20`
  (`sear:36,275,291`); the loop runs `asm10 asm20` (`sear:293`).
- **`-m` means the opposite in `sear2k` and `sear`.** `sear2k -m` skips minus
  rounds (`sear2k:54-56`); `sear -m` enables them (`sear:81-83`). `sear2k`
  runs depletion rounds by default.
- **Multi-record query files are mis-sized.** `len` sums all records
  (`sear:245`); a two-consensus file gets `flen` from the sum and no hit will
  pass. Use one record per file or `sear_multi`.
- **Header characters.** `_` and tabs in genome headers force a cleaned copy;
  `:`, `(`, `)` in headers are rewritten in minus banks (`sear:444`). ssearch36
  truncates identifiers (SINEderella `tools/length_variants.py:374` says 60
  characters) and cannot open paths longer than ~120 characters
  (`docs/SATELLITES.md:327-332`); `sear` runs in cwd with short names so it is
  not hit, but a replacement tool launched from deep paths would be.
- **Split reuse is directory-global.** `compgen -G "*.2k.part_*.bnk"`
  (`sear:369`) picks up parts of *any* bank in cwd; two genomes in one
  directory will mix.
- **`-s/--stop`** counts filtered hits across parts in file order and stops
  scanning (`sear:419-424`): the result is then not the whole genome. Not used
  by SINEderella.
- **MAFFT `--threadit $T`** makes `-b` alignments non-reproducible
  (SINEderella `PUBLISH_WORKFLOW.md:123`).
- **Lowercase / soft-masked genomes do nothing.** ssearch36 folds DNA
  libraries to uppercase; `-S` is a no-op for DNA in the shipped build
  (`docs/MASKING.md:36-42`). Mask with N (`tools/mask_genome.py`) if needed.

## Defaults that surprise

| Parameter | Default | Why it matters |
|---|---|---|
| Slf | 0.8 (`sear`), 0.9 (`sear2k`) | 97-99 % of RepeatMasker Alu copies missed on chr21 are fragments < 240 bp (`docs/RELEASE_CHECK_2026-10-05.md:30`); lower it (0.2 was used for satellite monomers, `docs/SATELLITES.md:163`) to see fragments |
| Homology | 65 % | identity of the best merged hit, not of the whole locus |
| Window / step | 2500 / 2000 | copies > 500 bp can be split across windows and are rescued only by the second merge |
| Parts | 20 | fixed; not tied to thread count |
| `-E` | 2 | per part, no `-Z` |
| Best flank | 50/50 | `--slop` changes both main and best flanks |
| `-b` N | 100 (ssearch36), 10 / 5 (minimap2, < 5 kb / >= 5 kb) | `sear:52,278-285` |

## Known bugs / open upstream issues (as of pinned commit)

- taxname path bug and cascade text (above); both fixed only in the
  SINEderella copy.
- No tests, no issue tracker entries visible in the shallow clone.
