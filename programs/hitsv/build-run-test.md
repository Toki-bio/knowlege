---
title: "hitsv: Build, run, test"
type: reference
generated: true
---

> Part of the [hitsv](../) digest · [source](https://github.com/Toki-bio/knowlege/blob/main/.claude/skills/hitsv/references/build-run-test.md)

{% raw %}
# Build, run, test

Verified on Linux x86_64 at commit `93880367` (2026-10-16).

## Build

```sh
git clone --depth 1 https://github.com/hitbc/HitSV
cd HitSV/Release
make clean          # required: 70 stale .o/.d files are tracked in git under Release/src/
make all -j 8       # ~1-2 min; g++ -O3, links -lpthread -lz -static (Release/makefile:68)
./HitSV             # prints the command list
```

- Only zlib and pthreads are needed; htslib is vendored and compiled in (`src/htslib/`).
- Static linking prints a harmless `getaddrinfo` warning from `knetfile.c`.
- The self-built binary differs bytewise from the shipped `HitSV_2.0.4-release`, which is expected
  (different toolchain). Version strings say 2.0.0 (`src/clib/desc.h:2`), even at v2.0.4.
- No test suite. `src/function_test.cpp` is developer scratch code, not a runner.

## Minimal runs that exercise the main paths

```sh
D=../demo
# LRS (HiFi) on the 1 Mb demo slice: ~0.1 s
./HitSV call -S 0 -E 0 -s 0 -F 1000000 -p HIFI \
  -r $D/hs37d5_1_0_1000000.fa -l $D/LRS_HG002_1_869000_870000_10X_demo.bam -o lrs.vcf 2> call.log

# SRS prerequisites (the demo SRS BAM and its JSON are NOT in the repo, only the .bai)
./HitSV srs_fa_stat ref.fa > ref.stat.txt            # ~1.1 s per Mb single-threaded on the demo slice
./HitSV ngs_bam_stat ref.fa srs.bam > srs.stat.json   # NOT "srs_bam_stat"
./HitSV srs_trans_reads ref.fa srs.bam TL.bam && samtools sort -o TL.sort.bam TL.bam && samtools index TL.sort.bam
./HitSV call -n srs.bam -I ref.stat.txt -T srs.stat.json -L TL.sort.bam -P asm10 -r ref.fa -o srs.vcf

# Hybrid
./HitSV call -l lrs.bam -n srs.bam -I ref.stat.txt -T srs.stat.json -L TL.sort.bam -p HIFI -P asm10 -r ref.fa -o hyb.vcf
```

Expected demo LRS output (tested): 16 records in `1:868891-872162`, including an 839 bp DEL in a
VNTR at 869406 and a 39 bp DEL at 871644, the rest SNP/INDEL small variants (`-v 1` default).
All PASS, QUAL 253/353.

Size the `srs_fa_stat` cost before running it on a whole genome. The demo slice suggests roughly an
hour per Gb on one core, but this was measured on a single 1 Mb human slice and is not a benchmark.
There are no threads, so split the FASTA by record if you need it faster. Keep the record order,
because the file stores chromosome indices, not names.

## Parallel runs

- The built-in `HitSV_multy_process.sh` works for `LRS` mode. In `SRS`/`Hybrid` modes, fix line 145
  (`srs_bam_stat` → `ngs_bam_stat`) first.
- Its `#SBATCH` header hard-codes the author's log path and partition (`HitSV_multy_process.sh:2-5`).
  Replace these for Monsoon or KIT.
- The script iterates `seq 0 MAX_CHR` (default 2000) and merges with `grep -v "#"` per part
  (`:192-198`). Pass `-c <n_contigs-1>` for scaffold-level assemblies.
- An alternative is your own SLURM array with `-S i -E i -H`, then `cat` the parts.

## Debugging tips

- Logs go to stderr. The README examples send them to `/dev/null`, which hides useful messages such
  as "Assembly will be skipped because over depth in this region" (`SV_core_SRS_ins_del.cpp:907-914`).
  Keep stderr while testing.
- Use `-S/-E/-s/-F` to reproduce a single locus. Indices are 0-based BAM header indices.
- `-b region.bed -B PURE_STR` dumps the local contigs for a region (LRS/ASM only). This is the
  quickest way to see what the assembler built at a SINE or satellite locus.
- First breakpoints: `SV_core.hpp:2355` (block pipeline), `SV_core_LRS.cpp:2042` (clip signal),
  `SV_core_SRS_ins_del.cpp:2921` (ksw2 realignment setup).
{% endraw %}
