---
name: hpsmc
description: What hPSMC (pseudo-diploid PSMC, Cahill et al. 2016) measures, how to run it, how to read its split time, what it was and was not validated on (apes and bears yes, human populations no direct test found), its limits in divergence time and phylogenetic distance, and a checklist for judging a run, including the Sicista betulina vs S. trizona case. Use when planning, running, reviewing or interpreting an hPSMC split-time estimate, or when comparing it with BPP, MSMC2 cross-coalescence or fossil dates.
---

# hPSMC: what it dates and where it breaks

> Type: notes (literature research with evidence labels). Written 2026-10-14. Full report with all links: `references/research-report-2026-10-14.md`.

Labels used: **[direct]** read in full by the researcher; **[summary]** tool-generated summary only, check the PDF before quoting; **[inference]** reasoning, not a published result; **[unverified]** not confirmed.

## What it is

- Cahill, Soares, Green, Shapiro 2016, Phil Trans R Soc B 371(1699):20150138 (not MBE, not Mol Ecol). Code: github.com/jacahill/hPSMC. [direct]
- Idea: haploidize one genome per lineage, merge into a fake diploid, run PSMC. Ne rises sharply where the two lineages can no longer coalesce, which is read as the split. The idea predates it: Li and Durbin 2011 used Yoruba plus non-African male X chromosomes. [direct]
- README pipeline [direct]: `samtools mpileup -q30 -Q60` piped to `pu2fa`; merge with `psmcfa_from_2_fastas.py -b10 -m5`; `psmc -N25 -t15 -r5 -p "4+25*2+4+6"`. The command-line flags for bin size, mutation rate and generation time can silently fall back to defaults.
- The time is not read by eye. Simulate clean splits with ms (`hPSMC_quantify_split_time.py`), keep simulations whose curve crosses the data between 1.5x and 10x the pre-divergence Ne, report the interval. [direct]
- It dates the **end of gene flow**, not the first split. A review calls it a maximum bound on divergence. [direct]

## Validated on

- Simulations (split 100 ka to 5 Ma, pre-divergence Ne 1,000-50,000, 0.1-1,000 migrants per generation at 1 Ma) and great apes plus *Ursus* bears. Details beyond the abstract are [summary].
- **Human populations: no direct test found** (limited search, about 25 calls; not exhaustive) [unverified negative]. Closest: Li and Durbin 2011 X pseudo-diploid (gene exchange until about 20-40 ka rather than a clean split) [direct]; Prüfer 2014 two-chromosome PSMC for modern human vs Neandertal/Denisovan [summary]; MSMC/MSMC2 cross-coalescence shows gradual human separation [summary].
- No benchmark of hPSMC against MSMC2, SMC++, BPP or G-PhoCS was found.

## Limits

- No published numeric ceiling in years, per-bp divergence or heterozygous sites per window [direct, absence in sources read].
- Simulated range 100 ka to 5 Ma [summary]. Applications reach 7-14 Ma (sea turtles) and about 10 Ma (birds-of-paradise) without simulation validation [direct]; treat those as the oldest signal, not shown accuracy [inference].
- **Gene flow:** one or more migrants per generation erases the signal; 0.1 detected less precisely [summary].
- **Ancestral Ne:** above about 10,000 pushes dates older [summary]. ILS and gene flow are not separable by hPSMC alone (reviewer opinion).
- **Bins:** cobraa README: theta above 0.01 makes the 100-bp bin inappropriate; use 5-10 bp at theta 0.05 [direct]. Fraction of heterozygous 100-bp bins is 1-(1-d)^100: 63% at d=1%, 87% at 2% [my calculation]. Use `-b10` for old splits.
- **Scaling:** dates scale as 1/mu and linearly with generation time.
- **Mapping bias:** PSMC effects beyond about 0.14% reference difference (preprint) [direct]; no source gives a "1-5% tolerance" [unverified]; no hPSMC-specific mapping-bias test.
- **Alternatives and their stated ranges:** MSMC2/MSMC-IM about 10 ka to 3 Myr [direct]; cobraa about 10 ka to 5 Ma [direct]; diCal2 validated on 10-70 ka only.

## Checklist before trusting a split time

1. Same documented mu and generation time as any method it is compared with; report a sensitivity grid, not one placeholder-based date.
2. `-b10` used; test other `-p` patterns for false peaks.
3. ms simulation bracketing at your own Ne, mu, generation time, with and without low migration (0.1 per generation). Report the 1.5x-10x interval.
4. Autosomes only. Remove X (males are hemizygous; X/Y heterozygote calls are errors) and short contigs (PSMC tolerated about 100 kb scaffolds, lost resolution at 10 kb) [direct, one preprint].
5. Mask ROH and repeats as sensitivity runs; check ONT or WGA artefact sources of false differences [inference, no quantification published].
6. State the mapping reference; compare mapping to each species' assembly if both exist.
7. Expect hPSMC younger than a no-gene-flow multispecies-coalescent date when there was early gene flow; that is consistent with gene flow but does not prove it.

## Case: Sicista betulina vs S. trizona (state 2026-10-13)

- Their numbers: hPSMC 0.7-0.85 Myr, BPP 1.19 Myr, with placeholder mu 5.4e-9 per generation and 1 yr per generation (C57BL/6 lab value, Uchimura 2015) [direct].
- "Published 2.6 Mya" is a fixed BEAST calibration prior (mean 2.6, SD 0.1) on the subtilis-group vs betulina-group node, Deng et al. 2025 [direct] and Zhu et al. 2025 [summary]. Not a betulina-vs-trizona estimate, no CI, the two papers cite different provenance (Kowalski 2001; Lebedev 2019), neither opened. It cannot judge either method.
- If the real generation time is 1.5-2 yr, 0.7-0.85 Myr becomes about 1.05-1.7 Myr [inference; no sourced Sicista generation time or mutation rate exists].
- At the placeholder rate a 0.8 Myr split gives d about 0.9%, so about 58% of 100-bp bins are heterozygous [my calculation]; confirm `-b10`.
- Unquantified for this dataset: ONT errors, ROH (about 30% of betulina), possible WGA in trizona (a project page mentions planned HiFi and no WGA, conflicting with the stated WGA), X chromosome, mapping both to one reference.

## Not verified

Full text of Cahill 2016; any hPSMC run on human populations; the origin of "reference choice does not matter"; Sicista mutation rate and generation time; whether trizona was WGA'd. See the report's last section.
