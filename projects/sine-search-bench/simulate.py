#!/usr/bin/env python3
"""Plant diverged copies of a SINE consensus into a background sequence.

Writes genome.fa (one or more contigs), truth.bed (planted copies with their
realised divergence and length) and query.fa (the consensus used).

Divergence model per copy: target identity drawn from a grid or uniform range;
substitutions with 2:1 transition bias; indels at a fraction of the
substitution count (geometric lengths); optional 5' truncation; variable
poly-A tail. Ground truth records the genomic span of the planted copy.
"""
import argparse, random, sys

COMP = str.maketrans("ACGT", "TGCA")
TRANSITION = {"A": "G", "G": "A", "C": "T", "T": "C"}


def read_fasta(path):
    name, seq = None, []
    for line in open(path):
        if line.startswith(">"):
            if name:
                yield name, "".join(seq)
            name, seq = line[1:].split()[0], []
        else:
            seq.append(line.strip().upper())
    if name:
        yield name, "".join(seq)


def markov_background(n, rng, gc=0.42):
    """Order-0 background with given GC; crude but no accidental repeats."""
    at = (1 - gc) / 2
    g = gc / 2
    return "".join(rng.choices("ACGT", weights=[at, g, g, at], k=n))


def mutate(seq, divergence, rng, indel_frac=0.12):
    """Return (mutated_seq, n_sub, n_indel_events). divergence = fraction of
    positions substituted; indel events = indel_frac * n_sub."""
    s = list(seq)
    L = len(s)
    n_sub = round(divergence * L)
    for i in rng.sample(range(L), min(n_sub, L)):
        b = s[i]
        if b not in "ACGT":
            continue
        if rng.random() < 2 / 3:
            s[i] = TRANSITION[b]
        else:
            s[i] = rng.choice([x for x in "ACGT" if x != b and x != TRANSITION[b]])
    n_indel = round(indel_frac * n_sub)
    for _ in range(n_indel):
        pos = rng.randrange(len(s))
        ln = 1
        while rng.random() < 0.5 and ln < 8:
            ln += 1
        if rng.random() < 0.5:
            del s[pos:pos + ln]
        else:
            s[pos:pos] = rng.choices("ACGT", k=ln)
    return "".join(s), n_sub, n_indel


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--consensus", required=True, help="FASTA with one consensus")
    ap.add_argument("--bg-length", type=int, default=20_000_000)
    ap.add_argument("--contigs", type=int, default=4)
    ap.add_argument("--background", default=None, help="real FASTA to use as background instead of random; "
                    "its first --contigs records are used, each cut to --bg-length/--contigs bp (0 = whole)")
    ap.add_argument("--copies", type=int, default=2000)
    ap.add_argument("--div-min", type=float, default=0.0)
    ap.add_argument("--div-max", type=float, default=0.40)
    ap.add_argument("--trunc5-frac", type=float, default=0.2, help="fraction of copies 5'-truncated")
    ap.add_argument("--tail-max", type=int, default=30, help="max extra poly-A tail length")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default=".")
    a = ap.parse_args()
    rng = random.Random(a.seed)

    (qname, cons), = list(read_fasta(a.consensus))
    cons = cons.replace("N", "A")
    per_contig = a.bg_length // a.contigs
    if a.background:
        contigs, names = [], []
        for name, seq in read_fasta(a.background):
            contigs.append(list(seq if per_contig == 0 else seq[:per_contig]))
            names.append(name)
            if len(contigs) == a.contigs:
                break
        per_contig = min(len(c) for c in contigs)
    else:
        contigs = [list(markov_background(per_contig, rng)) for _ in range(a.contigs)]
        names = [f"ctg{i+1}" for i in range(a.contigs)]

    # choose insertion sites: random non-N positions, at least 1 kb apart
    sites = []
    for ci in range(a.contigs):
        n = a.copies // a.contigs + (1 if ci < a.copies % a.contigs else 0)
        c = contigs[ci]
        chosen, tries = [], 0
        while len(chosen) < n and tries < 50 * n:
            tries += 1
            pos = rng.randrange(100, len(c) - 100)
            if "N" in c[pos - 50:pos + 50]:
                continue
            chosen.append(pos)
        chosen.sort()
        kept = []
        for pos in chosen:
            if not kept or pos - kept[-1] >= 1000:
                kept.append(pos)
        sites.extend((ci, pos) for pos in kept)
    # insert from the end of each contig backwards so coordinates stay valid
    sites.sort(key=lambda t: (t[0], -t[1]))
    truth = []
    for cid, (ci, pos) in enumerate(sites):
        div = rng.uniform(a.div_min, a.div_max)
        src = cons
        trunc = 0
        if rng.random() < a.trunc5_frac:
            trunc = rng.randrange(10, len(cons) // 2)
            src = cons[trunc:]
        src = src + "A" * rng.randrange(0, a.tail_max + 1)
        mut, nsub, nind = mutate(src, div, rng)
        strand = "+" if rng.random() < 0.5 else "-"
        if strand == "-":
            mut = mut.translate(COMP)[::-1]
        contigs[ci][pos:pos] = list(mut)
        truth.append((names[ci], pos, pos + len(mut), f"copy{cid}", div, strand, trunc, nsub, nind, len(mut)))

    with open(f"{a.out}/genome.fa", "w") as f:
        for ci, c in enumerate(contigs):
            f.write(f">{names[ci]}\n")
            s = "".join(c)
            for i in range(0, len(s), 80):
                f.write(s[i:i+80] + "\n")
    # truth coordinates: insertions done high-to-low per contig, so each
    # recorded pos is already correct for later (higher) insertions? No:
    # inserting at a lower position shifts higher ones. Recompute by scanning.
    truth.sort(key=lambda t: (t[0], t[1]))
    with open(f"{a.out}/truth.bed", "w") as f:
        f.write("#chrom\tstart\tend\tname\tdivergence\tstrand\ttrunc5\tnsub\tnindel\tlen\n")
        shift = {nm: 0 for nm in names}
        for t in truth:
            chrom, s, e, name, div, strand, trunc, nsub, nind, ln = t
            s2 = s + shift[chrom]
            f.write(f"{chrom}\t{s2}\t{s2+ln}\t{name}\t{div:.3f}\t{strand}\t{trunc}\t{nsub}\t{nind}\t{ln}\n")
            shift[chrom] += ln
    with open(f"{a.out}/query.fa", "w") as f:
        f.write(f">{qname}\n{cons}\n")
    print(f"wrote {a.out}/genome.fa ({sum(len(c) for c in contigs):,} bp), {len(truth)} copies", file=sys.stderr)


if __name__ == "__main__":
    main()
