> Type: explanation

# Twin-Strand Duplex Sequencing Explained

## Summary

**Twin-strand duplex sequencing** (commonly called **Duplex Sequencing**) is an ultra-high-accuracy DNA sequencing method that distinguishes **true biological mutations** from **technical errors** by independently sequencing **both original strands** of a DNA molecule.

The central idea is simple:

> A real mutation exists in both strands of the original DNA duplex. Technical errors usually appear on only one strand after library preparation or PCR.

This reduces sequencing error rates from approximately **10⁻³** to **10⁻⁷–10⁻⁸ errors per base**, enabling reliable detection of extremely rare variants.

---

# Why Standard Sequencing Produces False Mutations

DNA sequencing involves several error-prone steps:

* Library preparation
* PCR amplification
* Sequencing chemistry
* Base calling

Each of these can introduce incorrect nucleotides.

For example, the original molecule may be

```
A G C T G A
```

but the sequencer reports

```
A A C T G A
```

If only one observation exists, there is no way to determine whether this is:

* a genuine mutation, or
* a sequencing artifact.

---

# The Key Insight

DNA is naturally double-stranded.

```
Strand A:  A G C T G A
           | | | | | |
Strand B:  T C G A C T
```

These strands contain identical genetic information.

A **real mutation** is present in both strands.

A **technical error** introduced after the strands are separated usually affects only one strand.

Duplex sequencing exploits this property.

---

# How Duplex Sequencing Works

## Step 1 — Fragment DNA

Genomic DNA is randomly sheared into fragments.

```
Adapter ─ DNA fragment ─ Adapter
```

---

## Step 2 — Attach Molecular Tags

Before PCR begins, adapters containing **Unique Molecular Identifiers (UMIs)** are ligated onto every DNA fragment.

```
[Adapter + UMI]──DNA Fragment──[Adapter + UMI]
```

These UMIs uniquely identify each original DNA molecule.

Importantly, **the adapters—not PCR primers—carry the tags.**

PCR primers later bind only to these adapter sequences.

---

## Step 3 — Sequence Both Original Strands

Each original strand receives its own molecular identity.

Conceptually:

```
Original duplex

Strand A  ← Tag X

Strand B  ← Tag Y
```

After sequencing, software reconstructs:

* consensus of reads from Strand A
* consensus of reads from Strand B

These are called **Single-Strand Consensus Sequences (SSCSs).**

---

## Step 4 — Compare Both Strands

Only mutations confirmed by **both** original strands are accepted.

Example:

```
Strand A consensus:

A G C T

Strand B consensus:

T C G A
```

If the mutation is complementary on both strands, it is considered genuine.

Otherwise, it is discarded.

---

# Why Sequencing Errors Are Removed

Suppose the true DNA is

```
A A C T
```

During sequencing, one read is incorrectly called:

```
Observed:

A G C T
```

The opposite strand is sequenced correctly.

Result:

```
Strand A consensus:
A G C T

Strand B consensus:
T T G A
```

The two strands no longer agree.

The mutation is rejected.

---

# What About PCR Errors?

A common question is:

> PCR creates complementary strands. Doesn't the error simply get copied onto the opposite strand?

No.

The important distinction is between:

* **original DNA strands**
* **PCR-generated copies**

PCR copies are **not independent evidence**.

Instead, every read is traced back to one of the two **original strands** using the molecular tags.

Suppose PCR makes an error while copying Strand A:

```
Original Strand A

A G C T

↓

PCR error

A A C T
```

Now every descendant of this PCR product carries the error.

However, the descendants of the original **Strand B** remain correct.

At the end:

```
Consensus of Strand A:
contains error

Consensus of Strand B:
does not
```

Since the two original strands disagree, the mutation is discarded.

---

# Why Matching Errors Are Extremely Unlikely

For a technical error to fool duplex sequencing,

the same incorrect mutation would need to occur independently on **both original strands**.

For example:

```
Strand A

A → G

and simultaneously

Strand B

T → C
```

These are independent events.

If the per-strand error rate is roughly

```
10⁻³
```

then identical complementary errors occur with probability

```
10⁻³ × 10⁻³ = 10⁻⁶
```

After consensus building and quality filtering, practical error rates reach approximately

```
10⁻⁷–10⁻⁸ per base
```

---

# How Are UMIs Added?

A common misconception is that the molecular tags come from PCR primers.

They do not.

The workflow is:

1. Fragment genomic DNA.
2. Repair fragment ends.
3. Ligate adapters containing random UMIs.
4. Perform PCR using primers that bind the adapters.

Diagram:

```
          Adapter
        +---------+
        |  UMI    |
        +---------+

           │
           ▼

[UMI]──DNA Fragment──[UMI]

        │
        ▼

PCR primers bind here
```

Because the UMIs are attached **before amplification**, every PCR product inherits the identity of the original DNA molecule.

This allows software to determine:

* which reads originated from the same DNA fragment,
* which reads came from Strand A,
* which reads came from Strand B.

---

# Key Takeaways

* Duplex sequencing independently analyzes both original DNA strands.
* True mutations exist in both strands.
* PCR and sequencing errors almost always occur after strand separation and therefore affect only one original strand.
* UMIs are attached through adapter ligation before PCR, not by genomic primers.
* Comparing consensus sequences from both original strands reduces sequencing error rates by several orders of magnitude, enabling reliable detection of extremely rare mutations.
