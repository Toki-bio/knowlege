# Knowledge base

**Browse the site:** [https://toki-bio.github.io/knowlege/](https://toki-bio.github.io/knowlege/)

Personal knowledge base — durable, readable write-ups worth keeping, mostly distilled from AI chat discussions.

## Structure

Content lives under topic folders (`phylogenetics/`, `sequencing/`, `methods/`, `sines/`, `projects/`, `genomics/`, `population-genetics/`, etc.), one Markdown file per piece. New topic folders get added as they're needed — no fixed taxonomy up front.

Every file starts with YAML front matter (for the site) and a one-line type tag:

```markdown
---
title: Short title for navigation
type: explanation
---

> Type: reference | explanation | notes
```

- **reference** — short, factual, textbook-level: "what is X."
- **explanation** — longer write-ups that connect several concepts into a narrative (e.g. how a whole pipeline works end to end).
- **notes** — personal research notes, hypotheses, or interpretations. Not established fact — clearly the author's own thinking, kept separate so it's never mistaken for the other two categories.

The GitHub Pages site builds a sidebar from topic folders automatically. New topic folders need a matching entry under `defaults` in `_config.yml` (see the `phylogenetics` / `sequencing` examples) so notes get a `topic` and appear in the nav.

## Why this structure

Git history already tracks how an explanation improved over time, so there's no need to version files manually. Folder structure is kept shallow and is expected to be reorganized (`git mv`) as real topics emerge, rather than designed upfront for content that doesn't exist yet.
