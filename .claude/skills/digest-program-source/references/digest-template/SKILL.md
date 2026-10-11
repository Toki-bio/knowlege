---
name: PROGRAM-NAME
description: Internals of PROGRAM-NAME (one line on what it does). Use when reading, modifying, debugging, or explaining PROGRAM-NAME source, parameters, or output.
---

# PROGRAM-NAME internals

> Type: reference (facts verified against the pinned commit) with marked `notes`.

## Pinned source

| Field | Value |
|---|---|
| Repo | https://github.com/OWNER/REPO |
| Version / tag | vX.Y.Z |
| Commit | `abcdef0` (full SHA in references/architecture.md) |
| Digested on | YYYY-MM-DD |
| Language / build | C, Makefile |
| License | MIT |

## What it does in three sentences

...

## Where to look (read on demand)

| Question | File |
|---|---|
| How is the code organized, what calls what | `references/architecture.md` |
| The 3 to 6 structs that everything revolves around | `references/data-structures.md` |
| Step-by-step trace of each command / mode | `references/control-flow.md` |
| Input and output formats, flags that change them | `references/io-formats.md` |
| Gotchas, non-obvious defaults, known bugs | `references/pitfalls.md` |
| Build, run, test, debug locally | `references/build-run-test.md` |

Read only the file that answers the current question. Start with
`architecture.md` when the task is unfamiliar.

## Ten facts worth knowing before touching the code

1. ...
2. ...

## Open questions / unverified

- ...

## Change log

- YYYY-MM-DD: initial digest at `abcdef0`.
