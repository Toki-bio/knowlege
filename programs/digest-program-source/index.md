---
title: "digest-program-source — digest"
type: reference
topic: Program digests
generated: true
---

> Type: reference · Claude Code skill `digest-program-source` · [source on GitHub](https://github.com/Toki-bio/knowlege/tree/main/.claude/skills/digest-program-source)

**When to use:** Deep-read a third-party program's source code (a GitHub repo, a local clone, or a package) and write a reusable digest of it as a new skill in the Toki-bio/knowlege repo. Use when asked to "learn", "digest", "study", or "deeply read" a program's source so it can be reused in later sessions, or to update an existing program digest after upstream changed.

---

{% raw %}
# Digest a program's source into a reusable skill

The goal is a digest that a future Claude session, with no memory of this one,
can read in a minute and then work on the program as if it had already spent
an hour in the source. Write for that reader.

## 0. Check for an existing digest first

Look in `.claude/skills/` of the knowledge repo (Toki-bio/knowlege) for a
folder named after the program. If one exists, read its `SKILL.md`, compare the
pinned commit with upstream, and **update** rather than start over. Record what
changed in the "Change log" section.

## 1. Pin the source

- Clone or locate the source. Record: repo URL, tag/version, commit SHA, date,
  language(s), build system, license.
- Every `file:line` reference in the digest must be valid at that SHA. Never
  cite line numbers from a different checkout.

## 2. Read in this order

1. README, docs, CHANGELOG: what the authors say it does.
2. Build files and entry points (`main`, CLI parser, `setup.py`/`Cargo.toml`/`CMakeLists.txt`).
3. The core data structures. Most programs are understood through 3 to 6 structs or classes.
4. The main control flow for each top-level command or mode, as a numbered trace with `file:line`.
5. Tests: they show intended behavior and edge cases faster than the code does.
6. Only then, details: algorithms, I/O formats, parameters, performance tricks.

Use the Explore agent for broad sweeps and read key files yourself. Do not
paste large source blocks into the digest; cite `file:line` and summarize.

## 3. Write the digest

Copy `references/digest-template/` to
`<knowledge repo>/.claude/skills/<program-name>/` and fill it in. Rules:

- `SKILL.md` stays under 150 lines. It is the map, not the territory.
- The `description` frontmatter must name the program and say when to use the
  skill (for example "Use when reading, modifying, or debugging minimap2
  internals, or when its output or parameters need explaining"). This line is
  what triggers automatic reuse in later sessions, so make it specific.
- Bulk content goes in `references/*.md`, one file per concern:
  `architecture.md`, `data-structures.md`, `control-flow.md`, `io-formats.md`,
  `pitfalls.md`, `build-run-test.md`. Add or drop files as the program demands.
- Prefer tables and numbered traces over prose. Every claim about behavior
  should carry a `file:line` or a test name.
- Record open questions and things you could not verify, marked as such.
- Record "Type: notes" style content (your own interpretations) separately
  from verified facts, matching the repo's reference/explanation/notes split.

## 4. Commit and push to the knowledge repo

- Commit message: `Add <program> digest (upstream <short-sha>)` or
  `Update <program> digest to <version>`.
- Push to the branch the session was given. Do not open a PR unless asked.
- If this session is not running inside the knowledge repo, clone it to a
  scratch directory, add the skill there, commit, and push.

## 5. Report back

Tell the user: skill path, pinned SHA, what the digest covers, what it does
not, and open questions. One short list.
{% endraw %}
