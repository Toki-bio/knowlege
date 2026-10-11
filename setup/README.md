# Setup: using this repo's skills in every Claude Code session

The program digests live in `.claude/skills/<program>/`. This page explains how
to make them show up in sessions outside this repo.

## Local: VS Code extension or CLI (one-time)

The VS Code extension and the CLI share `~/.claude`, so a plugin installed at
user scope in one is available in the other.

**Option A, in the VS Code panel or CLI prompt:** type `/plugins`, open the
Marketplaces tab, add `Toki-bio/knowlege`, then install the `knowledge` plugin.

**Option B, one line in the Claude prompt (CLI or VS Code panel):**

```
/plugin install knowledge --marketplace Toki-bio/knowlege
```

**Option B2, from any terminal (no session needed):**

```sh
claude plugin marketplace add Toki-bio/knowlege
claude plugin install knowledge@toki-knowledge --scope user
```

**Option C, one click in VS Code:**

```
vscode://anthropic.claude-code/install-plugin?plugin=knowledge&marketplace=Toki-bio/knowlege
```

Then turn on auto-update for the marketplace (`/plugins` > Marketplaces >
Enable auto-update), or run `claude plugin update knowledge` after you push
new digests. Updates are detected by the plugin `version` in
`.claude-plugin/marketplace.json`, so bump it when adding a digest.

Also paste `claude-md-snippet.md` into your user-level instructions
(VS Code: Claude Code panel > Customize > Instructions, which edits
`~/.claude/CLAUDE.md`). That is the "write-back" rule: skills handle the
reading side automatically through their descriptions, the rule handles the
writing side.

## Cloud sessions (claude.ai/code)

Cloud containers do not see your local `~/.claude`. Paste `cloud-setup.sh`
into the environment's setup script. It clones this repo, copies the skills
into `~/.claude/skills/`, and appends the snippet to `~/.claude/CLAUDE.md`
inside the container. Sessions started *in* this repo need nothing: project
skills under `.claude/skills/` load automatically.

## Writing a new digest

In any session: `/digest-program-source https://github.com/owner/repo`. The
skill reads the source, fills the template, and commits a new skill folder to
this repo. See `.claude/skills/digest-program-source/SKILL.md` for the format.

## Caveat

When the plugin is installed *and* you work inside this repo, the same skills
load twice (once as project skills, once as `knowledge:<name>`). Harmless,
just noisy in `/skills`.
