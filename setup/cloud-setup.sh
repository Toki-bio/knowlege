#!/usr/bin/env bash
# Paste this into the "Setup script" of a claude.ai/code cloud environment
# (Settings > Environments > your environment). It runs before Claude starts,
# so every cloud session in that environment, in ANY repo, gets the knowledge
# skills and the write-back rule.
#
# Local machines (CLI or VS Code extension) do NOT need this: install the
# plugin once instead, see setup/README.md.
set -euo pipefail

KB_REPO="https://github.com/Toki-bio/knowlege.git"
KB_DIR="$HOME/.knowledge-repo"

# 1. Fetch the knowledge repo (shallow; it only needs the current content).
if [ -d "$KB_DIR/.git" ]; then
  git -C "$KB_DIR" pull --ff-only --quiet || true
else
  git clone --depth 1 --quiet "$KB_REPO" "$KB_DIR"
fi

# 2. Expose its skills as user-level skills for this container.
#    (Account-synced skills already live in ~/.claude/skills/, so this
#    directory is read by cloud sessions.)
mkdir -p "$HOME/.claude/skills"
for skill in "$KB_DIR"/.claude/skills/*/; do
  name="$(basename "$skill")"
  rm -rf "$HOME/.claude/skills/$name"
  cp -r "$skill" "$HOME/.claude/skills/$name"
done

# 3. Add the write-back rule to the user-level CLAUDE.md (idempotent).
mkdir -p "$HOME/.claude"
if ! grep -q "knowledge repo (Toki-bio/knowlege)" "$HOME/.claude/CLAUDE.md" 2>/dev/null; then
  cat "$KB_DIR/setup/claude-md-snippet.md" >> "$HOME/.claude/CLAUDE.md"
fi
