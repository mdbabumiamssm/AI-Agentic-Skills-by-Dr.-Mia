# Codex Repo Refresh

This repository now includes a runnable automation wrapper for Codex so you can repeat the same workflow used in the manual provider/framework refresh.

## What it does

The wrapper:

1. fetches the latest `origin`
2. fast-forwards the chosen base branch
3. creates a fresh working branch
4. runs `codex exec` non-interactively with a repository-specific prompt
5. validates changed `SKILL.md` frontmatter
6. commits the result
7. optionally pushes the branch or updates `main`

## Files

- `scripts/run_codex_repo_refresh.sh`
- `automation/prompts/refresh_provider_framework_skills.prompt.md`
- `automation/validate_skill_frontmatter.py`

## Prerequisites

- `codex` CLI installed and authenticated
- `git` and `python3` available
- write access to the GitHub remote
- optional: `gh` if you want the script to open a pull request automatically

## Recommended usage

Create a branch, run Codex, commit, and push the branch:

```bash
bash scripts/run_codex_repo_refresh.sh
```

Run Codex, commit locally, but do not push:

```bash
bash scripts/run_codex_repo_refresh.sh --publish none
```

Run Codex and update `main` directly if the result can be fast-forwarded cleanly:

```bash
bash scripts/run_codex_repo_refresh.sh --publish main
```

Use a specific model:

```bash
bash scripts/run_codex_repo_refresh.sh --model gpt-5.4
```

Push a branch and open a PR with GitHub CLI:

```bash
bash scripts/run_codex_repo_refresh.sh --open-pr
```

## Safety defaults

- refuses to start from a dirty worktree unless `--allow-dirty` is set
- does not update `main` by default
- uses a fresh timestamped branch
- refuses to fast-forward `main` if `origin/main` moved during the run

## Customization

Edit the prompt template if you want Codex to refresh a different area of the repository:

- `automation/prompts/refresh_provider_framework_skills.prompt.md`

The current prompt is tuned for:

- `Skills/AI_Providers`
- `Skills/Agentic_AI`
- `README.md`

## Notes

- The validator only checks changed `SKILL.md` files for required frontmatter fields.
- The wrapper handles git publication. The Codex prompt explicitly tells the agent not to push or merge on its own.
