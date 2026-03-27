You are running non-interactively inside the repository at `{{REPO_ROOT}}` on `{{TODAY}}`.

Your task is to refresh this repository's current AI provider and agent framework skill coverage in the same spirit as the manual work previously completed in this repo.

Requirements:

1. Inspect the current repository state before editing.
2. Focus on:
   - `README.md`
   - `Skills/AI_Providers`
   - `Skills/Agentic_AI`
3. Coordinate with both current and older repository content:
   - update or extend existing first-party skills instead of creating duplicates
   - preserve upstream changes already present on `{{BASE_BRANCH}}`
   - avoid touching `Skills/External_Collections` unless absolutely necessary
4. Use live web search when needed, but only rely on official documentation and official GitHub repositories for current provider/framework facts.
5. Preserve local repository conventions:
   - copyright header blocks
   - lean `SKILL.md` files with frontmatter fields:
     - `name`
     - `description`
     - `measurable_outcome`
     - `allowed-tools`
   - supporting `references/` files
   - `agents/openai.yaml`
6. Prioritize first-party AI platform and agent framework coverage.
7. Add or refresh skills only where the update is justified by current official sources.
8. After editing, run:
   - `python3 automation/validate_skill_frontmatter.py` for any changed `SKILL.md` files
9. Do not push branches, open PRs, or merge to `main`. The wrapper script handles git publication after you finish.

When done, provide a concise summary of:

- files changed
- why those changes were made
- any follow-up areas still worth refreshing later
