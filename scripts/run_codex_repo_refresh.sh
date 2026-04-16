#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'HELP'
Usage:
  bash scripts/run_codex_repo_refresh.sh [options]

Purpose:
  Run Codex non-interactively to refresh this repository's current AI provider and
  framework skill coverage, validate changed SKILL.md files, perform a repo-wide
  agentic/provider skill audit, commit the result, and optionally publish it.

Options:
  --base-branch <name>     Base branch to refresh from. Default: main
  --branch <name>          Working branch name. Default: chore/codex-repo-refresh-YYYYMMDD-HHMMSS
  --publish <mode>         One of: none, branch, main. Default: branch
  --model <model>          Codex model override. Example: gpt-5.4
  --no-search              Disable Codex live web search
  --allow-dirty            Allow starting from a dirty worktree
  --open-pr                After pushing a branch, open a PR with gh if available
  --commit-message <msg>   Commit message to use
  --help                   Show this help

Examples:
  bash scripts/run_codex_repo_refresh.sh
  bash scripts/run_codex_repo_refresh.sh --publish main
  bash scripts/run_codex_repo_refresh.sh --model gpt-5.4 --open-pr
  bash scripts/run_codex_repo_refresh.sh --publish none
HELP
}

require_cmd() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "error: required command not found: $1" >&2
    exit 1
  fi
}

die() {
  echo "error: $*" >&2
  exit 1
}

repo_slug_from_origin() {
  local remote_url
  remote_url="$(git -C "${REPO_ROOT}" config --get remote.origin.url)"
  remote_url="${remote_url#https://github.com/}"
  remote_url="${remote_url#git@github.com:}"
  remote_url="${remote_url%.git}"
  printf '%s\n' "${remote_url}"
}

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
PROMPT_TEMPLATE="${REPO_ROOT}/automation/prompts/refresh_provider_framework_skills.prompt.md"
VALIDATOR="${REPO_ROOT}/automation/validate_skill_frontmatter.py"

BASE_BRANCH="main"
BRANCH_NAME="chore/codex-repo-refresh-$(date +%Y%m%d-%H%M%S)"
PUBLISH_MODE="branch"
MODEL=""
SEARCH_ENABLED=true
ALLOW_DIRTY=false
OPEN_PR=false
COMMIT_MESSAGE="chore: codex refresh provider and framework skills"

while [[ $# -gt 0 ]]; do
  case "$1" in
    --base-branch)
      BASE_BRANCH="${2:?missing value for --base-branch}"
      shift 2
      ;;
    --branch)
      BRANCH_NAME="${2:?missing value for --branch}"
      shift 2
      ;;
    --publish)
      PUBLISH_MODE="${2:?missing value for --publish}"
      shift 2
      ;;
    --model)
      MODEL="${2:?missing value for --model}"
      shift 2
      ;;
    --no-search)
      SEARCH_ENABLED=false
      shift
      ;;
    --allow-dirty)
      ALLOW_DIRTY=true
      shift
      ;;
    --open-pr)
      OPEN_PR=true
      shift
      ;;
    --commit-message)
      COMMIT_MESSAGE="${2:?missing value for --commit-message}"
      shift 2
      ;;
    --help|-h)
      usage
      exit 0
      ;;
    *)
      die "unknown argument: $1"
      ;;
  esac
done

case "${PUBLISH_MODE}" in
  none|branch|main) ;;
  *)
    die "--publish must be one of: none, branch, main"
    ;;
esac

require_cmd git
require_cmd codex
require_cmd python3

[[ -f "${PROMPT_TEMPLATE}" ]] || die "missing prompt template: ${PROMPT_TEMPLATE}"
[[ -f "${VALIDATOR}" ]] || die "missing validator: ${VALIDATOR}"

git -C "${REPO_ROOT}" rev-parse --is-inside-work-tree >/dev/null 2>&1 || die "repository root is not a git repository: ${REPO_ROOT}"
git -C "${REPO_ROOT}" remote get-url origin >/dev/null 2>&1 || die "git remote 'origin' is not configured"

if [[ "${ALLOW_DIRTY}" == false ]] && [[ -n "$(git -C "${REPO_ROOT}" status --porcelain)" ]]; then
  die "worktree is dirty. Commit or stash changes first, or rerun with --allow-dirty"
fi

if [[ "${OPEN_PR}" == true ]]; then
  require_cmd gh
fi

echo "==> Fetching latest remote state"
git -C "${REPO_ROOT}" fetch origin

echo "==> Syncing ${BASE_BRANCH} with origin/${BASE_BRANCH}"
git -C "${REPO_ROOT}" checkout "${BASE_BRANCH}"
git -C "${REPO_ROOT}" merge --ff-only "origin/${BASE_BRANCH}"

if git -C "${REPO_ROOT}" rev-parse --verify "${BRANCH_NAME}" >/dev/null 2>&1; then
  die "branch already exists: ${BRANCH_NAME}"
fi

echo "==> Creating working branch ${BRANCH_NAME}"
git -C "${REPO_ROOT}" checkout -b "${BRANCH_NAME}"

RUN_DIR="$(mktemp -d)"
trap 'rm -rf "${RUN_DIR}"' EXIT

PROMPT_FILE="${RUN_DIR}/codex_prompt.md"
TODAY="$(date +%Y-%m-%d)"

sed \
  -e "s#{{REPO_ROOT}}#${REPO_ROOT}#g" \
  -e "s#{{TODAY}}#${TODAY}#g" \
  -e "s#{{BASE_BRANCH}}#${BASE_BRANCH}#g" \
  "${PROMPT_TEMPLATE}" > "${PROMPT_FILE}"

echo "==> Running Codex"
EXEC_ARGS=(exec -C "${REPO_ROOT}" --full-auto -o "${RUN_DIR}/codex_last_message.txt")
if [[ -n "${MODEL}" ]]; then
  EXEC_ARGS+=(-m "${MODEL}")
fi

if [[ "${SEARCH_ENABLED}" == true ]]; then
  CODEX_CMD=(codex --search "${EXEC_ARGS[@]}")
else
  CODEX_CMD=(codex "${EXEC_ARGS[@]}")
fi

"${CODEX_CMD[@]}" - < "${PROMPT_FILE}"

if [[ -f "${RUN_DIR}/codex_last_message.txt" ]]; then
  echo
  echo "==> Codex final message"
  cat "${RUN_DIR}/codex_last_message.txt"
  echo
fi

STATUS_OUTPUT="$(git -C "${REPO_ROOT}" status --porcelain)"
if [[ -z "${STATUS_OUTPUT}" ]]; then
  echo "==> Codex produced no repository changes"
  exit 0
fi

mapfile -t CHANGED_SKILLS < <(
  git -C "${REPO_ROOT}" status --porcelain | awk '{print $2}' | awk '/SKILL\.md$/ {print}'
)

if (( ${#CHANGED_SKILLS[@]} > 0 )); then
  echo "==> Validating changed SKILL.md files"
  ABS_SKILLS=()
  for rel_path in "${CHANGED_SKILLS[@]}"; do
    ABS_SKILLS+=("${REPO_ROOT}/${rel_path}")
  done
  python3 "${VALIDATOR}" "${ABS_SKILLS[@]}"
fi

echo "==> Running repo-wide skill audit for Agentic_AI and AI_Providers"
python3 "${VALIDATOR}" \
  "${REPO_ROOT}/Skills/Agentic_AI" \
  "${REPO_ROOT}/Skills/AI_Providers"

echo "==> Staging and committing changes"
git -C "${REPO_ROOT}" add .
git -C "${REPO_ROOT}" commit -m "${COMMIT_MESSAGE}"

if [[ "${PUBLISH_MODE}" == "none" ]]; then
  echo "==> Done. Changes committed locally on ${BRANCH_NAME}"
  exit 0
fi

echo "==> Pushing ${BRANCH_NAME}"
git -C "${REPO_ROOT}" push -u origin "${BRANCH_NAME}"

if [[ "${OPEN_PR}" == true && "${PUBLISH_MODE}" == "branch" ]]; then
  echo "==> Opening pull request"
  gh -R "$(repo_slug_from_origin)" pr create \
    --base "${BASE_BRANCH}" \
    --head "${BRANCH_NAME}" \
    --title "Codex refresh provider and framework skills" \
    --body "Automated Codex refresh of provider/framework skill coverage."
fi

if [[ "${PUBLISH_MODE}" == "branch" ]]; then
  echo "==> Done. Branch pushed: ${BRANCH_NAME}"
  exit 0
fi

echo "==> Preparing fast-forward update of ${BASE_BRANCH}"
git -C "${REPO_ROOT}" fetch origin
COUNTS="$(git -C "${REPO_ROOT}" rev-list --left-right --count "origin/${BASE_BRANCH}...${BRANCH_NAME}")"
LEFT_COUNT="${COUNTS%%$'\t'*}"

if [[ "${LEFT_COUNT}" != "0" ]]; then
  die "origin/${BASE_BRANCH} moved after the Codex run. Review branch ${BRANCH_NAME} and merge manually."
fi

git -C "${REPO_ROOT}" checkout "${BASE_BRANCH}"
git -C "${REPO_ROOT}" merge --ff-only "${BRANCH_NAME}"
git -C "${REPO_ROOT}" push origin "${BASE_BRANCH}"
