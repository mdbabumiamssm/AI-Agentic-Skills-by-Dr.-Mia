#!/usr/bin/env python3
from __future__ import annotations

import re
import sys
from pathlib import Path


REQUIRED_FIELDS = ["name", "description", "measurable_outcome", "allowed-tools"]


def parse_frontmatter(text: str) -> dict[str, object] | None:
    match = re.search(r"(?ms)^---\s*\n(.*?)\n---\s*\n", text)
    if not match:
        return None

    yaml_text = match.group(1)
    metadata: dict[str, object] = {}
    current_key: str | None = None

    for raw_line in yaml_text.splitlines():
        line = raw_line.rstrip()
        if not line or line.lstrip().startswith("#"):
            continue

        if line.strip().startswith("- "):
            if current_key is None:
                continue
            value = line.strip()[2:].strip()
            metadata.setdefault(current_key, [])
            assert isinstance(metadata[current_key], list)
            metadata[current_key].append(value)
            continue

        if ":" not in line:
            continue

        key, value = line.split(":", 1)
        key = key.strip()
        value = value.strip()
        current_key = key

        if not value:
            metadata[key] = []
            continue

        if len(value) >= 2 and (
            (value.startswith('"') and value.endswith('"'))
            or (value.startswith("'") and value.endswith("'"))
        ):
            value = value[1:-1]

        metadata[key] = value

    return metadata


def iter_skill_files(inputs: list[str]) -> list[Path]:
    paths: list[Path] = []
    for item in inputs:
        path = Path(item)
        if not path.exists():
            raise FileNotFoundError(item)
        if path.is_file():
            paths.append(path)
            continue
        paths.extend(sorted(path.rglob("SKILL.md")))
    return paths


def validate_skill_file(path: Path) -> list[str]:
    errors: list[str] = []
    text = path.read_text(encoding="utf-8")
    metadata = parse_frontmatter(text)
    if metadata is None:
        return [f"{path}: missing frontmatter"]

    for field in REQUIRED_FIELDS:
        if field not in metadata:
            errors.append(f"{path}: missing required field '{field}'")

    allowed_tools = metadata.get("allowed-tools")
    if allowed_tools is not None and not isinstance(allowed_tools, list):
        errors.append(f"{path}: 'allowed-tools' must be a YAML list")

    if "references/sources.md" in text:
        reference_file = path.parent / "references" / "sources.md"
        if not reference_file.exists():
            errors.append(f"{path}: references/sources.md is referenced but missing")
        else:
            ref_text = reference_file.read_text(encoding="utf-8", errors="ignore")
            if "http" not in ref_text:
                errors.append(f"{reference_file}: expected at least one source URL")

    if "agents/openai.yaml" in text or "openai.yaml" in text:
        ui_file = path.parent / "agents" / "openai.yaml"
        if not ui_file.exists():
            errors.append(f"{path}: agents/openai.yaml is referenced but missing")

    return errors


def main() -> int:
    if len(sys.argv) < 2:
        print("usage: python3 automation/validate_skill_frontmatter.py <SKILL.md or directory> [...]", file=sys.stderr)
        return 2

    try:
        skill_files = iter_skill_files(sys.argv[1:])
    except FileNotFoundError as exc:
        print(f"error: path not found: {exc}", file=sys.stderr)
        return 2

    if not skill_files:
        print("no SKILL.md files found")
        return 0

    errors: list[str] = []
    for path in skill_files:
        errors.extend(validate_skill_file(path))

    if errors:
        print("frontmatter validation failed:")
        for error in errors:
            print(f"- {error}")
        return 1

    print(f"validated {len(skill_files)} SKILL.md file(s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
