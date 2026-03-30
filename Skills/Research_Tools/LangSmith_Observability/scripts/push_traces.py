"""Replay stored trace JSON into LangSmith projects."""

import argparse
import json
from pathlib import Path

from langsmith import Client


def push_trace(trace_path: Path, project: str):
    client = Client()
    payload = json.loads(trace_path.read_text())
    client.runs.create(**payload, project_name=project)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("trace", type=Path)
    parser.add_argument("--project", default="Clinical-Ops-Agent")
    args = parser.parse_args()
    push_trace(args.trace, args.project)


if __name__ == "__main__":
    main()
