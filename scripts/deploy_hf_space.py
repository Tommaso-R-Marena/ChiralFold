#!/usr/bin/env python3
"""Deploy hf_space/ to https://huggingface.co/spaces/The-Philosopher/ChiralFold

Usage:
    export HF_TOKEN=hf_...   # write token from https://huggingface.co/settings/tokens
    python scripts/deploy_hf_space.py
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

from huggingface_hub import HfApi, create_repo

REPO_ID = "The-Philosopher/ChiralFold"
ROOT = Path(__file__).resolve().parents[1]
SPACE_DIR = ROOT / "hf_space"


def main() -> int:
    token = os.environ.get("HF_TOKEN") or os.environ.get("HUGGINGFACE_HUB_TOKEN")
    if not token:
        print("Set HF_TOKEN (Hugging Face write token) and re-run.", file=sys.stderr)
        return 1
    if not SPACE_DIR.is_dir():
        print(f"Missing {SPACE_DIR}", file=sys.stderr)
        return 1

    api = HfApi(token=token)
    who = api.whoami()
    print(f"Authenticated as: {who.get('name')}")
    create_repo(
        repo_id=REPO_ID,
        repo_type="space",
        space_sdk="gradio",
        exist_ok=True,
        token=token,
    )
    api.upload_folder(
        folder_path=str(SPACE_DIR),
        repo_id=REPO_ID,
        repo_type="space",
        commit_message="Deploy ChiralFold web Space",
        ignore_patterns=["__pycache__", "*.pyc", ".git*"],
    )
    print(f"Live at https://huggingface.co/spaces/{REPO_ID}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
