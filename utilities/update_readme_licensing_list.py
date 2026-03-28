#!/usr/bin/env python3
"""
This script updates the version and licensing information of thirdparty components
from `docs/develop/licensing.md`, so that we do not need to maintain two copies.
"""
import argparse
from pathlib import Path
import re
import sys


def extract_section(text: str, heading: str) -> str:
    pattern = re.compile(
        rf"^(?:{re.escape(heading)})[ \t]*\n(.*?)(?=^#{{1,6}}[ \t]+|\Z)",
        re.MULTILINE | re.DOTALL,
    )
    m = pattern.search(text)
    if not m:
        raise ValueError(f"Cannot find section heading: {heading}")
    return m.group(1)


def extract_first_bullet_list(section_text: str) -> str:
    pattern = re.compile(
        r"(^- .*(?:\n(?!\n|#)(?:  .*|    .*|[-*+] .*)|\n)*)(?=\n(?:\n|[^ \t-]|\Z)|\Z)",
        re.MULTILINE,
    )
    m = pattern.search(section_text)
    if not m:
        raise ValueError("Cannot find bullet list in source section")
    return m.group(1).rstrip() + "\n"


def replace_first_bullet_list(section_text: str, new_list: str) -> str:
    pattern = re.compile(
        r"(^- .*(?:\n(?!\n|#)(?:  .*|    .*|[-*+] .*)|\n)*)(?=\n(?:\n|[^ \t-]|\Z)|\Z)",
        re.MULTILINE,
    )
    m = pattern.search(section_text)
    if not m:
        raise ValueError("Cannot find bullet list in destination section")
    return section_text[: m.start()] + new_list + section_text[m.end() :]


def replace_section(text: str, heading: str, new_section_body: str) -> str:
    pattern = re.compile(
        rf"^(?P<head>{re.escape(heading)}[ \t]*\n)(?P<body>.*?)(?=^#{{1,6}}[ \t]+|\Z)",
        re.MULTILINE | re.DOTALL,
    )
    m = pattern.search(text)
    if not m:
        raise ValueError(f"Cannot find section heading: {heading}")
    return text[: m.start("body")] + new_section_body + text[m.end("body") :]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-s", "--src", type=str, default=None, help="licensing.md")
    ap.add_argument("-d", "--dst", type=str, default=None, help="README.md")
    args = ap.parse_args()

    repo_root = Path(__file__).resolve().parents[1]
    src = repo_root / "docs" / "develop" / "licensing.md"
    dst = repo_root / "README.md"
    if args.src is not None:
        src = Path(args.src)
    if args.dst is not None:
        dst = Path(args.dst)

    src_text = src.read_text(encoding="utf-8")
    dst_text = dst.read_text(encoding="utf-8")

    src_section = extract_section(src_text, "# Licensing")
    dst_section = extract_section(dst_text, "## Licensing")

    new_list = extract_first_bullet_list(src_section)
    new_dst_section = replace_first_bullet_list(dst_section, new_list)

    updated_text = replace_section(dst_text, "## Licensing", new_dst_section)

    if updated_text != dst_text:
        dst.write_text(updated_text, encoding="utf-8")

    return 0


if __name__ == "__main__":
    sys.exit(main())
