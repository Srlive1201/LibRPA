#!/usr/bin/env python3
from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence

try:
    import yaml
except ImportError as exc:  # pragma: no cover - exercised only on missing dependency
    raise SystemExit("PyYAML is required; install docs/requirements.txt") from exc


@dataclass
class OptionDoc:
    name: str
    raw_type: str
    line: int
    description: str
    details: str
    default: str
    status: str
    since: str
    deprecated: str


def clean_doc_line(line: str) -> str:
    text = line.strip()
    for prefix in ("//!", "///", "/*!", "/**", "*/"):
        if text.startswith(prefix):
            text = text[len(prefix):].strip()
    if text.startswith("*"):
        text = text[1:].strip()
    return text


def iter_doc_lines(lines: Sequence[str], start: int) -> tuple[List[str], int]:
    text = lines[start].lstrip()
    if text.startswith(("//!", "///")):
        out: List[str] = []
        i = start
        while i < len(lines) and lines[i].lstrip().startswith(("//!", "///")):
            out.append(clean_doc_line(lines[i]))
            i += 1
        return out, i

    out = [clean_doc_line(lines[start])]
    i = start + 1
    while i < len(lines):
        out.append(clean_doc_line(lines[i]))
        if "*/" in lines[i]:
            return out, i + 1
        i += 1
    return out, i


def is_field_declaration(line: str) -> bool:
    text = line.strip()
    return bool(text.endswith(";") and "(" not in text and ")" not in text and not text.startswith("static "))


def parse_field(line: str) -> tuple[str, str] | None:
    # Field documentation also covers driver members with in-class defaults.
    declaration = line.strip().removesuffix(";").split("=", 1)[0].strip()
    match = re.match(
        r"(?P<type>(?:std::)?[A-Za-z_][A-Za-z0-9_:<>]*(?:\s+[A-Za-z_][A-Za-z0-9_:<>]*)*)"
        r"\s+(?P<name>[A-Za-z_][A-Za-z0-9_]*)(?:\s*\[[^\]]+\])?$",
        declaration,
    )
    if not match:
        return None
    return match.group("name"), match.group("type").strip()


def parse_doc(lines: Iterable[str]) -> tuple[str, str, Dict[str, str]]:
    description: List[str] = []
    details: List[str] = []
    meta: Dict[str, List[str]] = {}
    current = "description"

    for raw in lines:
        line = raw.strip()
        if line.startswith("@brief "):
            description.append(line.removeprefix("@brief ").strip())
            current = "details"
            continue
        if line.startswith("@par "):
            rest = line.removeprefix("@par ").strip()
            label, _, value = rest.partition(" ")
            current = label.lower()
            meta.setdefault(current, [])
            if value:
                meta[current].append(value.strip())
            continue
        if line.startswith("@since "):
            meta["since"] = [line.removeprefix("@since ").strip()]
            current = "details"
            continue
        if line.startswith("@deprecated"):
            meta["deprecated"] = [line.removeprefix("@deprecated").strip()]
            current = "details"
            continue
        if re.match(r"@[A-Za-z]+", line):
            details.append(line)
            current = "details"
            continue

        if current == "description":
            if line:
                description.append(line)
                current = "details"
            continue
        if current in {"default", "status", "since", "deprecated"}:
            if line:
                meta.setdefault(current, []).append(line)
            continue
        if line:
            details.append(line)

    meta_text = {key: normalize_markdown(" ".join(value).strip()) for key, value in meta.items()}
    return (
        normalize_markdown(" ".join(description).strip()),
        normalize_markdown("\n".join(details).strip()),
        meta_text,
    )


def normalize_markdown(text: str) -> str:
    return (
        text.replace(r"\f$", "$")
        .replace(r"\f[", "$$")
        .replace(r"\f]", "$$")
        .replace(r"\c ", "`")
    )


def parse_struct(path: Path, struct_name: str) -> Dict[str, OptionDoc]:
    lines = path.read_text(encoding="utf-8").splitlines()
    in_struct = False
    pending_doc: List[str] = []
    options: Dict[str, OptionDoc] = {}
    i = 0

    while i < len(lines):
        text = lines[i].strip()
        if not in_struct:
            if re.search(rf"\bstruct\s+{re.escape(struct_name)}\b", text) or re.search(r"\btypedef\s+struct\b", text):
                in_struct = True
            i += 1
            continue

        if re.match(rf"}}\s*{re.escape(struct_name)}\s*;|}};", text):
            break

        if text.startswith(("//!", "///", "/*!", "/**")):
            doc, i = iter_doc_lines(lines, i)
            pending_doc.extend(doc)
            continue

        if is_field_declaration(text):
            parsed = parse_field(text)
            if parsed:
                name, raw_type = parsed
                description, details, meta = parse_doc(pending_doc)
                options[name] = OptionDoc(
                    name=name,
                    raw_type=raw_type,
                    line=i + 1,
                    description=description,
                    details=details,
                    default=meta.get("default", ""),
                    status=meta.get("status", ""),
                    since=meta.get("since", ""),
                    deprecated=meta.get("deprecated", ""),
                )
            pending_doc = []
            i += 1
            continue

        if text and not text.startswith(("/*", "*", "//")):
            pending_doc = []
        i += 1

    return options


def display_type(raw_type: str) -> str:
    if raw_type == "char" or raw_type == "std::string":
        return "string"
    if raw_type in {"bool", "LibrpaSwitch"}:
        return "bool"
    if raw_type in {"double", "int", "long long"}:
        return raw_type
    if raw_type.startswith("Librpa"):
        return "enum/string"
    return raw_type


def table_cell(text: str) -> str:
    return text.replace("\n", "<br>").replace("|", "\\|")


def code_parameter_names(text: str, parameter_names: set[str]) -> str:
    for name in sorted(parameter_names, key=len, reverse=True):
        text = re.sub(rf"(?<![`A-Za-z0-9_]){re.escape(name)}(?![`A-Za-z0-9_])", f"`{name}`", text)
    return text


def status_text(option: OptionDoc) -> str:
    status = option.status
    if option.deprecated:
        status = f"Deprecated: {option.deprecated}" if not status else f"{status}; Deprecated: {option.deprecated}"
    return status


def option_anchor(name: str) -> str:
    return "runtime-parameter-" + option_short_anchor(name)


def option_short_anchor(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", name.lower()).strip("-")


def option_link(option: OptionDoc) -> str:
    return f"[`{option.name}`](#{option_short_anchor(option.name)})"


def source_label(source: str) -> str:
    return {"api": "API", "driver": "Driver"}.get(source, source)


def option_ref(field: str) -> tuple[str, str]:
    if "." not in field:
        raise ValueError(f"field must be SOURCE.NAME: {field}")
    return tuple(field.split(".", 1))  # type: ignore[return-value]


def get_option(field: str, sources: Dict[str, Dict[str, OptionDoc]]) -> OptionDoc:
    source, name = option_ref(field)
    try:
        return sources[source][name]
    except KeyError as exc:
        raise KeyError(f"unknown field {field}") from exc


def render_table(fields: Sequence[str], sources: Dict[str, Dict[str, OptionDoc]], parameter_names: set[str]) -> str:
    rows = [
        "| Parameter Name | Description | Default Value | Since |",
        "|----------------|-------------|---------------|-------|",
    ]

    for field in sorted(fields, key=lambda field: (option_ref(field)[1].lower(), option_ref(field)[0])):
        option = get_option(field, sources)

        rows.append(
            "| "
            + " | ".join(
                table_cell(value)
                for value in (
                    option_link(option),
                    code_parameter_names(option.description, parameter_names),
                    option.default,
                    option.since,
                )
            )
            + " |"
        )

    return "\n".join(rows)


def option_detail_lines(source: str, option: OptionDoc, parameter_names: set[str]) -> List[str]:
    lines = [
        f"({option_anchor(option.name)})=",
        f"({option_short_anchor(option.name)})=",
        "",
        f":::{{dropdown}} ({source_label(source)}) `{option.name}`",
        "",
    ]
    if option.description:
        lines.extend([code_parameter_names(option.description, parameter_names), ""])
    lines.extend(
        [
            f"- Type: `{display_type(option.raw_type)}`",
            f"- Default: {option.default or 'not documented'}",
        ]
    )
    if status_text(option):
        lines.append(f"- Status: {code_parameter_names(status_text(option), parameter_names)}")
    if option.since:
        lines.append(f"- Since: {option.since}")
    if option.details:
        lines.extend(["", code_parameter_names(option.details, parameter_names)])
    lines.append(":::")
    return lines


def render_details(fields: Sequence[str], sources: Dict[str, Dict[str, OptionDoc]]) -> str:
    seen: set[tuple[str, str]] = set()
    parameter_names = {name for options in sources.values() for name in options}
    sorted_fields = sorted(fields, key=lambda field: (option_ref(field)[1].lower(), option_ref(field)[0]))
    parts = ["## Parameter Details"]
    for field in sorted_fields:
        source, _ = option_ref(field)
        option = get_option(field, sources)
        key = (source, option.name)
        if key in seen:
            continue
        seen.add(key)
        parts.extend(["", *option_detail_lines(source, option, parameter_names)])
    return "\n".join(parts)


def render(config: dict, sources: Dict[str, Dict[str, OptionDoc]]) -> str:
    out = [
        "<!-- Generated by docs/user_guide/generate_runtime_parameters.py; edit runtime_parameters.yml or header docstrings instead. -->",
        "",
        f"# {config['title']}",
        "",
    ]
    if config.get("intro"):
        out.extend([config["intro"].rstrip(), ""])

    out.extend(["## Overview", ""])
    parameter_names = {name for options in sources.values() for name in options}
    detail_fields: List[str] = []
    for block in config["blocks"]:
        if "heading" in block:
            out.extend([f"### {block['heading']}", ""])
        if block.get("text"):
            out.extend([block["text"].rstrip(), ""])
        if block.get("fields"):
            out.extend([render_table(block["fields"], sources, parameter_names), ""])
            detail_fields.extend(block["fields"])
        if block.get("markdown"):
            out.extend([block["markdown"].rstrip(), ""])

    out.extend([render_details(detail_fields, sources), ""])
    return "\n".join(line for line in out if line is not None).rstrip() + "\n"


def load_sources(config: dict, config_path: Path) -> Dict[str, Dict[str, OptionDoc]]:
    root = config_path.parent
    sources: Dict[str, Dict[str, OptionDoc]] = {}
    for name, source in config["sources"].items():
        path = (root / source["path"]).resolve()
        sources[name] = parse_struct(path, source["struct"])
    return sources


def extract_function_body(text: str, name: str) -> str:
    match = re.search(rf"\b{name}\s*\([^)]*\)\s*\{{", text)
    if not match:
        return ""
    start = match.end()
    depth = 1
    i = start
    while i < len(text) and depth:
        depth += text[i] == "{"
        depth -= text[i] == "}"
        i += 1
    return text[start:i - 1]


def strip_quotes(value: str) -> str:
    value = value.strip()
    if len(value) >= 2 and value[0] == value[-1] == '"':
        return value[1:-1]
    return value


def normalize_default(value: str) -> Optional[str]:
    value = value.strip().strip("`")
    low = value.lower()
    if low.startswith(("required", "build-dependent")):
        return None
    if " (" in value:
        value = value.split(" (", 1)[0]
    value = value.strip().strip("`")
    if value.lower() == "empty":
        return ""
    value = {
        "LIBRPA_SWITCH_ON": "true",
        "LIBRPA_SWITCH_OFF": "false",
        "LIBRPA_VERBOSE_INFO": "info",
    }.get(value, value)
    return re.sub(r"\s+", " ", value)


def defaults_equal(left: str, right: str) -> bool:
    if left == right:
        return True
    try:
        return float(left) == float(right)
    except ValueError:
        return left.rstrip("/") == right.rstrip("/")


def parse_api_defaults(path: Path) -> Dict[str, Optional[str]]:
    body = extract_function_body(path.read_text(encoding="utf-8"), "librpa_init_options")
    values: Dict[str, List[str]] = {}
    for match in re.finditer(r"\blibrpa_set_([A-Za-z_][A-Za-z0-9_]*)\s*\(\s*opts\s*,\s*([^)]+?)\s*\)", body):
        field, raw = match.groups()
        value = strip_quotes(raw)
        if field == "output_dir" and value and not value.endswith("/"):
            value += "/"
        values.setdefault(field, []).append(value)
    for match in re.finditer(r"\bopts\s*->\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*([^;]+);", body):
        field, raw = match.groups()
        value = normalize_default(strip_quotes(raw))
        if value is not None:
            values.setdefault(field, []).append(value)

    defaults: Dict[str, Optional[str]] = {}
    for field, field_values in values.items():
        unique = sorted(set(field_values))
        defaults[field] = unique[0] if len(unique) == 1 else None
    return defaults


def parse_driver_defaults(path: Path, header_path: Path) -> Dict[str, str]:
    text = path.read_text(encoding="utf-8")
    header = header_path.read_text(encoding="utf-8")
    constants = dict(re.findall(r"\bstatic\s+constexpr\s+int\s+(\w+)\s*=\s*([0-9]+)\s*;", header))
    match = re.search(r"DriverParams::DriverParams\(\)\s*:\s*(.*?)\n\{", text, re.S)
    if not match:
        return {}
    defaults: Dict[str, str] = {}
    for field, raw in re.findall(r"\b([A-Za-z_][A-Za-z0-9_]*)\s*\(([^()]*)\)\s*,?", match.group(1)):
        value = strip_quotes(raw)
        value = constants.get(value, value)
        normalized = normalize_default(value)
        if normalized is not None:
            defaults[field] = normalized
    return defaults


def check_default_consistency(config: dict, config_path: Path, sources: Dict[str, Dict[str, OptionDoc]]) -> List[str]:
    root = config_path.parent
    api_defaults = parse_api_defaults((root / "../../src/api/options.cpp").resolve())
    driver_header = (root / config["sources"]["driver"]["path"]).resolve()
    driver_defaults = parse_driver_defaults((root / "../../driver/driver.cpp").resolve(), driver_header)
    impl = {"api": api_defaults, "driver": driver_defaults}
    issues: List[str] = []

    for block in config["blocks"]:
        for field in block.get("fields", []):
            source, name = option_ref(field)
            doc_default = normalize_default(sources[source][name].default)
            if doc_default is None:
                continue
            impl_default = impl[source].get(name)
            if impl_default is None:
                continue
            if not defaults_equal(doc_default, impl_default):
                issues.append(f"{field}: header default `{sources[source][name].default}` != implementation `{impl_default}`")
    return issues


def parse_args(argv: Sequence[str]) -> argparse.Namespace:
    here = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Generate runtime_parameters.md from LibRPA option headers.")
    parser.add_argument("--config", type=Path, default=here / "runtime_parameters.yml")
    parser.add_argument("--output", type=Path, default=here / "runtime_parameters.md")
    parser.add_argument("--check", action="store_true", help="fail if output is not up to date")
    parser.add_argument("--check-defaults", action="store_true", help="compare @par Default values with implementation defaults")
    return parser.parse_args(argv)


def main(argv: Sequence[str]) -> int:
    args = parse_args(argv)
    config = yaml.safe_load(args.config.read_text(encoding="utf-8"))
    sources = load_sources(config, args.config)
    rendered = render(config, sources)

    if args.check_defaults:
        issues = check_default_consistency(config, args.config, sources)
        if issues:
            print("Runtime parameter default mismatch(es):", file=sys.stderr)
            for issue in issues:
                print(f"  - {issue}", file=sys.stderr)
            return 1

    if args.check:
        old = args.output.read_text(encoding="utf-8") if args.output.exists() else ""
        if old != rendered:
            print(f"{args.output} is out of date; run {Path(__file__).name}", file=sys.stderr)
            return 1
        return 0

    args.output.write_text(rendered, encoding="utf-8")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
