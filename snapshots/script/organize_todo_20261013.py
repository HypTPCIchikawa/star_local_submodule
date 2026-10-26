#!/usr/bin/env python3
"""10/13 Todoだけを年月日階層へ移す、一回限りの機械的移行。

既定は計画表示のみ。--applyは元ファイルを消去せず移動し、Markdownの
ナビゲーションだけを修正する。ログ・設定snapshot・TXT/JSON/C等の内容は
変更しない。--checkは保存済みmanifestと現在のファイルを照合する。
"""

import argparse
import datetime
import hashlib
import json
import os
import pathlib
import re
import shutil


SOURCE = pathlib.PurePosixPath("Todo/20261013")
DESTINATION = pathlib.PurePosixPath("Todo/2026/10/13")
ROOT_TEXT = pathlib.PurePosixPath("Todo/STAR_Todo_20261013.txt")
MANIFEST = pathlib.PurePosixPath("Todo/2026/10/todo13_migration_20261020.json")
EXCLUDED_PARTS = {"validation", "source", "snapshots", ".git", "__pycache__"}
LINK = re.compile(r"(?P<start>!?\[[^\]\n]*\]\()(?P<target><[^>\n]+>|[^\s)\n]+)(?P<tail>(?:\s+[\"'][^\n]*?[\"'])?\))")
REFERENCE = re.compile(r"(?P<start>^\s*\[[^\]\n]+\]:\s*)(?P<target><[^>\n]+>|\S+)", re.MULTILINE)
INLINE_PATH = re.compile(r"`(?P<path>[^`\s]*Todo/(?:20261013(?:/[^`\s]*)?|STAR_Todo_20261013\.txt))`")


def digest(path):
    value = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            value.update(chunk)
    return value.hexdigest()


def relative(path, root):
    return path.relative_to(root).as_posix()


def map_path(path, root):
    source = root / SOURCE
    if path == root / ROOT_TEXT:
        return root / DESTINATION / ROOT_TEXT.name
    try:
        return root / DESTINATION / path.relative_to(source)
    except ValueError:
        return path


def raw_path_update(value):
    # User-facing absolute /star paths retain their original mount spelling.
    return value.replace("Todo/20261013", "Todo/2026/10/13").replace(
        "Todo/STAR_Todo_20261013.txt", "Todo/2026/10/13/STAR_Todo_20261013.txt")


def destination_target(value, old_file, new_file, root):
    angled = value.startswith("<") and value.endswith(">")
    raw = value[1:-1] if angled else value
    if not raw or raw.startswith("#") or re.match(r"^[A-Za-z][A-Za-z0-9+.-]*:", raw):
        return value, None
    # Query/fragment suffixes do not participate in path resolution.
    match = re.match(r"([^?#]*)(.*)", raw)
    path, suffix = match.groups()
    if not path:
        return value, None
    if path.startswith("/"):
        revised = raw_path_update(path) + suffix
        return ("<" + revised + ">" if angled else revised), None
    old_target = pathlib.Path(os.path.normpath(str(old_file.parent / path)))
    new_target = map_path(old_target, root)
    revised = pathlib.PurePath(os.path.relpath(new_target, new_file.parent)).as_posix() + suffix
    # Preserve an explicit ./ spelling when there is no relocation-induced change.
    if old_file == new_file and new_target == old_target:
        return value, None
    existence = {"oldTarget": str(old_target), "newTarget": str(new_target),
                 "oldTargetExisted": old_target.exists()}
    return ("<" + revised + ">" if angled else revised), existence


def rewrite_markdown(text, old_file, new_file, root):
    """Rewrite navigation only; fenced historical commands remain byte-for-byte."""
    checks = []
    output = []
    fence = None
    for line in text.splitlines(keepends=True):
        fence_match = re.match(r"^\s*(`{3,}|~{3,})", line)
        if fence_match:
            marker = fence_match.group(1)
            if fence is None:
                fence = marker[0]
            elif marker[0] == fence:
                fence = None
            output.append(line)
            continue
        if fence is not None:
            output.append(line)
            continue

        def update(match):
            target, check = destination_target(match.group("target"), old_file, new_file, root)
            if check:
                checks.append(check)
            return match.group("start") + target + (match.groupdict().get("tail") or "")

        line = LINK.sub(update, line)
        line = REFERENCE.sub(update, line)
        line = INLINE_PATH.sub(lambda match: "`" + raw_path_update(match.group("path")) + "`", line)
        output.append(line)
    return "".join(output), checks


def files_under(directory):
    result = []
    for base, directories, names in os.walk(str(directory), followlinks=False):
        for name in directories:
            if pathlib.Path(base, name).is_symlink():
                raise RuntimeError("Directory symlink requires explicit review: " + str(pathlib.Path(base, name)))
        for name in names:
            path = pathlib.Path(base, name)
            if path.is_symlink() or not path.is_file():
                raise RuntimeError("Non-regular migration input requires explicit review: " + str(path))
            result.append(path)
    return sorted(result)


def document_paths(root):
    paths = {root / "README.md", root / "docs/REFERENCE.md"}
    for directory in ("Todo", "mdfiles", "analysisnote"):
        for base, directories, names in os.walk(str(root / directory), followlinks=False):
            directories[:] = [name for name in directories if name not in EXCLUDED_PARTS]
            for name in names:
                path = pathlib.Path(base, name)
                if name.endswith(".md") and not path.is_symlink():
                    paths.add(path)
    return sorted(path for path in paths if path.is_file())


def plan(root):
    source = root / SOURCE
    target = root / DESTINATION
    root_text = root / ROOT_TEXT
    if not source.is_dir() or source.is_symlink():
        raise RuntimeError("Expected original Todo directory is unavailable: " + str(source))
    if target.exists() or target.is_symlink():
        raise RuntimeError("Destination already exists; refusing to overwrite: " + str(target))
    if not root_text.is_file() or root_text.is_symlink():
        raise RuntimeError("Expected original Todo text is unavailable: " + str(root_text))
    if (root / MANIFEST).exists():
        raise RuntimeError("Migration manifest already exists; use --check")

    inputs = files_under(source) + [root_text]
    inventory = [{"before": relative(path, root), "after": relative(map_path(path, root), root),
                  "sha256Before": digest(path), "bytesBefore": path.stat().st_size}
                 for path in inputs]
    changes = []
    checks = []
    for old_file in document_paths(root):
        new_file = map_path(old_file, root)
        content = old_file.read_text(encoding="utf-8")
        revised, link_checks = rewrite_markdown(content, old_file, new_file, root)
        checks.extend(dict(item, document=relative(new_file, root)) for item in link_checks)
        if revised != content:
            changes.append({"before": relative(old_file, root), "after": relative(new_file, root),
                            "sha256Before": digest(old_file), "content": revised,
                            "sha256After": hashlib.sha256(revised.encode("utf-8")).hexdigest()})
    return inventory, changes, checks


def check_manifest(root, manifest):
    failures = []
    if (root / SOURCE).exists():
        failures.append("Original directory still exists")
    if (root / ROOT_TEXT).exists():
        failures.append("Original top-level Todo text still exists")
    for item in manifest["movedFiles"]:
        path = root / item["after"]
        if not path.is_file() or digest(path) != item["sha256After"]:
            failures.append("Moved file mismatch: " + item["after"])
    for item in manifest["markdownRewrites"]:
        path = root / item["after"]
        if not path.is_file() or digest(path) != item["sha256After"]:
            failures.append("Rewritten document mismatch: " + item["after"])
    for item in manifest["navigationChecks"]:
        if item["oldTargetExisted"] and not pathlib.Path(item["newTarget"]).exists():
            failures.append("Relocated navigation target missing: " + item["newTarget"])
    if failures:
        raise RuntimeError("Migration verification failed:\n" + "\n".join(failures))
    return {"status": "PASS", "movedFiles": len(manifest["movedFiles"]),
            "markdownRewrites": len(manifest["markdownRewrites"]),
            "navigationChecks": len(manifest["navigationChecks"])}


def apply(root, inventory, changes, checks):
    # Recheck every content hash immediately before the non-overwriting relocation.
    for item in inventory + changes:
        if digest(root / item["before"]) != item["sha256Before"]:
            raise RuntimeError("Concurrent input change detected: " + item["before"])
    if (root / DESTINATION).exists() or (root / MANIFEST).exists():
        raise RuntimeError("Target appeared after planning; refusing to overwrite")
    (root / DESTINATION).parent.mkdir(parents=True, exist_ok=True)
    shutil.move(str(root / SOURCE), str(root / DESTINATION))
    shutil.move(str(root / ROOT_TEXT), str(root / DESTINATION / ROOT_TEXT.name))
    # Bulk mechanical Markdown migration, not an analysis/source/configuration edit.
    for item in changes:
        (root / item["after"]).write_text(item["content"], encoding="utf-8")
    rewritten_hashes = {item["after"]: item["sha256After"] for item in changes}
    for item in inventory:
        item["sha256After"] = rewritten_hashes.get(item["after"], item["sha256Before"])
        item["bytePreserved"] = item["sha256After"] == item["sha256Before"]
    manifest = {"schemaVersion": 1, "dateJst": "2026-10-20",
                "createdUtc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
                "source": str(SOURCE), "destination": str(DESTINATION),
                "rootTodoTextDestination": str(DESTINATION / ROOT_TEXT.name),
                "policy": "Only live Markdown navigation rewritten; non-Markdown, validation snapshots and logs byte-preserved; rootfile/job provenance untouched",
                "movedFiles": inventory,
                "markdownRewrites": [{key: value for key, value in item.items() if key != "content"} for item in changes],
                "navigationChecks": checks}
    manifest["verification"] = check_manifest(root, manifest)
    with (root / MANIFEST).open("x", encoding="utf-8") as stream:
        json.dump(manifest, stream, ensure_ascii=False, indent=2)
        stream.write("\n")
    return manifest["verification"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo-root", type=pathlib.Path, default=pathlib.Path.cwd())
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--apply", action="store_true")
    mode.add_argument("--check", action="store_true")
    args = parser.parse_args()
    root = args.repo_root.resolve()
    if not (root / "AGENTS.md").is_file():
        raise RuntimeError("--repo-root must be the star-analyzer repository")
    if args.check:
        with (root / MANIFEST).open(encoding="utf-8") as stream:
            result = check_manifest(root, json.load(stream))
    else:
        inventory, changes, checks = plan(root)
        if args.apply:
            result = apply(root, inventory, changes, checks)
        else:
            result = {"status": "PLAN_ONLY", "source": str(SOURCE), "destination": str(DESTINATION),
                      "movedFiles": len(inventory), "unchangedNonMarkdownAndSnapshots": sum(
                          1 for item in inventory if item["after"] not in {change["after"] for change in changes}),
                      "markdownRewrites": [item["after"] for item in changes],
                      "navigationChecks": len(checks)}
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
