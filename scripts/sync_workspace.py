#!/usr/bin/env python3
"""Version personal notes directly and snapshot private workspace files safely."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

REPO = Path(__file__).resolve().parents[1]
DENIED = {'.root', '.o', '.so', '.dylib', '.pcm', '.d', '.log', '.out',
          '.err', '.gz', '.tar', '.zip', '.pyc', '.pyo', '.aux', '.toc'}
IMAGES = {'.pdf', '.png', '.jpg', '.jpeg'}
PRUNE = {'.git', '__pycache__', 'build', 'lib', '.aws', '.ssh', '.codex',
         'subjobs', 'input_cache', 'node_modules'}


def git(where, *args):
    return subprocess.check_output(['git', '-C', str(where)] + list(args))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def policy():
    value = json.loads((REPO / 'policy.json').read_text())
    if value['schemaVersion'] != 1:
        raise ValueError('Unsupported policy version')
    for name in value['snapshotScopes'] + value['extraPaths'] + list(value['extraLargeReports']):
        safe_relative(name)
    return value


def safe_relative(name):
    p = Path(name)
    if p.is_absolute() or '..' in p.parts or not p.parts or '.git' in p.parts:
        raise ValueError('Unsafe relative path: ' + name)
    return p


def permitted(path, relative, cfg, output=False):
    if path.is_symlink() or not path.is_file():
        return False, 'symlink_or_nonfile'
    ext = path.suffix.lower()
    if ext in DENIED or path.name.startswith('AutoDict') or '_C_ACLiC_dict' in path.name:
        return False, 'generated'
    if ext not in cfg['sourceExtensions'] and relative not in cfg['extraPaths']:
        return False, 'extension'
    maximum = cfg['extraLargeReports'].get(relative, cfg['maxFileBytes'])
    if output and ext in {'.json', '.tsv', '.csv'}:
        maximum = min(maximum, cfg['maxOutputMetadataBytes'])
    if path.stat().st_size > maximum:
        return False, 'size'
    if output and ext in {'.tsv', '.csv', '.list'}:
        return False, 'numeric_output'
    if output and ext == '.json' and path.name not in cfg['outputMetadataNames']:
        return False, 'numeric_output'
    if ext not in IMAGES and b'\0' in path.read_bytes():
        return False, 'binary'
    return True, ''


def files_under(path):
    if path.is_symlink():
        return
    if path.is_file():
        yield path
        return
    if not path.is_dir():
        return
    for base, dirs, names in os.walk(str(path), followlinks=False):
        dirs[:] = sorted(d for d in dirs if d not in PRUNE and not (Path(base) / d).is_symlink())
        for name in sorted(names):
            yield Path(base) / name


def collect(workspace, cfg):
    tracked = {x.decode() for x in git(workspace, 'ls-files', '-z').split(b'\0') if x}
    modified = {x.decode() for x in git(workspace, 'diff', '--name-only', '-z', 'HEAD').split(b'\0') if x}
    candidates = {}
    for scope in cfg['snapshotScopes']:
        for p in files_under(workspace / scope):
            relative = str(p.relative_to(workspace))
            if relative not in tracked or relative in modified:
                candidates[relative] = ('local_override' if relative in tracked else 'private_source', False)
    for relative in cfg['extraPaths']:
        candidates[relative] = ('personal_document', False)
    for p in files_under(workspace / 'job/run'):
        if p.suffix.lower() in {'.py', '.sh', '.c', '.h'}:
            relative = str(p.relative_to(workspace))
            if relative not in tracked or relative in modified:
                candidates[relative] = ('private_job_helper', False)
    outputs = workspace / 'rootfile'
    if outputs.is_dir():
        for ana in sorted(outputs.iterdir()):
            if not ana.is_dir() or ana.is_symlink():
                continue
            for p in ana.iterdir():
                if p.is_file() and p.suffix.lower() == '.md':
                    candidates[str(p.relative_to(workspace))] = ('output_document', True)
            for folder in ['pdf', 'mdfiles']:
                for p in files_under(ana / folder):
                    candidates[str(p.relative_to(workspace))] = ('small_output', True)
    def inspect(item):
        relative, (kind, output) = item
        p = workspace / safe_relative(relative)
        ok, reason = permitted(p, relative, cfg, output)
        if ok:
            return relative, {'sha256': sha(p), 'bytes': p.stat().st_size, 'kind': kind}, None
        elif reason in {'size', 'binary'}:
            return relative, None, {'path': relative, 'reason': reason, 'bytes': p.stat().st_size}
        return relative, None, None
    accepted, skipped = {}, []
    # Files are independent; overlap metadata/read latency on the shared filesystem.
    with ThreadPoolExecutor(max_workers=8) as pool:
        for i, (relative, record, rejected) in enumerate(pool.map(inspect, sorted(candidates.items())), 1):
            if record:
                accepted[relative] = record
            if rejected:
                skipped.append(rejected)
            if i % 1000 == 0:
                print('Inspected {}/{} candidate files'.format(i, len(candidates)), file=sys.stderr, flush=True)
    return accepted, skipped


def capture(workspace, cfg, write=False):
    records, skipped = collect(workspace, cfg)
    summary = {'files': len(records), 'bytes': sum(v['bytes'] for v in records.values()),
               'oversizeOrBinaryFiles': len(skipped)}
    if write:
        previous_path = REPO / 'metadata/snapshots.json'
        previous = json.loads(previous_path.read_text())['files'] if previous_path.exists() else {}
        for relative, record in records.items():
            source, dest = workspace / relative, REPO / 'snapshots' / relative
            if dest.is_symlink():
                raise ValueError('Refuse snapshot symlink: ' + str(dest))
            if dest.exists() and relative in previous and sha(dest) not in {previous[relative]['sha256'], record['sha256']}:
                raise ValueError('Snapshot also edited; resolve before capture: ' + relative)
        for relative in previous.keys() - records.keys():
            dest = REPO / 'snapshots' / safe_relative(relative)
            if dest.exists() and sha(dest) != previous[relative]['sha256']:
                raise ValueError('Refuse to remove edited snapshot: ' + relative)
        for relative, record in records.items():
            source, dest = workspace / relative, REPO / 'snapshots' / relative
            if not dest.exists() or sha(dest) != record['sha256']:
                dest.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(source, dest)
                if sha(dest) != record['sha256']:
                    raise ValueError('Source changed during capture; retry before committing: ' + relative)
        for relative in previous.keys() - records.keys():
            dest = REPO / 'snapshots' / safe_relative(relative)
            if dest.exists():
                dest.unlink()
        previous_path.parent.mkdir(parents=True, exist_ok=True)
        previous_path.write_text(json.dumps({'schemaVersion': 1,
            'parentBaseCommit': git(workspace, 'rev-parse', 'HEAD').decode().strip(),
            'files': records}, indent=2, ensure_ascii=False) + '\n')
        (REPO / 'metadata/excluded_large_files.json').write_text(json.dumps(skipped, indent=2, ensure_ascii=False) + '\n')
    print(json.dumps(summary, sort_keys=True))
    return records


def link_missing(path, target):
    if path.is_symlink() and path.resolve() == target.resolve():
        return
    if path.exists() or path.is_symlink():
        raise ValueError('Existing path is preserved; resolve manually: ' + str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.symlink_to(os.path.relpath(str(target), str(path.parent)), target_is_directory=True)


def mount(workspace, cfg):
    for name in cfg['directDirectories']:
        (REPO / name).mkdir(exist_ok=True)
        link_missing(workspace / name, REPO / name)
    for name in cfg['bridgeDirectories']:
        if (workspace / name).exists():
            link_missing(REPO / name, workspace / name)


def restore(workspace, cfg):
    records = json.loads((REPO / 'metadata/snapshots.json').read_text())['files']
    pending = []
    for relative, record in records.items():
        source, dest = REPO / 'snapshots' / safe_relative(relative), workspace / safe_relative(relative)
        if sha(source) != record['sha256']:
            raise ValueError('Snapshot checksum mismatch: ' + relative)
        # Resolve parents as well: never restore through an escaping directory link.
        if os.path.commonpath([str(dest.parent.resolve()), str(workspace.resolve())]) != str(workspace.resolve()):
            raise ValueError('Destination escapes workspace: ' + relative)
        if dest.exists() or dest.is_symlink():
            if dest.is_symlink() or not dest.is_file() or sha(dest) != record['sha256']:
                raise ValueError('Existing different file preserved: ' + relative)
        else:
            pending.append((source, dest))
    for source, dest in pending:
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, dest)
    print('Restored {} missing files; existing identical files kept'.format(len(pending)))


def check_index(cfg):
    staged = [x.decode() for x in git(REPO, 'diff', '--cached', '--name-only', '--diff-filter=ACMR', '-z').split(b'\0') if x]
    entries = {}
    for row in git(REPO, 'ls-files', '--stage', '-z').split(b'\0'):
        if row:
            info, name = row.split(b'\t', 1)
            mode, oid, stage = info.decode().split()
            if stage != '0':
                raise ValueError('Unmerged index: ' + name.decode())
            entries[name.decode()] = (mode, oid)
    oids = sorted({entries[name][1] for name in staged})
    query = ('\n'.join(oids) + '\n').encode() if oids else b''
    sizes = {}
    if query:
        result = subprocess.run(['git', '-C', str(REPO), 'cat-file', '--batch-check'],
                                input=query, stdout=subprocess.PIPE, check=True)
        for row in result.stdout.decode().splitlines():
            oid, kind, size = row.split()
            if kind != 'blob':
                raise ValueError('Only regular blobs may be staged')
            sizes[oid] = int(size)
    text_oids = set()
    for name in staged:
        safe_relative(name)
        mode, oid = entries[name]
        if mode not in {'100644', '100755'}:
            raise ValueError('Only regular files may be committed: ' + name)
        ext = Path(name).suffix.lower()
        maximum = cfg['maxManifestBytes'] if name.startswith('metadata/') else cfg['maxFileBytes']
        original = name[len('snapshots/'):] if name.startswith('snapshots/') else name
        maximum = cfg['extraLargeReports'].get(original, maximum)
        if ext in DENIED or sizes[oid] > maximum:
            raise ValueError('Generated or oversized file in index: ' + name)
        if ext not in IMAGES:
            text_oids.add(oid)
    if text_oids:
        result = subprocess.run(['git', '-C', str(REPO), 'cat-file', '--batch'],
            input=('\n'.join(sorted(text_oids)) + '\n').encode(), stdout=subprocess.PIPE, check=True)
        offset = 0
        for _ in text_oids:
            end = result.stdout.index(b'\n', offset)
            oid, kind, size = result.stdout[offset:end].decode().split()
            start, size = end + 1, int(size)
            if b'\0' in result.stdout[start:start + size]:
                raise ValueError('Binary content in a text file: ' + oid)
            offset = start + size + 1
    print('PASS staged-file policy: {} files'.format(len(staged)))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['plan', 'capture', 'mount', 'restore', 'check-index'])
    parser.add_argument('--workspace', type=Path, default=REPO.parent)
    args = parser.parse_args()
    cfg = policy()
    workspace = args.workspace.resolve()
    if args.action == 'check-index':
        check_index(cfg)
    elif args.action in {'plan', 'capture'}:
        capture(workspace, cfg, write=args.action == 'capture')
    elif args.action == 'mount':
        mount(workspace, cfg)
    else:
        restore(workspace, cfg)


if __name__ == '__main__':
    try:
        main()
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        print('ERROR: ' + str(error), file=sys.stderr)
        sys.exit(1)
