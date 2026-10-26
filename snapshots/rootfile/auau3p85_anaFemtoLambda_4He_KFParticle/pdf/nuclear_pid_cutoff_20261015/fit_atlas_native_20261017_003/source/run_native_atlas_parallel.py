#!/usr/bin/env python3
"""目視確認済みの描画計画を最大2workerで実行し、全冊後に読戻し監査する。"""
import argparse
import concurrent.futures
from pathlib import Path
import shutil
import sys

import render_native_atlas as atlas


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--workers', required=True, type=int, choices=(1, 2))
    args = parser.parse_args()
    output = atlas.resolve(args.output_dir)
    plan = atlas.read(output / 'plan.json')
    atlas.require(plan['presentationOnly'] and not plan['refitted'], 'Only prepared display plans are supported')
    driver = Path(atlas.__file__).resolve()
    atlas.require(atlas.sha(driver) == atlas.sha(output / 'source' / driver.name),
                  'Prepared driver snapshot differs from the current rendering driver')
    source = Path(__file__).resolve()
    snapshot = output / 'source' / source.name
    atlas.require(not snapshot.exists(), 'Refusing to overwrite a prior parallel-launch source')
    shutil.copy2(str(source), str(snapshot))
    atlas.save(output / 'parallel_launch.json', dict(argv=sys.argv, cwd=str(Path.cwd()),
         workers=args.workers, sourceSha256=atlas.sha(source), createdJST=atlas.stamp(),
         presentationOnly=True, refitted=False))
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        list(executor.map(lambda book: atlas.run(output, book['id']), plan['books']))
    atlas.check(output)


if __name__ == '__main__':
    main()
