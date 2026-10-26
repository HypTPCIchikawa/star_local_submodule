import importlib.util
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest import mock

MODULE = Path(__file__).resolve().parents[1] / 'scripts/sync_workspace.py'
spec = importlib.util.spec_from_file_location('personal_sync', MODULE)
sync = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sync)


class PersonalSync(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.base = Path(self.temp.name)
        self.repo, self.workspace = self.base / 'private', self.base / 'shared'
        self.repo.mkdir()
        self.workspace.mkdir()
        self.cfg = sync.policy()

    def test_root_and_generated_files_rejected_even_when_small(self):
        for name in ['data.root', 'build.so', 'job.log', 'AutoDict_x.cxx']:
            p = self.base / name
            p.write_bytes(b'data')
            self.assertFalse(sync.permitted(p, name, self.cfg)[0])

    def test_small_pdf_and_png_allowed_and_large_pdf_rejected(self):
        p = self.base / 'report.pdf'
        p.write_bytes(b'%PDF-1.4\0image')
        self.assertTrue(sync.permitted(p, p.name, self.cfg, True)[0])
        p.write_bytes(b'x' * (self.cfg['maxFileBytes'] + 1))
        self.assertEqual(sync.permitted(p, p.name, self.cfg, True), (False, 'size'))

    def test_binary_and_symlink_rejected(self):
        p = self.base / 'bad.py'
        p.write_bytes(b'a\0b')
        self.assertEqual(sync.permitted(p, p.name, self.cfg), (False, 'binary'))
        link = self.base / 'link.py'
        link.symlink_to(p)
        self.assertFalse(sync.permitted(link, link.name, self.cfg)[0])

    def test_output_numeric_arrays_excluded_but_summary_saved(self):
        p = self.base / 'summary.json'
        p.write_text('{}')
        self.assertTrue(sync.permitted(p, p.name, self.cfg, True)[0])
        for name in ['spectra.json', 'native_counts.tsv']:
            p = self.base / name
            p.write_text('1')
            self.assertEqual(sync.permitted(p, p.name, self.cfg, True)[1], 'numeric_output')

    def test_mount_refuses_existing_directory(self):
        original = self.workspace / 'Todo'
        original.mkdir()
        (original / 'keep.md').write_text('keep')
        with self.assertRaises(ValueError):
            sync.link_missing(original, self.repo / 'Todo')
        self.assertEqual((original / 'keep.md').read_text(), 'keep')

    def test_restore_preflights_conflict_before_writing_any_file(self):
        records = {}
        for name in ['a.py', 'b.py']:
            p = self.repo / 'snapshots' / name
            p.parent.mkdir(exist_ok=True)
            p.write_text('saved')
            records[name] = {'sha256': sync.sha(p)}
        (self.repo / 'metadata').mkdir()
        (self.repo / 'metadata/snapshots.json').write_text(json.dumps({'files': records}))
        (self.workspace / 'b.py').write_text('new work')
        with mock.patch.object(sync, 'REPO', self.repo), self.assertRaises(ValueError):
            sync.restore(self.workspace, self.cfg)
        self.assertFalse((self.workspace / 'a.py').exists())
        self.assertEqual((self.workspace / 'b.py').read_text(), 'new work')

    def test_restore_missing_root_file_and_reject_path_escape(self):
        p = self.repo / 'snapshots/AGENTS.md'
        p.parent.mkdir()
        p.write_text('saved')
        (self.repo / 'metadata').mkdir()
        (self.repo / 'metadata/snapshots.json').write_text(json.dumps({'files': {'AGENTS.md': {'sha256': sync.sha(p)}}}))
        with mock.patch.object(sync, 'REPO', self.repo):
            sync.restore(self.workspace, self.cfg)
        self.assertEqual((self.workspace / 'AGENTS.md').read_text(), 'saved')
        for name in ['../outside.txt', '/absolute.txt', '.git/config']:
            with self.assertRaises(ValueError):
                sync.safe_relative(name)

    def test_recursive_update_none_is_optional(self):
        parent, checkout = self.base / 'parent', self.base / 'checkout'
        subprocess.run(['git', 'init', '-q', str(parent)], check=True)
        subprocess.run(['git', '-C', str(parent), 'config', 'user.name', 'test'], check=True)
        subprocess.run(['git', '-C', str(parent), 'config', 'user.email', 'test@example.invalid'], check=True)
        (parent / '.gitmodules').write_text('[submodule "personal"]\n path = personal\n url = /missing/private-repo\n update = none\n')
        subprocess.run(['git', '-C', str(parent), 'add', '.gitmodules'], check=True)
        subprocess.run(['git', '-C', str(parent), 'commit', '-qm', 'base'], check=True)
        revision = sync.git(parent, 'rev-parse', 'HEAD').decode().strip()
        subprocess.run(['git', '-C', str(parent), 'update-index', '--add', '--cacheinfo', '160000,' + revision + ',personal'], check=True)
        subprocess.run(['git', '-C', str(parent), 'commit', '-qm', 'optional module'], check=True)
        result = subprocess.run(['git', '-c', 'protocol.file.allow=always', 'clone', '-q', '--recurse-submodules', str(parent), str(checkout)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        self.assertEqual(result.returncode, 0, result.stderr.decode())
        self.assertFalse((checkout / 'personal/.git').exists())

    def test_batch_index_guard_rejects_root_and_accepts_small_pdf(self):
        subprocess.run(['git', 'init', '-q', str(self.repo)], check=True)
        (self.repo / 'report.pdf').write_bytes(b'%PDF\0content')
        subprocess.run(['git', '-C', str(self.repo), 'add', 'report.pdf'], check=True)
        with mock.patch.object(sync, 'REPO', self.repo):
            sync.check_index(self.cfg)
        (self.repo / 'data.root').write_bytes(b'root')
        subprocess.run(['git', '-C', str(self.repo), 'add', 'data.root'], check=True)
        with mock.patch.object(sync, 'REPO', self.repo), self.assertRaises(ValueError):
            sync.check_index(self.cfg)


if __name__ == '__main__':
    unittest.main()
