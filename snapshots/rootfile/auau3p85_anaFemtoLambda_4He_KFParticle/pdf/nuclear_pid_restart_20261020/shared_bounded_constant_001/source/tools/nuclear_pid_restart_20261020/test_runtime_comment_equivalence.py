from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import render_tag_study as view


class RuntimeCommentEquivalenceTests(unittest.TestCase):
    def fixture(self, root):
        main = root/'config/mainconf/example.yaml'
        saved = root/'study/source/config/mainconf/example.yaml'
        for path in (main, saved):
            path.parent.mkdir(parents=True)
            path.write_text('# Original comment\nnuclearid: cuts/nuclear.yaml\n')
        expected = view.digest(saved)
        plan = dict(sourceSha256={str(main): expected}, inputs={str(saved): expected})
        return main, saved, plan

    def test_comments_are_not_physics_changes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            main, _, plan = self.fixture(root)
            main.write_text('# Git PID restored; no value changes\nnuclearid: cuts/nuclear.yaml\n')
            with patch.object(view, 'PROJECT', root):
                view.verify_runtime_yaml(plan, main)

    def test_actual_reference_change_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            main, _, plan = self.fixture(root)
            main.write_text('nuclearid: cuts/different.yaml\n')
            with patch.object(view, 'PROJECT', root), self.assertRaises(ValueError):
                view.verify_runtime_yaml(plan, main)

    def test_changed_snapshot_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            main, saved, plan = self.fixture(root)
            main.write_text('# Changed comment\nnuclearid: cuts/nuclear.yaml\n')
            saved.write_text('# Invalidated snapshot\nnuclearid: cuts/nuclear.yaml\n')
            with patch.object(view, 'PROJECT', root), self.assertRaises(ValueError):
                view.verify_runtime_yaml(plan, main)

    def test_byte_identical_reference(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            main, _, plan = self.fixture(root)
            with patch.object(view, 'PROJECT', root):
                view.verify_runtime_yaml(plan, main)


if __name__ == '__main__':
    unittest.main()
