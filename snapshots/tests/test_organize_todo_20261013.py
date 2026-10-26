"""No ROOT dependency: one-off Todo navigation migration regression checks."""

import importlib.util
import pathlib
import tempfile
import unittest


SCRIPT = pathlib.Path(__file__).resolve().parents[1] / "script/organize_todo_20261013.py"
SPEC = importlib.util.spec_from_file_location("organize_todo_20261013", SCRIPT)
MIGRATION = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MIGRATION)


class NavigationTest(unittest.TestCase):
    def setUp(self):
        self.root = pathlib.Path("/repo")
        self.old = self.root / "Todo/20261013/README.md"
        self.new = self.root / "Todo/2026/10/13/README.md"

    def rewrite(self, value, old=None, new=None):
        return MIGRATION.rewrite_markdown(value, old or self.old, new or self.new, self.root)[0]

    def test_external_root_link_gains_two_levels(self):
        self.assertEqual(self.rewrite("[ROOT](../../rootfile/a.pdf)\n"),
                         "[ROOT](../../../../rootfile/a.pdf)\n")

    def test_sibling_stays_sibling(self):
        self.assertEqual(self.rewrite("[plan](task3_plan.md#scope)"), "[plan](task3_plan.md#scope)")

    def test_top_text_moves_into_day(self):
        self.assertEqual(self.rewrite("[Todo](../STAR_Todo_20261013.txt)"),
                         "[Todo](STAR_Todo_20261013.txt)")

    def test_inbound_day_link(self):
        old = self.root / "Todo/2026/10/14/README.md"
        self.assertEqual(self.rewrite("[prior](../../../20261013/README.md)", old, old),
                         "[prior](../13/README.md)")

    def test_analysis_day_is_not_the_todo_day(self):
        old = self.root / "analysisnote/20261014/summary20261014.md"
        value = "[log](../20261013/summary20261013.md)"
        self.assertEqual(self.rewrite(value, old, old), value)

    def test_absolute_user_mount_kept(self):
        value = "[Todo](/star/u/name/star-analyzer/Todo/20261013/README.md)"
        self.assertEqual(self.rewrite(value),
                         "[Todo](/star/u/name/star-analyzer/Todo/2026/10/13/README.md)")

    def test_fenced_history_kept(self):
        value = "```bash\ncat Todo/20261013/a.md\n[old](../../rootfile/a.pdf)\n```\n"
        self.assertEqual(self.rewrite(value), value)

    def test_inline_navigation_only(self):
        value = "See `Todo/20261013/a.md`; historical command `cat Todo/20261013/a.md`."
        self.assertEqual(self.rewrite(value),
                         "See `Todo/2026/10/13/a.md`; historical command `cat Todo/20261013/a.md`.")

    def test_reference_link(self):
        self.assertEqual(self.rewrite("[output]: ../../rootfile/a.pdf\n"),
                         "[output]: ../../../../rootfile/a.pdf\n")

    def test_angle_destination_preserved(self):
        self.assertEqual(self.rewrite("[output](<../../rootfile/my report.pdf>)"),
                         "[output](<../../../../rootfile/my report.pdf>)")


class EndToEndTest(unittest.TestCase):
    def test_binary_and_text_preserved_and_checked(self):
        with tempfile.TemporaryDirectory() as directory:
            root = pathlib.Path(directory)
            source = root / MIGRATION.SOURCE
            (source / "validation").mkdir(parents=True)
            (root / "analysisnote/20261013").mkdir(parents=True)
            (root / "analysisnote/20261013/summary.md").write_text("history", encoding="utf-8")
            (source / "README.md").write_text("[log](../../analysisnote/20261013/summary.md)\n[Todo](../STAR_Todo_20261013.txt)\n", encoding="utf-8")
            snapshot = b"\x00\xffsource snapshot\r\n"
            (source / "validation/source.C").write_bytes(snapshot)
            (root / MIGRATION.ROOT_TEXT).write_bytes(b"Todo original\r\n")
            inventory, changes, checks = MIGRATION.plan(root)
            result = MIGRATION.apply(root, inventory, changes, checks)
            self.assertEqual(result["status"], "PASS")
            self.assertEqual((root / MIGRATION.DESTINATION / "validation/source.C").read_bytes(), snapshot)
            self.assertEqual((root / MIGRATION.DESTINATION / MIGRATION.ROOT_TEXT.name).read_bytes(), b"Todo original\r\n")
            self.assertFalse((root / MIGRATION.SOURCE).exists())
            self.assertFalse((root / MIGRATION.ROOT_TEXT).exists())
            with self.assertRaises(RuntimeError):
                MIGRATION.plan(root)


if __name__ == "__main__":
    unittest.main()
