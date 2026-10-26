from pathlib import Path
import re
import unittest

class ReportLinks(unittest.TestCase):
    def test_current_study_links_when_raw_data_are_present(self):
        project=Path(__file__).resolve().parents[2]
        base=project/'rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020'
        if not (base/'shared_bounded_linear_001/views_001/completed.json').exists():
            self.skipTest('Analysis products are intentionally not version-controlled')
        report=project/'Todo/2026/10/20/summary_polynomial_all_purity_20261020.md'
        targets=re.findall(r'\]\(([^)]+)\)',report.read_text())
        self.assertGreater(len(targets),10)
        for target in targets:
            if '://' not in target and not target.startswith('#'):
                self.assertTrue((report.parent/target.split('#')[0]).resolve().exists(),target)
        for out in ('shared_purity_summary_002','shared_bounded_purity_summary_001','transfer_comparison_001','transfer_comparison_bounded_001'):
            self.assertTrue((base/out/'completed.json').is_file(),out)

if __name__=='__main__':unittest.main()
