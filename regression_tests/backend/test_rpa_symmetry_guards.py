import pathlib
import tempfile
import unittest
import xml.etree.ElementTree as ET

from backend.validate import Validate


ROOT = pathlib.Path(__file__).resolve().parents[1]
CASE = "rpa_abacus_BN_headwing_sym_shrink_v1_libri"


class TestRpaSymmetryGuards(unittest.TestCase):
    def _validator(self, name):
        case = ET.parse(ROOT / "testsuite.xml").find(
            ".//testcase[@directory='{}']".format(CASE))
        entry = case.find("validate[@name='{}']".format(name))
        self.assertIsNotNone(entry, "Missing per-q RPA guard: " + name)
        return Validate(name, "librpa.out", entry.get("comparison"), None, None,
                        entry.get("regex"), None, None)

    def _compare(self, name, mutate=lambda text: text):
        reference_dir = ROOT / "refs" / CASE / "librpa"
        reference = (reference_dir / "librpa.out").read_text()
        with tempfile.TemporaryDirectory() as tmp:
            (pathlib.Path(tmp) / "librpa.out").write_text(mutate(reference))
            return self._validator(name).evaluate(tmp, str(reference_dir))[0]

    def test_unmodified_reference_passes(self):
        for name in ("Per-q EcRPA real parts", "Per-q EcRPA imaginary parts"):
            self.assertTrue(self._compare(name))

    def test_nonzero_imaginary_part_fails_with_unchanged_total(self):
        self.assertFalse(self._compare("Per-q EcRPA imaginary parts",
                                      lambda text: text.replace("7.34051e-35", "0.31", 1)))

    def test_cancelling_q_errors_fail_with_unchanged_total(self):
        self.assertFalse(self._compare("Per-q EcRPA real parts",
                                      lambda text: text.replace("-0.0886164", "-0.0786164", 1)
                                      .replace("-0.163253", "-0.173253", 1)))

    def test_nonfinite_imaginary_part_fails(self):
        self.assertFalse(self._compare("Per-q EcRPA imaginary parts",
                                      lambda text: text.replace("7.34051e-35", "nan", 1)))

    def test_missing_q_contribution_fails(self):
        self.assertFalse(self._compare("Per-q EcRPA real parts",
                                      lambda text: "\n".join(line for line in text.splitlines()
                                                              if "7.34051e-35" not in line)))


if __name__ == "__main__":
    unittest.main()
