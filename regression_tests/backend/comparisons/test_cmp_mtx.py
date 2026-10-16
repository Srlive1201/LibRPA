import unittest

import cmp_mtx


REAL_MTX = """%%MatrixMarket matrix coordinate real general
%
2 2 2
2 1 3.0
1 2 1.0
"""

ZERO_MTX = """%%MatrixMarket matrix coordinate real general
2 3 0
"""


class TestCmpMtx(unittest.TestCase):

    def _compare(self, matrix1, matrix2, **kwargs):
        compare = cmp_mtx.abs_diff(1e-4, **kwargs)
        return compare({"matrix.mtx": matrix1}, {"matrix.mtx": matrix2})

    def test_no_extracted_files_fails(self):
        compare = cmp_mtx.abs_diff(1e-4)
        passed, msg = compare({}, {})
        self.assertFalse(passed)
        self.assertIn("no files found", msg)

    def test_no_extracted_matrices_fails(self):
        compare = cmp_mtx.abs_diff(1e-4)
        passed, msg = compare({"matrix.mtx": []}, {"matrix.mtx": []})
        self.assertFalse(passed)
        self.assertIn("no matrices found", msg)

    def test_zero_entry_matrix_passes(self):
        passed, _ = self._compare(ZERO_MTX, ZERO_MTX)
        self.assertTrue(passed)

    def test_coordinate_order_does_not_matter(self):
        reordered = """%%MatrixMarket matrix coordinate real general
2 2 2
1 2 1.00001
2 1 3.00001
"""
        passed, _ = self._compare(REAL_MTX, reordered)
        self.assertTrue(passed)

    def test_complex_values_use_magnitude_difference(self):
        matrix1 = """%%MatrixMarket matrix coordinate complex general
1 1 1
1 1 1.0 2.0
"""
        matrix2 = """%%MatrixMarket matrix coordinate complex general
1 1 1
1 1 1.00003 2.00004
"""
        passed, _ = self._compare(matrix1, matrix2)
        self.assertTrue(passed)

    def test_header_mismatch_fails(self):
        other_shape = REAL_MTX.replace("2 2 2", "3 2 2")
        passed, msg = self._compare(REAL_MTX, other_shape)
        self.assertFalse(passed)
        self.assertIn("matrix header mismatch", msg)

    def test_comments_can_be_compared(self):
        other_comment = REAL_MTX.replace("%\n", "% generated elsewhere\n")
        passed, msg = self._compare(REAL_MTX, other_comment, comments="true")
        self.assertFalse(passed)
        self.assertIn("matrix comments mismatch", msg)

    def test_nan_in_one_matrix_fails(self):
        matrix1 = REAL_MTX.replace("3.0", "nan")
        passed, msg = self._compare(matrix1, REAL_MTX)
        self.assertFalse(passed)
        self.assertIn("nan mismatch", msg)

    def test_matching_nonfinite_entries_fail(self):
        for value in ("nan", "inf", "-inf"):
            with self.subTest(value=value):
                matrix = REAL_MTX.replace("3.0", value)
                passed, _ = self._compare(matrix, matrix)
                self.assertFalse(passed)

    def test_nonfinite_tolerance_is_rejected(self):
        for tolerance in ("nan", "inf", "-1"):
            with self.subTest(tolerance=tolerance):
                with self.assertRaises(ValueError):
                    cmp_mtx.abs_diff(tolerance)


if __name__ == "__main__":
    unittest.main()
