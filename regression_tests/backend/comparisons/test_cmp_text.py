import unittest

import cmp_text


class TestCmpText(unittest.TestCase):
    def test_equal_runtime_flags_pass(self):
        values = {"flags": ["exx true", "gw true", "rpa true"]}
        self.assertTrue(cmp_text.equal(values, values)[0])

    def test_disabled_runtime_flag_fails(self):
        self.assertFalse(cmp_text.equal({"flags": ["gw false"]},
                                         {"flags": ["gw true"]})[0])

    def test_missing_or_empty_flags_fail(self):
        for values in ({}, {"flags": []}, {"flags": [""]}):
            with self.subTest(values=values):
                self.assertFalse(cmp_text.equal(values, values)[0])

    def test_missing_one_flag_fails(self):
        self.assertFalse(cmp_text.equal({"flags": ["exx true"]},
                                         {"flags": ["exx true", "gw true"]})[0])


if __name__ == "__main__":
    unittest.main()
