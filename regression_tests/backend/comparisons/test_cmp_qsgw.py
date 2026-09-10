import unittest

import cmp_qsgw


CONTRACT_HEADER = (
    "# qsgw_contract_version 6\n"
    "# fixed_basis immutable_mf0\n"
    "# live_update eigenvalues_wfc\n"
    "# velocity disabled_stage1\n"
    "# head disabled_stage1\n"
    "# wing disabled_stage1\n"
    "# symmetry exx_off_gw_off_rpa_off\n"
    "# hartree disabled\n"
    "# band disabled_stage1\n"
    "# h_qsgw_cut disabled_non_band\n"
    "# qsgw_input_contract qsgw_input.contract\n"
    "# qsgw_mixer linear\n"
    "# qsgw_mixing_beta 0.2\n"
)
EIGENVALUE_HEADER = CONTRACT_HEADER + (
    "# iter channel spin kpoint kx ky kz band energy_eV\n"
)
SUMMARY_HEADER = CONTRACT_HEADER + (
    "# iter max_delta_eV residual_l2_Ha residual_max_Ha "
    "efermi_eV gap_eV electron_count requested_mode applied_mode beta "
    "fallback rcond coefficient_l1 coefficient_count converged "
    "coefficients fallback_reason\n"
)


class TestQsgwEigenvalueTrace(unittest.TestCase):

    def _trace(self, energy_ev, kx=0.0):
        return EIGENVALUE_HEADER + (
            "0 0 0 0 {:.17e} 0.0 0.0 0 5.00000000000000000e-1\n"
            "1 0 0 0 {:.17e} 0.0 0.0 0 {:.17e}\n"
        ).format(kx, kx, energy_ev)

    def _compare(self, test, reference, **kwargs):
        compare = cmp_qsgw.eigenvalue_trace(**kwargs)
        return compare(
            {"qsgw_eigenvalues.dat": test},
            {"qsgw_eigenvalues.dat": reference},
        )

    def test_energy_difference_within_hartree_tolerance_passes(self):
        reference = self._trace(1.0)
        test = self._trace(
            1.0 + 0.5 * cmp_qsgw.HA2EV * 1.0e-6)

        passed, message = self._compare(
            test, reference, tolerance_ha="1e-6")

        self.assertTrue(passed, message)
        self.assertIn("max abs eigenvalue diff", message)

    def test_energy_difference_above_hartree_tolerance_fails(self):
        reference = self._trace(1.0)
        test = self._trace(
            1.0 + 2.0 * cmp_qsgw.HA2EV * 1.0e-6)

        passed, message = self._compare(
            test, reference, tolerance_ha="1e-6")

        self.assertFalse(passed)
        self.assertIn("max abs eigenvalue diff", message)

class TestQsgwIterationSummary(unittest.TestCase):

    def _trace(self, gap=1.0, applied_mode=0):
        return SUMMARY_HEADER + (
            "0 0.0 0.0 0.0 -1.0 {:.17e} 8.0 "
            "-1 -1 2.0e-1 0 1.0 0.0 0 0 none none\n"
            "1 1.0e-3 2.0e-4 1.0e-4 -1.0 {:.17e} 8.0 "
            "0 {} 2.0e-1 0 1.0 1.0 1 0 1.0 none\n"
        ).format(gap, gap, applied_mode)

    def _compare(self, test, reference, **kwargs):
        compare = cmp_qsgw.iteration_summary(**kwargs)
        return compare(
            {"qsgw_iterations.dat": test},
            {"qsgw_iterations.dat": reference},
        )

    def test_summary_within_field_tolerances_passes(self):
        passed, message = self._compare(
            self._trace(gap=1.0 + 5.0e-6),
            self._trace(gap=1.0),
            energy_tolerance_ev="1e-5",
        )

        self.assertTrue(passed, message)
        self.assertIn("max gap diff", message)

    def test_gap_above_tolerance_fails(self):
        passed, message = self._compare(
            self._trace(gap=1.0 + 2.0e-5),
            self._trace(gap=1.0),
            energy_tolerance_ev="1e-5",
        )

        self.assertFalse(passed)
        self.assertIn("gap", message)

if __name__ == "__main__":
    unittest.main()
