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
    "# qsgw_input_contract_sha256 " + "a" * 64 + "\n"
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

    def test_kpoint_coordinate_mismatch_fails(self):
        passed, message = self._compare(
            self._trace(1.0, kx=1.0e-6),
            self._trace(1.0),
            coordinate_tolerance="1e-12",
        )

        self.assertFalse(passed)
        self.assertIn("k-point coordinate", message)

    def test_band_contract_requires_band_trace_channel(self):
        header = EIGENVALUE_HEADER.replace(
            "# band disabled_stage1\n"
            "# h_qsgw_cut disabled_non_band",
            "# band fixed_reference_rotation_live\n"
            "# h_qsgw_cut band_postprocess\n"
            "# qsgw_band0_unoccupied_keep 2\n"
            "# qsgw_band0_cut_mode 0\n"
            "# qsgw_band0_cut_shift_ha 20",
        )
        trace = header + (
            "0 0 0 0 0.0 0.0 0.0 0 0.5\n"
            "1 0 0 0 0.0 0.0 0.0 0 1.0\n"
        )

        passed, message = self._compare(trace, trace)

        self.assertFalse(passed)
        self.assertIn("channel set", message)


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

    def test_discrete_mixing_metadata_mismatch_fails(self):
        passed, message = self._compare(
            self._trace(applied_mode=1),
            self._trace(applied_mode=0),
        )

        self.assertFalse(passed)
        self.assertIn("applied_mode", message)

    def test_contract_mismatch_fails(self):
        test = self._trace().replace(
            "# qsgw_mixing_beta 0.2",
            "# qsgw_mixing_beta 0.3",
        )

        passed, message = self._compare(test, self._trace())

        self.assertFalse(passed)
        self.assertIn("contract differs", message)

    def test_missing_contract_fails(self):
        trace = self._trace().replace(CONTRACT_HEADER, "")

        passed, message = self._compare(trace, trace)

        self.assertFalse(passed)
        self.assertIn("missing QSGW contract", message)

    def test_iteration_zero_is_required(self):
        trace = self._trace().replace(
            "0 0.0 0.0 0.0 -1.0 1.00000000000000000e+00 8.0 "
            "-1 -1 2.0e-1 0 1.0 0.0 0 0 none none\n",
            "",
        )

        passed, message = self._compare(trace, trace)

        self.assertFalse(passed)
        self.assertIn("iteration zero", message)

    def test_supported_head_only_contract_passes(self):
        trace = (
            self._trace()
            .replace(
                "# velocity disabled_stage1",
                "# velocity fixed_reference",
            )
            .replace(
                "# head disabled_stage1",
                "# head scf_grid_analytic_live",
            )
        )

        passed, message = self._compare(trace, trace)

        self.assertTrue(passed, message)

    def test_invalid_symmetry_or_head_velocity_contract_fails(self):
        invalid_symmetry = self._trace().replace(
            "# symmetry exx_off_gw_off_rpa_off",
            "# symmetry crystal_reduction",
        )
        passed, message = self._compare(
            invalid_symmetry, invalid_symmetry)
        self.assertFalse(passed)
        self.assertIn("symmetry", message)

        inconsistent_head = self._trace().replace(
            "# head disabled_stage1",
            "# head scf_grid_analytic_live",
        )
        passed, message = self._compare(
            inconsistent_head, inconsistent_head)
        self.assertFalse(passed)
        self.assertIn("head/velocity", message)

    def test_unsupported_hartree_wing_or_legacy_contract_fails(self):
        hartree = self._trace().replace(
            "# hartree disabled", "# hartree delta_density")
        passed, message = self._compare(hartree, hartree)
        self.assertFalse(passed)
        self.assertIn("hartree", message)

        wing = self._trace().replace(
            "# wing disabled_stage1",
            "# wing scf_grid_analytic_live",
        )
        passed, message = self._compare(wing, wing)
        self.assertFalse(passed)
        self.assertIn("wing", message)

        legacy = self._trace().replace(
            "# qsgw_contract_version 6",
            "# qsgw_contract_version 5",
        )
        passed, message = self._compare(legacy, legacy)
        self.assertFalse(passed)
        self.assertIn("version", message)


if __name__ == "__main__":
    unittest.main()
