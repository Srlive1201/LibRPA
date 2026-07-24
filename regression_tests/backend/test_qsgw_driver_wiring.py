from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
QSGW_DRIVER = REPO_ROOT / "driver" / "tasks" / "qsgw.cpp"


def function_body(source: str, signature: str) -> str:
    start = source.index(signature)
    opening = source.index("{", start)
    depth = 0
    for position in range(opening, len(source)):
        if source[position] == "{":
            depth += 1
        elif source[position] == "}":
            depth -= 1
            if depth == 0:
                return source[opening : position + 1]
    raise AssertionError(f"unterminated function: {signature}")


def test_qsgw_reuses_upstream_gw_in_an_immutable_reference_basis() -> None:
    source = QSGW_DRIVER.read_text()
    runner = function_body(source, "void run_qsgw_stage_one(")
    normalized = " ".join(runner.split())

    assert "const MeanField reference = dataset->mf;" in runner
    assert "ScopedReferenceEigenvectors fixed_basis_projection(" in runner
    assert "h.build_g0w0_sigma(opts);" in runner
    assert "dataset->p_exx->build_KS_kgrid_blacs(" in runner
    assert "dataset->p_g0w0->build_sigc_matrix_KS_kgrid_blacs(" in runner
    assert (
        "diagonalize_in_reference_basis( dataset->mf, reference, "
        "mixed_hamiltonian" in normalized
    )


def test_qsgw_supports_only_same_grid_analytic_head_updates() -> None:
    source = QSGW_DRIVER.read_text()
    runner = function_body(source, "void run_qsgw_stage_one(")

    assert "opts.option_dielect_func == 4" in runner
    assert "QSGW independent PyATB head updates are unsupported" in runner
    assert "read_headwing_input(driver_params.input_dir, false);" in runner
    assert "if (update_head && iteration > 1)" in runner
    assert "dataset->velocity_matrix = reference_velocity;" in runner
    assert "dataset->p_headwing.reset();" in runner
    assert runner.index("dataset->velocity_matrix = reference_velocity;") < (
        runner.index("dataset->p_headwing.reset();")
    ) < runner.index("h.build_g0w0_sigma(opts);")
    assert (
        "initialize_ds_headwing(*dataset, opts, false);" not in runner
    )


def test_qsgw_band_supplies_actual_r_legacy_bvk_mapping_to_upstream() -> None:
    source = QSGW_DRIVER.read_text()
    assert '#include "../../src/qsgw/band_bvk_remap.h"' in source
    runner = function_body(source, "void run_qsgw_stage_one(")
    assert "build_legacy_band_bvk_remap(" in runner
    assert runner.count("dataset->kfrac_band_list, actual_r_bvk_remap,") == 2


def test_qsgw_band_uses_the_fixed_band_reference_and_writes_each_iteration() -> None:
    source = QSGW_DRIVER.read_text()
    runner = function_body(source, "void run_qsgw_stage_one(")
    normalized = " ".join(runner.split())

    assert "band_reference = dataset->mf_band;" in runner
    assert "dataset->p_exx->build_KS_band_blacs(" in runner
    assert "dataset->p_g0w0->build_sigc_matrix_KS_band_blacs(" in runner
    assert (
        "diagonalize_in_reference_basis( dataset->mf_band, "
        "*band_reference, mixed_band_hamiltonian)" in normalized
    )
    assert "write_qsgw_band_spin_tables(" in runner
    assert '"QSGW_band_spin_"' in runner


def test_clean_driver_contains_no_hartree_or_hamiltonian_export_path() -> None:
    source = QSGW_DRIVER.read_text()
    forbidden = {
        "hartree_workflow.h",
        "hartree_route.h",
        "operator_fourier.h",
        "abacus_csr.h",
        "qsgw_update_hartree",
        "qsgw_export_hamiltonian_for_pyatb",
        "project_grid_operator_to_band(",
        "write_abacus_hamiltonian_csr(",
    }
    for token in forbidden:
        assert token not in source
