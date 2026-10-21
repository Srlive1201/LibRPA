#!/usr/bin/env python3
"""Compare two native current-HR QSGW/cRPA updates from full and IBZ SOC inputs.

Reuse the existing full-BZ and paired-IBZ fixtures, with all magnetic symmetry
routes enabled, fractional endpoint weights and parent tracking. Compare the
complete optical/self-energy/Hamiltonian matrices in common eigenframes as
well as QP roots and local interactions. This is a finite-model control.
"""

import argparse
import json
from pathlib import Path
import re
import tarfile

import numpy as np

from verify_crpa_workflows import HA_TO_EV, compare_interactions, crpa_states
from verify_qsgw_headwing import TRACE_DTYPE
from verify_qsgw_restart import quasiparticles, run, stage_input


def coefficients(directory):
    result = {}
    for path in sorted((directory / "librpa.d/crpa").glob("state_*")):
        with np.load(path / "crpa_model.npz") as archive:
            result[int(path.name.removeprefix("state_"))] = archive["C"].copy()
    return result


def matrices(directory):
    records = np.loadtxt(directory / "librpa.d/qsgw_matrices.dat", dtype=TRACE_DTYPE)
    keys = ["iteration", "channel", "component", "spin", "k", "frequency_index"]
    result = {}
    for key in set(zip(*(records[name] for name in keys))):
        selected = np.ones(records.shape, dtype=bool)
        for name, value in zip(keys, key):
            selected &= records[name] == value
        block = records[selected]
        shape = int(block["row"].max()) + 1, int(block["column"].max()) + 1
        assert len(block) == shape[0] * shape[1], key
        value = np.full(shape, np.nan + 1j*np.nan)
        value[block["row"], block["column"]] = block["real"] + 1j*block["imag"]
        assert np.all(np.isfinite(value)), key
        result[key] = value
    return result


def verify_feedback(values, states):
    report = {}
    for component in ("head_tensor", "wing_auxiliary", "head_velocity_0",
                      "head_velocity_1", "head_velocity_2"):
        magnitude, change = 0.0, 0.0
        for key, value in values.items():
            if key[0] != 1 or key[2] != component:
                continue
            later = values[(2, *key[1:])]
            if component.startswith("head_velocity_"):
                # Compare the changed operator in a common entering eigenframe,
                # rather than counting arbitrary phases as velocity feedback.
                transform = np.linalg.solve(states[0][key[3], key[4]],
                                            states[1][key[3], key[4]])
                later = transform @ later @ transform.conj().T
            response = value - np.eye(len(value)) if component == "head_tensor" else value
            magnitude = max(magnitude, float(np.max(abs(response))))
            change = max(change, float(np.max(abs(later - value))))
        # Same nontrivial-feedback gate as verify_qsgw_headwing.py.
        assert magnitude > 1e-12 and change > 1e-12, (component, magnitude, change)
        report[component] = {"magnitude": magnitude, "change": change}
    return report


def compare(reference, candidate):
    a, b = crpa_states(reference), crpa_states(candidate)
    assert set(a) == set(b) == {0, 1, 2}
    for states in (a, b):
        assert np.max(abs(states[2]["A_weights"] - states[0]["A_weights"])) > 1e-6
    report = {"local_interactions": {state: compare_interactions(a[state], b[state])
                                     for state in a}}
    ca, cb = coefficients(reference), coefficients(candidate)
    frame = {state: np.linalg.solve(ca[state], cb[state]) for state in ca}
    for matrices_at_state in frame.values():
        for spin in matrices_at_state:
            for matrix in spin:
                np.testing.assert_allclose(matrix.conj().T @ matrix, np.eye(len(matrix)),
                                           rtol=0, atol=1e-10)

    qa = quasiparticles(reference / "librpa.d", 1, 2)
    qb = quasiparticles(candidate / "librpa.d", 1, 2)
    assert qa.keys() == qb.keys()
    qp_error = 0.0
    for key in qa:
        assert qa[key][0] == qb[key][0], key
        np.testing.assert_allclose(qa[key][1], qb[key][1], rtol=0, atol=1e-7)
        qp_error = max(qp_error, float(np.max(abs(qa[key][1] - qb[key][1]))))
    report.update(qp_rows=len(qa), maximum_qp_difference_ev=qp_error)

    ma, mb = matrices(reference), matrices(candidate)
    report["reference_feedback"] = verify_feedback(ma, ca)
    report["candidate_feedback"] = verify_feedback(mb, cb)
    wanted = lambda key: key[2] != "rotation_u" and not key[2].startswith("wfc")
    keys = {key for key in ma if wanted(key)}
    assert keys == {key for key in mb if wanted(key)}
    fixed_basis = {"exx", "sigma_c_iw", "raw_h", "mixed_h", "vc", "vxc_dft", "h0"}
    errors = {}
    for key in keys:
        iteration, channel, component, spin, k, frequency = key
        value = mb[key]
        if component in fixed_basis or component.startswith("head_velocity_"):
            state = iteration - 1 if component.startswith("head_velocity_") else 0
            transform = frame[state][spin, k]
            value = transform @ value @ transform.conj().T
        assert value.shape == ma[key].shape, key
        error = float(np.max(abs(value - ma[key])))
        assert error < 1e-9, (key, error)
        errors[component] = max(errors.get(component, 0.0), error)
    assert {"head_tensor", "wing_auxiliary", "sigma_c_iw", "mixed_h",
            "head_velocity_0", "head_velocity_1", "head_velocity_2"} <= errors.keys()
    report["aligned_matrix_errors_au"] = errors
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--mpiexec", default="mpirun")
    parser.add_argument("--ranks", type=int, action="append")
    args = parser.parse_args()
    root = args.work_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    fixtures = Path(__file__).resolve().parent / "testcases"
    controls = {
        "task": "qsgw", "qsgw_min_iter": "2", "qsgw_max_iter": "2",
        "qsgw_cut_mode": "0", "qsgw_headwing_velocity": "current_hr",
        "qsgw_write_iteration_matrices": "true", "qsgw_output_qp": "true",
        "use_fullcoul_eps": "true", "use_fullcoul_wc": "true",
        "tfgrids_freq_min": "0.4", "tfgrids_freq_max": "4.0",
        "use_symmetry_exx": "true", "use_symmetry_rpa": "true", "use_symmetry_gw": "true",
        "crpa_enabled": "true", "crpa_species_labels": "X", "crpa_correlated_species": "X",
        "crpa_parent_orbitals": "s", "crpa_output_orbitals": "s",
        "crpa_overlap_file": "srs1_nao.csr", "crpa_orbital_bands": "0 1 2 3",
        "crpa_response_mode": "orbital_weights", "crpa_track_orbitals": "true",
        "crpa_output_analysis": "false", "crpa_output_model": "true",
        "crpa_output_raw": "true", "crpa_output_static": "false", "crpa_output_format": "npz",
    }
    datasets, templates = {}, {}
    for name, fixture in (("full", "g0w0_soc_full_native_headwing_libri"),
                          ("ibz", "g0w0_soc_ibz_headwing_libri")):
        home = root / name
        home.mkdir()
        with tarfile.open(fixtures / fixture / "dataset.tar.gz") as archive:
            archive.extractall(home, filter="data")
        datasets[name] = home / "dataset"
        templates[name] = fixtures / fixture / "librpa/librpa.in"
        # Native QSGW consumes Ry; convert only the extracted fixture's eV H.
        hamiltonian = datasets[name] / "hrs1_nao.csr"
        hamiltonian.write_text(re.sub(
            r"\(([-+0-9.eE]+),([-+0-9.eE]+)\)",
            lambda match: f"({float(match[1])*2/HA_TO_EV:.17e},{float(match[2])*2/HA_TO_EV:.17e})",
            hamiltonian.read_text()))
    report, reference = {}, None
    for ranks in args.ranks or (1, 4):
        runs = {}
        for name in ("full", "ibz"):
            runs[name] = root / name / f"np{ranks}"
            stage_input(templates[name], datasets[name], runs[name], controls)
            run(args.binary.resolve(), args.mpiexec, ranks, runs[name])
        report[ranks] = {"ibz_full": compare(runs["full"], runs["ibz"])}
        if reference is not None:
            report[ranks]["rank_comparison"] = compare(reference, runs["full"])
        else:
            reference = runs["full"]
        print(json.dumps({ranks: report[ranks]}, indent=2), flush=True)
    (root / "comparison.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
