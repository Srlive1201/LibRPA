#!/usr/bin/env python3
"""Verify three live head/wing updates against an independent split restart.

The bundled aims-rotate case needs no extra data. Native current-HR verification
requires a full-grid ABACUS Si dataset containing full Vxc, S(R), r(R), RI and
Coulomb inputs. Numerical trajectories do not establish material convergence.
"""

import argparse
import json
from pathlib import Path
import tarfile

import numpy as np

from verify_qsgw_restart import compare, run, stage_input


TRACE_DTYPE = [
    ("iteration", "i4"), ("channel", "i4"), ("component", "U32"),
    ("spin", "i4"), ("k", "i4"), ("frequency_index", "i4"),
    ("frequency", "f8"), ("row", "i4"), ("column", "i4"),
    ("real", "f8"), ("imag", "f8"),
]


def compare_feedback(continuous, resumed, initial):
    a = np.loadtxt(continuous / "librpa.d/qsgw_matrices.dat", dtype=TRACE_DTYPE)
    b = np.loadtxt(resumed / "librpa.d/qsgw_matrices.dat", dtype=TRACE_DTYPE)
    boundary = np.loadtxt(initial / "librpa.d/qsgw_matrices.dat", dtype=TRACE_DTYPE)
    errors = {}
    raw_velocity_errors = {}
    frame_error = 0.0
    changes = {}

    def block(data, iteration, name, channel, spin, k):
        rows = data[(data["iteration"] == iteration) & (data["component"] == name)
                    & (data["channel"] == channel) & (data["spin"] == spin)
                    & (data["k"] == k)]
        assert len(rows), (iteration, name, channel, spin, k)
        n = int(max(rows["row"].max(), rows["column"].max())) + 1
        assert len(rows) == n*n
        matrix = np.full((n, n), np.nan + 1j*np.nan)
        matrix[rows["row"], rows["column"]] = rows["real"] + 1j*rows["imag"]
        return matrix

    components = ["head_tensor", "wing_auxiliary",
                  "head_velocity_0", "head_velocity_1", "head_velocity_2"]
    for name in components:
        for iteration in (2, 3):
            x = a[(a["iteration"] == iteration) & (a["component"] == name)]
            y = b[(b["iteration"] == iteration) & (b["component"] == name)]
            assert len(x) == len(y) > 0, (name, iteration)
            for key, _ in TRACE_DTYPE[:9]:
                if key == "frequency":
                    # Nodes are regenerated from the live spectrum; MPI
                    # roundoff also perturbs their frequency scale.
                    np.testing.assert_allclose(x[key], y[key], rtol=1e-12, atol=1e-12)
                else:
                    np.testing.assert_array_equal(x[key], y[key])
            error = np.max(abs(x["real"] + 1j*x["imag"] - y["real"] - 1j*y["imag"]))
            if name.startswith("head_velocity_"):
                raw_velocity_errors[f"{iteration}:{name}"] = float(error)
                error = 0.0
                for channel, spin, k in set(zip(x["channel"], x["spin"], x["k"])):
                    # V enters update n in the eigenbasis produced by n-1.
                    # Compare it through the independently saved rotations:
                    # G = U_cont^dagger U_restart, V_restart = G^dagger V_cont G.
                    u = block(a, iteration-1, "rotation_u", channel, spin, k)
                    v = block(boundary if iteration == 2 else b, iteration-1,
                              "rotation_u", channel, spin, k)
                    frame = u.conj().T @ v
                    residual = np.max(abs(frame.conj().T @ frame - np.eye(len(frame))))
                    frame_error = max(frame_error, float(residual))
                    assert residual < 1e-10
                    va = block(a, iteration, name, channel, spin, k)
                    vb = block(b, iteration, name, channel, spin, k)
                    error = max(error, float(np.max(abs(frame.conj().T @ va @ frame - vb))))
            assert error < 1e-8, (name, iteration, error)
            errors[f"{iteration}:{name}"] = float(error)
        x = a[(a["iteration"] == 1) & (a["component"] == name)]
        y = a[(a["iteration"] == 2) & (a["component"] == name)]
        assert len(x) == len(y) > 0, name
        changes[name] = float(np.max(abs(x["real"] + 1j*x["imag"] -
                                        y["real"] - 1j*y["imag"])))
        assert changes[name] > 1e-12, (name, "feedback was not exercised")
    frequency_checks = []
    for iteration in (1, 2, 3):
        head = a[(a["iteration"] == iteration) & (a["component"] == "head_tensor")]
        wing = a[(a["iteration"] == iteration) & (a["component"] == "wing_auxiliary")]
        body = a[(a["iteration"] == iteration) & (a["component"] == "sigma_c_iw")]
        head_nodes, wing_nodes, body_nodes = (
            np.unique(component["frequency"]) for component in (head, wing, body)
        )
        assert len(head_nodes) == len(wing_nodes) == len(body_nodes) > 0, iteration
        np.testing.assert_allclose(head_nodes, body_nodes, rtol=0, atol=1e-12)
        np.testing.assert_allclose(wing_nodes, body_nodes, rtol=0, atol=1e-12)
        frequency_checks.append({
            "iteration": iteration,
            "head_nodes_Ha": head_nodes.tolist(),
            "wing_nodes_Ha": wing_nodes.tolist(),
            "body_nodes_Ha": body_nodes.tolist(),
            "maximum_node_difference_Ha": float(max(
                np.max(abs(head_nodes - body_nodes)),
                np.max(abs(wing_nodes - body_nodes)),
            )),
        })
    return {"restart_component_errors": errors, "raw_velocity_errors": raw_velocity_errors,
            "maximum_frame_unitarity_error": frame_error,
            "first_to_second_update_changes": changes,
            "frequency_node_checks": frequency_checks}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--case", choices=("aims-rotate", "native-current"), required=True)
    parser.add_argument("--dataset", type=Path)
    parser.add_argument("--ranks", type=int, default=2)
    parser.add_argument("--mpiexec", default="mpirun")
    args = parser.parse_args()
    if args.case == "native-current" and args.dataset is None:
        parser.error("native-current requires --dataset with full-grid native S/r inputs")
    root = args.work_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    source = Path(__file__).resolve().parent
    if args.case == "aims-rotate":
        fixture = source / "testcases/qsgw_aims_Si_k333_headonly_libri"
        with tarfile.open(fixture / "dataset.tar.gz") as archive:
            archive.extractall(root, filter="data")
        dataset = root / "dataset"
        template = fixture / "librpa/librpa.in"
        controls = {"qsgw_headwing_velocity": "rotate_reference"}
    else:
        dataset = args.dataset.resolve()
        template = source / "manual/Si_k444_live_headwing/librpa.in"
        controls = {"qsgw_headwing_velocity": "current_hr"}
    controls.update(option_dielect_func="3", replace_w_head="true",
                    qsgw_write_iteration_matrices="true",
                    qsgw_convergence_tolerance_ev="1e-15")
    continuous, initial, resumed = (root / name for name in ("continuous", "initial", "resumed"))
    for directory, maximum in ((continuous, 3), (initial, 1), (resumed, 3)):
        settings = {**controls, "qsgw_max_iter": str(maximum), "qsgw_min_iter": "1"}
        if directory == resumed:
            settings["qsgw_restart_file"] = str(initial / "librpa.d/qsgw_checkpoint.bin")
        stage_input(template, dataset, directory, settings)
        run(args.binary.resolve(), args.mpiexec, args.ranks, directory)
    report = {
        "case": args.case,
        "first_update": compare(continuous, initial, 1, 1),
        "restored_boundary": compare(initial, resumed, 1, 1, exact_boundary=True),
        "continued_updates": compare(continuous, resumed, 2, 3),
        "feedback": compare_feedback(continuous, resumed, initial),
    }
    (root / "comparison.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
