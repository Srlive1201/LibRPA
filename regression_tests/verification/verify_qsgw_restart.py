#!/usr/bin/env python3
"""Compare complete QSGW trajectories with iteration-boundary restarts.

Uses the bundled Si k333/head and H2O band/cut numerical cases. The original
datasets and reference tables are never modified. Requires NumPy and the same
MPI/compiler runtime environment as the supplied LibRPA executable. The fixed
mixers split three updates after one.
"""

import argparse
import csv
import json
import os
from pathlib import Path
import re
import subprocess
import tarfile

import numpy as np


def quasiparticles(directory, first_update, last_update):
    """Check entering-state energies and the QP equation, then return its roots."""
    with (directory / "qsgw_qp.csv").open() as source:
        rows = [row for row in csv.DictReader(source)
                if first_update <= int(row["gw_update"]) <= last_update]
    assert rows, (directory, first_update, last_update)
    eigenvalues = np.atleast_2d(np.loadtxt(directory / "qsgw_eigenvalues.dat"))
    energies = {tuple(int(row[i]) for i in (0, 1, 2, 3, 7)): row[8]
                for row in eigenvalues}
    roots = {}
    for row in rows:
        key = (int(row["gw_update"]), row["channel"], int(row["spin"]),
               int(row["kpoint"]), int(row["band"]))
        assert key not in roots, key
        values = np.array([float(row[name]) for name in (
            "occupation", "e0_eV", "vstatic_eV", "sigma_x_eV",
            "re_sigma_c_effective_eV", "im_sigma_c_eV", "eqp_eV")])
        assert np.all(np.isfinite(values)), key
        state = int(row["state_index"])
        eigenkey = (state, 0 if row["channel"] == "grid" else 1,
                    key[2] - 1, key[3] - 1, key[4] - 1)
        np.testing.assert_allclose(values[1], energies[eigenkey], rtol=0, atol=1e-9)
        np.testing.assert_allclose(values[6], values[1] - values[2] + values[3] + values[4],
                                   rtol=0, atol=1e-9)
        roots[key] = (int(row["solver_info"]), values)

    # The ordinary G0W0 writers must describe the same physical roots. Their
    # tolerances follow the printed precision, rather than the CSV precision.
    controls = (directory.parent / "librpa.in").read_text()
    aims_units = re.search(r"(?m)^\s*constants_choice\s*=\s*aims\s*$", controls)
    ha2ev = 27.2113845 if aims_units else 27.211386245988
    for update in range(first_update, last_update + 1):
        grid = sorted((key, value[1]) for key, value in roots.items()
                      if key[0] == update and key[1] == "grid")
        grid.sort(key=lambda item: (item[0][3], item[0][2], item[0][4]))
        table = []
        for line in (directory / f"energy_qp_{update}").read_text().splitlines():
            fields = line.split()
            if len(fields) == 4 and fields[0].isdigit():
                table.append([float(field) for field in fields])
        expected = np.array([[key[4], value[0], value[1] / ha2ev, value[6] / ha2ev]
                             for key, value in grid])
        table = np.array(table)
        np.testing.assert_allclose(table[:, :2], expected[:, :2], rtol=0, atol=5e-5)
        np.testing.assert_allclose(table[:, 2:], expected[:, 2:], rtol=0, atol=1e-8)
        for spin in sorted({key[2] for key in roots if key[0] == update and key[1] == "band"}):
            band = sorted((key, value[1]) for key, value in roots.items()
                          if key[:3] == (update, "band", spin))
            table = np.atleast_2d(np.loadtxt(directory / f"GW_band_spin_{spin}_{update}.dat"))
            nbands = len({key[4] for key, _ in band})
            expected = np.array([[value[0], value[6]] for _, value in band]).reshape(-1, 2 * nbands)
            np.testing.assert_allclose(table[:, 4:], expected, rtol=0, atol=5.1e-6)
    return roots


def stage_input(template, dataset, run_dir, overrides):
    run_dir.mkdir(parents=True, exist_ok=False)
    text = template.read_text()
    parameters = {"input_dir": str(dataset), **overrides}
    for name, value in parameters.items():
        pattern = rf"(?m)^\s*{re.escape(name)}\s*=.*$"
        text = re.sub(pattern, f"{name} = {value}", text)
        if not re.search(rf"(?m)^\s*{re.escape(name)}\s*=", text):
            text += f"\n{name} = {value}\n"
    (run_dir / "librpa.in").write_text(text)


def run(binary, mpiexec, ranks, directory):
    print(f"Running {directory.name}: {ranks} MPI rank(s)", flush=True)
    environment = dict(os.environ, OMP_NUM_THREADS="1", MKL_NUM_THREADS="1")
    with (directory / "run.log").open("w") as output:
        subprocess.run([mpiexec, "-n", str(ranks), str(binary)], cwd=directory,
                       env=environment, stdout=output, stderr=subprocess.STDOUT,
                       check=True)


def matrices(filename, minimum_iteration, maximum_iteration):
    # Eigenvector phases and rotations inside degenerate eigenspaces can vary
    # with MPI reductions. Compare reference-basis operators for the propagated
    # trajectory; the restored boundary additionally compares WFCs exactly.
    result = {}
    with filename.open() as source:
        for line in source:
            if line.startswith("#"):
                continue
            fields = line.split()
            iteration = int(fields[0])
            if minimum_iteration <= iteration <= maximum_iteration:
                key = tuple(fields[i] for i in (0, 1, 2, 3, 4, 5, 7, 8))
                result[key] = np.array([float(fields[i]) for i in (6, 9, 10)])
    return result


def compare(reference, restarted, first_iteration, last_iteration,
            exact_boundary=False):
    reference = reference / "librpa.d"
    restarted = restarted / "librpa.d"
    metrics = {}
    for name, columns, tolerance in (
        ("eigenvalues", slice(8, 9), 1e-7),
        ("iterations", slice(1, 9), 1e-7),
    ):
        a = np.atleast_2d(np.loadtxt(reference / f"qsgw_{name}.dat"))
        b = np.atleast_2d(np.loadtxt(restarted / f"qsgw_{name}.dat"))
        a = a[(a[:, 0] >= first_iteration) & (a[:, 0] <= last_iteration)]
        b = b[(b[:, 0] >= first_iteration) & (b[:, 0] <= last_iteration)]
        assert a.shape == b.shape, (name, a.shape, b.shape)
        if name == "eigenvalues":
            np.testing.assert_allclose(a[:, :8], b[:, :8], rtol=0, atol=1e-12)
        else:
            np.testing.assert_array_equal(a[:, 0], b[:, 0])
        delta = float(np.max(np.abs(a[:, columns] - b[:, columns])))
        assert delta <= (0.0 if exact_boundary else tolerance), (name, delta)
        metrics[f"maximum_{name}_difference"] = delta
        metrics[f"{name}_rows"] = len(a)

    a = matrices(reference / "qsgw_matrices.dat", first_iteration, last_iteration)
    b = matrices(restarted / "qsgw_matrices.dat", first_iteration, last_iteration)
    # A restored boundary has mixed H, WFC and occupations; it does not recompute
    # the previous Sigma, Vc or diagonalization rotation. Entering-eigenframe
    # velocities need a basis alignment, checked by the head/wing verifier.
    if exact_boundary:
        keys = set(b)
        assert keys <= set(a)
    else:
        def fixed_basis(key):
            return (key[2] != "rotation_u" and not key[2].startswith("wfc")
                    and not key[2].startswith("head_velocity_"))
        keys = {key for key in a if fixed_basis(key)}
        assert keys == {key for key in b if fixed_basis(key)}
    delta = max(float(np.max(np.abs(a[key] - b[key]))) for key in keys)
    assert delta <= (0.0 if exact_boundary else 1e-8), delta
    metrics["maximum_matrix_difference_au"] = delta
    metrics["matrix_rows"] = len(keys)
    if not exact_boundary:
        qa = quasiparticles(reference, first_iteration, last_iteration)
        qb = quasiparticles(restarted, first_iteration, last_iteration)
        assert qa.keys() == qb.keys()
        qp_delta = 0.0
        for key in qa:
            assert qa[key][0] == qb[key][0], key
            np.testing.assert_allclose(qa[key][1][0], qb[key][1][0], rtol=0, atol=1e-12)
            delta = float(np.max(abs(qa[key][1][1:] - qb[key][1][1:])))
            assert delta <= 1e-7, (key, delta)
            qp_delta = max(qp_delta, delta)
        metrics["maximum_qp_difference_eV"] = qp_delta
        metrics["qp_rows"] = len(qa)
        for filename in reference.glob(f"*band_spin_*_{last_iteration}.dat"):
            x = np.atleast_2d(np.loadtxt(filename))
            y = np.atleast_2d(np.loadtxt(restarted / filename.name))
            np.testing.assert_allclose(x, y, rtol=0, atol=1e-5)
    return metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", type=Path, required=True)
    parser.add_argument("--work-dir", type=Path, required=True)
    parser.add_argument("--mpiexec", default="mpirun")
    parser.add_argument("--ranks", type=int, default=4)
    parser.add_argument("--case", choices=("si-none", "h2o-band"),
                        action="append")
    args = parser.parse_args()
    binary = args.binary.resolve()
    root = args.work_dir.resolve()
    root.mkdir(parents=True, exist_ok=False)
    testcases = Path(__file__).resolve().parents[1] / "testcases"
    cases = {
        "si-none": ("qsgw_aims_Si_k333_headonly_libri", {}),
        "h2o-band": ("qsgw_band_aims_mole_H2O_libri",
                     {"qsgw_mixing_beta": "0.3", "qsgw_unoccupied_keep": "1",
                      "qsgw_cut_shift_ha": "0.05"}),
    }
    report = {}
    for case in args.case or cases:
        fixture, controls = cases[case]
        home = root / case
        home.mkdir()
        with tarfile.open(testcases / fixture / "dataset.tar.gz") as archive:
            archive.extractall(home, filter="data")
        template = testcases / fixture / "librpa" / "librpa.in"
        common = {"qsgw_write_iteration_matrices": "true", **controls}
        total, split = 3, 1
        continuous, initial, resumed = (home / name for name in ("continuous", "initial", "resumed"))
        for directory, end in ((continuous, total), (initial, split), (resumed, total)):
            overrides = {**common, "qsgw_min_iter": str(end), "qsgw_max_iter": str(end)}
            if directory == resumed:
                overrides["qsgw_restart_file"] = str(initial / "librpa.d" / "qsgw_checkpoint.bin")
            stage_input(template, home / "dataset", directory, overrides)
            run(binary, args.mpiexec, args.ranks, directory)
        report[case] = {
            "mpi_ranks": args.ranks,
            "first_update": compare(continuous, initial, 1, split),
            "restored_boundary": compare(initial, resumed, split, split, exact_boundary=True),
            "continued_updates": compare(continuous, resumed, split + 1, total),
        }
        print(json.dumps({case: report[case]}, indent=2), flush=True)
    (root / "comparison.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
