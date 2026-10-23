#!/usr/bin/env python3
"""Restore small ABACUS CSR export losses from matching dense full-mesh H/S.

Requires NumPy and the modern, single-ionic-step scalar ABACUS CSR format. Dense
hk*_nao.txt/sk*_nao.txt files and bz_sampling_out must belong to the same
complete, unshifted SCF mesh. H stays in ABACUS Ry units throughout.

Only IFT[dense(k) - FT(original(R))] is added: original real-space tails are
preserved. Even BvK boundaries are split equally between opposite images,
matching LibRPA's native QSGW reconstruction convention. This recovers mesh
precision; it cannot recover the discarded terms' unique off-mesh range.
The residual bound rejects incompatible inputs; no QSGW gate is relaxed.
Outputs are new files and never overwrite the original export.
"""

import argparse
import itertools
import json
from pathlib import Path
import re

import numpy as np


def read_csr(path):
    lines = path.read_text().splitlines()
    if sum('--- Ionic Step' in line for line in lines) != 1:
        raise ValueError('Expected one modern ABACUS ionic step')
    n = int(next(line.split()[0] for line in lines if 'number of localized basis' in line))
    nr = int(next(line.split()[0] for line in lines if 'number of Bravais lattice vector R' in line))
    start = next(i for i, line in enumerate(lines) if 'CSR Format' in line) + 1
    tokens = iter(re.findall(r'\([^()]*\)|[^\s]+', '\n'.join(
        line.split('#', 1)[0] for line in lines[start:])))
    blocks = {}
    for _ in range(nr):
        r = tuple(int(next(tokens)) for _ in range(3))
        nnz = int(next(tokens))
        values = []
        for _ in range(nnz):
            value = next(tokens)
            if value.startswith('('):
                real, imag = map(float, value[1:-1].split(','))
                values.append(complex(real, imag))
            else:
                values.append(complex(float(value)))
        columns = np.array([int(next(tokens)) for _ in range(nnz)], dtype=int)
        pointers = np.array([int(next(tokens)) for _ in range(n + 1)])
        if (r in blocks or pointers[0] != 0 or pointers[-1] != nnz
                or np.any(np.diff(pointers) < 0) or np.any(columns < 0)
                or np.any(columns >= n) or not np.all(np.isfinite(values))):
            raise ValueError('Invalid CSR block')
        matrix = np.zeros((n, n), dtype=complex)
        for row in range(n):
            col = columns[pointers[row]:pointers[row + 1]]
            if np.any(np.diff(col) <= 0):
                raise ValueError('Unordered or duplicate CSR columns')
            matrix[row, col] = values[pointers[row]:pointers[row + 1]]
        blocks[r] = matrix
    if next(tokens, None) is not None or not blocks:
        raise ValueError('Unexpected CSR payload')
    return blocks


def write_csr(path, blocks, real_values=False):
    n = next(iter(blocks.values())).shape[0]
    with path.open('x') as output:
        output.write(f'--- Ionic Step 1 ---\n{n} # number of localized basis\n'
                     f'{len(blocks)} # number of Bravais lattice vector R\n# CSR Format\n')
        for r, matrix in sorted(blocks.items()):
            rows, columns = np.nonzero(matrix)
            values = matrix[rows, columns]
            pointers = np.r_[0, np.cumsum(np.bincount(rows, minlength=n))]
            output.write(' '.join(map(str, (*r, len(values)))) + '\n# CSR values\n')
            output.write(' '.join(f'{x.real:.17e}' if real_values else
                                  f'({x.real:.17e},{x.imag:.17e})' for x in values) + '\n')
            output.write('# CSR column indices\n' + ' '.join(map(str, columns)) + '\n')
            output.write('# CSR row pointers\n' + ' '.join(map(str, pointers)) + '\n')


def read_dense(path, n):
    text = path.read_text()
    values = np.array([complex(float(a), float(b))
                       for a, b in re.findall(r'\(([^,()]+),([^,()]+)\)', text)])
    if values.size != n * (n + 1) // 2 or not np.all(np.isfinite(values)):
        raise ValueError(f'Expected finite triangular complex dense matrix: {path}')
    matrix = np.zeros((n, n), dtype=complex)
    matrix[np.triu_indices(n)] = values
    return matrix + np.triu(matrix, 1).conj().T


def fourier(blocks, kpoints):
    coordinates = np.array(list(blocks))
    matrices = np.array(list(blocks.values()))
    weights = np.exp(2j * np.pi * (kpoints @ coordinates.T))
    return (weights @ matrices.reshape(len(blocks), -1)).reshape(
        len(kpoints), *matrices.shape[1:])


def images(index, size):
    if size % 2 == 0 and index == size // 2:
        return (-size // 2, size // 2)
    return (index if index <= size // 2 else index - size,)


def restore(blocks, dense, kpoints, mesh, max_residual, real_values=False):
    sampled = fourier(blocks, kpoints)
    delta = dense - sampled
    residual = float(np.max(abs(delta)))
    if residual > max_residual:
        raise ValueError(f'Dense/CSR discrepancy {residual:.6e} exceeds {max_residual:.6e}')
    grid = np.empty((*mesh, *dense.shape[1:]), dtype=complex)
    indices = np.rint(kpoints * mesh).astype(int) % mesh
    grid[tuple(indices.T)] = delta
    correction = np.fft.fftn(grid, axes=(0, 1, 2)) / len(kpoints)
    result = {r: matrix.copy() for r, matrix in blocks.items()}
    for index in itertools.product(*(range(size) for size in mesh)):
        copies = list(itertools.product(*(images(i, size) for i, size in zip(index, mesh))))
        for r in copies:
            if r not in result:
                result[r] = np.zeros_like(dense[0])
            result[r] += correction[index] / len(copies)
    if real_values:
        if max(float(np.max(abs(matrix.imag))) for matrix in result.values()) > 1e-11:
            raise ValueError('Scalar ABACUS overlap requires real R-space entries')
        result = {r: matrix.real.copy() for r, matrix in result.items()}
    # Check the actual output operator, including scalar-S real projection.
    roundtrip = float(np.max(abs(fourier(result, kpoints) - dense)))
    hermiticity = max(float(np.max(abs(matrix - result.get(tuple(-x for x in r),
                                            np.zeros_like(matrix)).conj().T)))
                      for r, matrix in result.items())
    if roundtrip > 1e-11 or hermiticity > 1e-11:
        raise ValueError(f'Restored model failed roundtrip/Hermiticity: {roundtrip}, {hermiticity}')
    return result, dict(max_dense_minus_csr=residual, roundtrip=roundtrip,
                        real_space_hermiticity=hermiticity,
                        original_blocks=len(blocks), restored_blocks=len(result))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--hr', type=Path, required=True)
    parser.add_argument('--sr', type=Path, required=True)
    parser.add_argument('--dense-dir', type=Path, required=True)
    parser.add_argument('--bz', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--max-residual', type=float, default=1e-8,
                        help='Maximum allowed dense-minus-CSR entry, in native units (H: Ry)')
    args = parser.parse_args()
    if not np.isfinite(args.max_residual) or args.max_residual <= 0:
        parser.error('--max-residual must be positive and finite')
    lines = args.bz.read_text().splitlines()
    mesh = np.array(list(map(int, lines[0].split())))
    counts = list(map(int, lines[1].split()))
    sampling = np.loadtxt(lines[2:], ndmin=2)
    if (mesh.shape != (3,) or np.any(mesh <= 0) or sampling.shape[1] < 5
            or not np.all(np.isfinite(sampling))):
        raise ValueError('Requires a complete ordered unshifted full mesh')
    kpoints = sampling[:, 2:5]
    indices = np.rint(kpoints * mesh).astype(int) % mesh
    if (counts != [int(np.prod(mesh))] * 2
            or len(kpoints) != np.prod(mesh) or len(set(map(tuple, indices))) != len(kpoints)
            or np.max(abs(kpoints * mesh - np.rint(kpoints * mesh))) > 1e-6
            or not np.array_equal(sampling[:, 0], np.arange(1, len(kpoints) + 1))):
        raise ValueError('Requires a complete ordered unshifted full mesh')
    restored, report = {}, dict(mesh=mesh.tolist(), hamiltonian_units='Ry',
                                original_hr=str(args.hr.resolve()), original_sr=str(args.sr.resolve()),
                                dense_dir=str(args.dense_dir.resolve()),
                                max_allowed_residual=args.max_residual,
                                off_mesh_recovery='BvK residual only; original native tails retained')
    for kind, source, prefix in (('H', args.hr, 'hk'), ('S', args.sr, 'sk')):
        blocks = read_csr(source)
        n = next(iter(blocks.values())).shape[0]
        dense = np.array([read_dense(args.dense_dir / f'{prefix}{i}_nao.txt', n)
                          for i in range(1, len(kpoints) + 1)])
        restored[kind], report[kind] = restore(
            blocks, dense, kpoints, mesh, args.max_residual, real_values=(kind == 'S'))
    if next(iter(restored['H'].values())).shape != next(iter(restored['S'].values())).shape:
        raise ValueError('Hamiltonian and overlap dimensions differ')
    args.output_dir.mkdir(parents=True, exist_ok=False)
    for kind, source in (('H', args.hr), ('S', args.sr)):
        write_csr(args.output_dir / source.name, restored[kind], real_values=(kind == 'S'))
    (args.output_dir / 'restoration.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
