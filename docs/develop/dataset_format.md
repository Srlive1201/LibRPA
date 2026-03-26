# Dataset Format

This page documents the formats of the input data files required by the LibRPA driver.

## `stru_out`

The file `stru_out` contains structural information and k-point mesh data.

Its contents are arranged in the following order:

- **Lattice vectors**: 3 lines, each containing 3 floating-point numbers, in units of Bohr.
- **Reciprocal lattice vectors**: 3 lines, each containing 3 floating-point numbers, in units of Bohr${}^{-1}$.
- **Number of atoms**: 1 line containing the number of atoms in the unit cell, `n_atoms`.
- **Atomic coordinates and types**: `n_atoms` lines, each containing the Cartesian coordinates and the atom type of one atom.

The following entries describe the Brillouin-zone sampling.
This information is now stored in [`bz_sampling_out`](#bz-sampling-out) and is retained here for backward compatibility.

- **k-point grid dimensions**: 1 line containing three integers, `nkx`, `nky`, and `nkz`. The total number of k-points in the full grid is
  ```text
  nkpts = nkx * nky * nkz
  ```
- **Full k-point list**: `nkpts` lines, each containing 3 floating-point numbers giving the Cartesian coordinates of one k-point, in units of Bohr${}^{-1}$.
- **Mapping to irreducible k-points**: `nkpts` lines, each containing 1 integer.
  Suppose the integer on the n-th line is m.
  This means that the irreducible representative of the n-th k-point in the full k-point set is the m-th k-point in the full set.

## `basis_out`

This file describes the atomic basis sets used in the calculation,
including both the one-electron basis and the auxiliary basis.

Its structure is as follows.

The first line contains four entries:

1. total number of atom types, `n_atom_types`
2. total number of one-electron basis functions
3. total number of auxiliary basis functions
4. a string specifying the convention used for the ordering of azimuthal quantum numbers

For example:

```text
1        10        36    aims
```

The next `n_atom_types` lines provide a summary for each atom type. Each line contains:

1. atom type index
2. number of one-electron basis functions for this atom type
3. number of auxiliary basis functions for this atom type

For example:

```text
1         5        18
```

Next comes the description of the one-electron basis. There are `n_atom_types` blocks, one for each atom type.

In each block:

- the first line contains the atom type index and the number of radial functions
- the following lines list the angular momentum quantum number `l` for each radial function, one integer per line

For example:

```text
1       3
0
0
1
```

This means that atom type `1` has `5` radial functions in the one-electron basis, with angular momenta `0, 0, 0, 1, 1`.

After the one-electron basis blocks, the same block structure is repeated for the auxiliary basis.
Again, there are `n_atom_types` blocks. For each block:

- the first line contains the atom type index and the number of radial functions
- the following lines list the angular momentum quantum number `l` for each radial function in the auxiliary basis, one integer per line

For example:

```text
1      8
0
0
0
0
1
1
1
2
```

This means that atom type `1` has `8` radial functions in the auxiliary basis, with angular momenta `0, 0, 0, 0, 1, 1, 1, 2`.

(bz-sampling-out)=
## `bz_sampling_out`

This file describes the Brillouin-zone sampling used in the calculation,
including the full k-point grid and its reduction to the irreducible set.

Its structure is as follows.
The first line contains three integers:

1. `nk1`
2. `nk2`
3. `nk3`

These specify the number of k-point divisions along the three reciprocal lattice directions.
For example:

```text
3   3   3
```

The second line contains two integers:

1. total number of k-points in the full Brillouin-zone grid
2. number of k-points in the irreducible Brillouin zone

For example:

```text
27     14
```

This means that the full k-point grid contains `27` points, of which `14` belong to the irreducible set.

The next `n_k_points` lines describe all k-points in the full grid. Each line contains ten fields:

1. k-point index in the full set (1-based)
2. k-point weight
3. fractional coordinate `k1`
4. fractional coordinate `k2`
5. fractional coordinate `k3`
6. Cartesian coordinate `kx`
7. Cartesian coordinate `ky`
8. Cartesian coordinate `kz`
9. index of the corresponding irreducible k-point in the irreducible set
10. index of the corresponding irreducible k-point in the full k-point list

For example:

```text
2   0.37037037037E-01   0.00000000000E+00   0.00000000000E+00   0.33333333333E+00   0.00000000000E+00   0.00000000000E+00   0.35439508162E+00   2   2
```

This line indicates that full-grid k-point `2`

- has weight `0.037037037037`
- has fractional coordinates `(0, 0, 1/3)`
- has Cartesian coordinates `(0, 0, 0.35439508162)`
- maps to irreducible k-point `2`
- whose representative in the full k-point list is also point `2`

After the full k-point list, the file contains `n_irkpoints` lines summarizing the irreducible k-points. Each line contains three fields:

1. irreducible k-point index in the irreducible set (1-based)
2. index of its representative in the full k-point list
3. total weight of this irreducible k-point

For example:

```text
3      4   0.74074074074E-01
```

This means that irreducible k-point `3`

- is represented by full-grid k-point `4`
- carries total weight `0.074074074074`

A few remarks

- The weights in the full k-point list are the weights assigned to individual points in the full grid.
- The weights in the irreducible k-point summary are the accumulated weights of the corresponding symmetry-equivalent k-points.
- The representative index stored in the last field of the full k-point list can be used to identify which full-grid point serves as the representative of the irreducible class.

## `Cs_data_xxx.txt`

These files contain the localized RI triple coefficients.

In plain text format, each file has a header with two integers:
total number of atoms and number of periodic unit cells.
Then till the end of file, the data is formatted as blocks of RI coefficient $C$ on each pair of atoms and unit cell
```
i_atom_1  i_atom_2  n_1  n_2  n_3  n_basis_1  n_basis_2  n_aux_basis_1
C(1, 1, 1)
...
C(n_aux_basis_1, n_basis_2, n_basis_1)
```
Here `C` is the RI coefficients between the atom `i_atom_1` and `i_atom_2` in unit cells separated by
lattice vector $\mathbf{R} = n_1 \mathbf{a}_1 + n_2 \mathbf{a}_2 + n_3 \mathbf{a}_3$.
The auxiliary basis is located on `i_atom_1`. The number of basis functions on `i_atom_1` and `i_atom_2` 
are `n_basis_1` and `n_basis_2`, respectively. The number of auxiliary functions is `n_aux_basis_1`.
The indices of `C` runs in the Fortran order, i. e. the first index runs the fastest.

In binary format, the data is organized similarly in the plain text format, except for an extra integer
is included in the header, which is the number of atom pairs and lattice vectors included in the file.
The coefficients are saved in double precision. To better illustrate the format of binary file, the
following Python snippet could be helpful
```python
import struct
import numpy as np

# ensure that "Cs_data_0" exists and was generated with binary output mode in DFT code
cfile_path = "Cs_data_0.txt"

with open(cfile_path, 'rb') as h:
    n_atoms, n_cells, n_apcell_file = struct.unpack('iii', h.read(12))
    for _ in range(n_apcell_file):
        a1, a2, r1, r2, r3, nb1, nb2, nbb1 = struct.unpack('i' * 8, h.read(4 * 8))
        apcell = (a1, a2, r1, r2, r3)
        array_size = nb1 * nb2 * nbb1
        array = np.array(struct.unpack('d' * array_size, h.read(8 * array_size)))
        array = np.reshape(array, (nb1, nb2, nbb1))
        apcells[apcell] = array
```

## `band_out`

This file contains band energies and occupation numbers from the mean-field starting-point calculation.
It has a 5-line header
```
n_k_points
n_spins
n_states
n_basis
e_fermi
```
The first 4 lines contain an integer in each. The 5th line is a float number, which is the Fermi energy
in Hartree unit.

The remaining lines consists of `n_k_points*n_spins` blocks of `n_states+1` lines, in the format of
```
i_k_point    i_spin
1           f_1        e_1_ha      e_1_ev
2           f_2        e_2_ha      e_2_ev
3           f_3        e_3_ha      e_3_ev
...
n           f_n        e_n_ha      e_n_ev
...
```
This block contains the energies and occupation numbers of states $\left|\psi_{n,k\sigma}\right\rangle$
`i_k_point` marks the index of k-point $k$ in the full k-point set.
`i_spin` specify the spin channel $\sigma$.
In each of the following lines, the first integer species the index of state.
The 3 float numbers stand for the occupation number, the energy in Hartree unit and that in electronvolt
unit, respectively.
For spin-unpolarized calculation, `f_n` is a number from 0 to 2, otherwise it is from 0 to 1.

## `KS_eigenvector_xxx.txt`
These files contain the wave functions (eigenvectors) from the starting-point calculation expanded by orbital basis.
Each file can be divided in blocks of `n_states*n_basis*n_spins+1` lines,
where `n_states`, `n_basis` and `n_spins` will be extracted from `band_out` file.
Each block stores the data for a particular k-point, $c^i_{n,k\sigma}$:
```
i_k_point
c(1,1,1)_real c(1,1,1)_imag
...
c(i,n,s)_real c(i,n,s)_imag
...
```
The first line contains single integer, the index of the k-point of following data.
The remaining lines store the data with running index $i$, $n$, $\sigma$ in C-style row-major order,
i. e., spin index runs fastest, then state index and finally basis index.
Each line has two float numbers, which are the real and imaginary part of $c^i_{n,k\sigma}$.

## `coulomb_mat_xxx.txt`

These files contains the Coulomb matrices in auxiliary basis.
A single header line contains an integer, the number of irreducible k-point at
which the Coulomb matrices are computed.
The remaining part of the file is organized in blocks
```
n_aux_basis    row_start    row_end    col_start    col_end
i_k_point      k_weight
v(row_start, col_start  )_real       v(row_start, col_start  )_imag
v(row_start, col_start+1)_real       v(row_start, col_start+1)_imag
...
v(row_end, col_end)_real             v(row_end, col_end)_imag
```
where
- integer `n_aux_basis` is the total number of auxiliary basis functions.
- integer `row_start`, `row_end`, `col_start` and `col_end` mark the submatrix of the full Coulomb matrix
  that this block contain.
- integer `i_k_point` is the index of k-point of the current Coulomb matrix, in the full k-point list.
- float number `k_weight` is the weight of the irreducible k-points.

After the block header, there should be `(row_end-row_start+1)` times `(col_end-col_start+1)` lines
for the actual matrix element data. Each line contains two float numbers, which are the real and imaginary
parts of the element. The data is ordered in C-style row major.

## `coulomb_cut_xxx.txt`

There files are basically the same as `coulomb_mat_xxx.txt`, but store the truncated Coulomb to
be used in the GW calculation.

## `vxc_out`

The file `vxc_out` stores the exchange-correlation potential for electronic states on the SCF k-point grid.

The header consists of three lines:

```text
n_k_points
n_spins
n_states
```

The header is followed by n_k_points * n_spins * n_states data lines.
Each line contains two columns and corresponds to one state identified by the tuple `(i_k, i_spin, i_state)`.
The data are ordered such that `i_state` runs fastest, followed by `i_spin`, and then `i_k`. In other words, the lines are arranged as

```
# i_k       i_spin      i_state
    0            0            0
    0            0            1
    0            0            2
...
    0            0   n_states-1
    0            1            0
...
    1            0            0
...
```

The two columns contain the same exchange-correlation potential,
first in Hartree unit while the second in eV.

## Input files for band structure calculation

For band-structure calculations, LibRPA reads the following input files:

- `band_kpath_info`
- `band_KS_eigenvalue_k_{ik:05d}.txt`
- `band_KS_eigenvector_k_{ik:05d}.txt`
- `band_vxc_k_{ik:05d}.txt`

Here `ik` is the 1-based index of the k-point along the band path, written with five digits.

### `band_kpath_info`

The file `band_kpath_info` defines the k-point path used for the band-structure calculation.

The first line contains four integers:

1. number of basis functions
2. number of states
3. number of spin channels
4. number of k-points on the band path

For example:

```text
18    18     1    10
```

The remaining `n_kpath_points` lines each contain three floating-point numbers, giving the fractional coordinates of one k-point on the band path.

For example:

```text
0.500000000000000000E+00   0.500000000000000000E+00   0.500000000000000000E+00
```

Each such line represents one k-point in fractional reciprocal coordinates.

### `band_KS_eigenvalue_k_{ik:05d}.txt`

For each k-point on the band path, the file `band_KS_eigenvalue_k_{ik:05d}.txt` stores the Kohn-Sham eigenvalues used by LibRPA.
Each line corresponds to one state and contains five columns:

1. spin index
2. state index
3. occupation number
4. Kohn-Sham eigenvalue in Hartree
5. Kohn-Sham eigenvalue in eV

For example:

```text
1       3   0.200000000000000000E+01  -0.658097773108510893E+02  -0.179077515426494506E+04
```

This line indicates that, at the selected k-point,

- the spin index is `1`
- the state index is `3`
- the occupation number is `2.0`
- the eigenvalue is given both in Hartree and in eV

The data are ordered such that state index `i_state` runs fastest and followed by spin index `i_spin`.

### `band_vxc_k_{ik:05d}.txt`

For each k-point on the band path, the file `band_vxc_k_{ik:05d}.txt` stores the diagonal matrix elements of the exchange-correlation potential for the corresponding Kohn-Sham states.

Each line contains three columns:

1. spin index
2. state index
3. exchange-correlation potential in Hartree

For example:

```text
1       3  -0.562542321738239171E+01
```

This line gives the exchange-correlation potential for state `3` in spin channel `1` at the selected k-point.
The data are ordered such that state index `i_state` runs fastest and followed by spin index `i_spin`.

### `band_KS_eigenvector_k_{ik:05d}.txt`

For each k-point on the band path, the file `band_KS_eigenvector_k_{ik:05d}.txt` stores the Kohn-Sham eigenvectors at that k-point.
The file contains a complex array of shape `(n_spins, n_states, n_basis)`,
written in binary format using C-style ordering.
