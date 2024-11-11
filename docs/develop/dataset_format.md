# Dataset Format

## `stru_out`

This file contains structure and k-mesh data.
It should contain the following data in order
- lattice vectors: 3 lines, 3 float number in each line, unit: Bohr radius
- reciprocal lattice vectors: 3 lines, 3 float number in each line, unit: inverse Bohr radius
- number of k-grids along each lattice vectors: 1 line, `nkx`, `nky`, `nkz`. The
  total number of k-points `nkpts` equals to the product of `nkx`, `nky` and
  `nkz`
- Cartesian coordinates of each k-point: `nkpts` lines, 3 float number in each line, unit: inverse Bohr radius
- mapping of k-point to its irreducible counterpart: `nkpts` lines, 1 integer in each line.

The mapping should be considered as below:
suppose the number on the `n`-th line is `m`, it means that
the irreducible k-point corresponding to the `n`-th k-point in the full k-point set is the `m`-th
k-point in the full set.

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
