# Driver Usage

## Data preparation

To run the LibRPA driver for RPA/*GW* calculation, you need to export necessary
data from atomic-basis first-principles code.
The data include

- Atomistic structure
- one-particle orbital energies, occupation numbers and wave functions
- Bloch vectors used in periodic calculation
- Triple co-efficients of resolution of identity (RI)
- Coulomb matrix under the auxiliary basis for RI

### FHI-aims

Export of data required by LibRPA from FHI-aims is supported since the release version
231212 and also in the latest master branch. You can switch it on by adding to your `control.in` file:

```text
output librpa
# legacy switch before 240507 release
print_librpa_input .true.
```

It will dump necessary data files for many-body computation tasks.
Please refer to the FHI-aims manual for more details of the `output librpa` tag.

### ABACUS

Similar to the case with FHI-aims, the latest master branch of ABACUS supports the export of files required by LibRPA.
You only need to add the following line to the `INPUT` file of ABACUS:

```text
rpa 1
```

## LibRPA job

### Input file

Before triggering the LibRPA driver, an input file named `librpa.in` is required
at the working directory. An example is below:

```ini
task = rpa
input_dir = .
nfreq = 16
cs_threshold = 1e-4
```

The driver will read the parameters defined in the file and run the calculation
Please refer to the [guide page of runtime parameters](runtime_parameters)
for more information about the driver and API parameters.

By default, `input_preset = fhi-aims` preserves the historical driver filenames,
including `stru_out`, `basis_out`, `band_out`, `mommat_ks_kpt_*.dat`,
`Cs_data_*`, and `coulomb_mat_*` under `input_dir`. For ABACUS output, select:

```ini
input_preset = abacus
```

This changes `fn_stru`, `fn_eigocc_scf`, and `fn_vxc_scf` to `stru_out.txt`,
`band_out.txt`, and `vxc_out.txt`, and changes `prefix_velocity` from
`mommat_ks_kpt_` to `velocity_matrix`. Other filename and prefix defaults
remain unchanged.

If a host code exports the same data under different names, set the corresponding `fn_*` or `prefix_*` driver parameters in `librpa.in`:

```ini
fn_eigocc_scf = band_out_alt
prefix_lri_coeff = Cs_data_alt
prefix_coul_full = coulomb_mat_alt
prefix_coul_cut = coulomb_cut_alt
```

If `stru_out` contains symmetry operations, the driver builds a symmetry context from
the structure, basis metadata, and k-point grids. Runtime `use_symmetry_*` switches
then decide whether calculation paths use that context.

### Run the calculation

After setting up the input file `librpa.in`, the LibRPA driver can be called by issuing:

```bash
/path/to/LibRPA/build/chi0_main.exe
```

To run with multiple processes using MPI, you can simply invoke the relevant MPI driver,
for example

```bash
mpirun -np <nprocs> /path/to/LibRPA/build/chi0_main.exe
```

where `<nprocs>` is the number of MPI processes.
In addition, you want to run with multiple threads, you need to specify the
environment variable `OMP_NUM_THREADS`. For example

```bash
export OMP_NUM_THREADS=4
mpirun -np 4 /path/to/LibRPA/build/chi0_main.exe
```

This will run the LibRPA calculation using 4 MPI processes, each with 4 OpenMP threads.
