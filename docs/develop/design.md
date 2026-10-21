# Code Design

LibRPA is designed to provide efficient and scalable implementations of various beyond-DFT methods, including but not limited to random-phase approximation (RPA) correlation energies. The library is developed in C++ and Python, utilizing a hybrid parallelism model with MPI and OpenMP to maximize performance on modern high-performance computing systems.

## Framework of RPA calculation

![Framework of RPA calculation](../assets/librpa_flowchat.jpg)

## RPA interface to FHI-aims

![RPA interface to FHI-aims](../assets/aims_interface.drawio.jpg)

## RPA interface to ABACUS

![RPA interface to ABACUS](../assets/abacus_interface.drawio.jpg)

## Numerical tasks and the driver

G0W0 and QSGW are library tasks. Their numerical workflows live in
`src/core/g0w0_task` and `src/core/qsgw/task`.
The core `Dataset` owns the mean fields, response/screening state and parallel
layouts; the API registry resolves public handlers to those states. Numerical
tasks assemble lower-level math, MPI and I/O services. Those lower modules do
not include core physics headers. Atomic-basis adapters and physics-specific
output assembly live in core; generic serializers accept value buffers.

`driver/tasks` parses host input, invokes the public C/C++ APIs and writes
requested output. Its QSGW entry point submits numerical Vxc and native H/S/r
blocks to an opaque public task. The library owns iteration, fixed/live frames,
convergence, state invalidation and checkpoints. A rank-zero output callback
receives borrowed text or localized-operator views for presentation; it must
not mutate or recursively call the task/handler. The handler must remain alive until task deletion. After preparation keep
the SCF eigenframe/grid fixed, finish band/Coulomb/velocity input before running,
and do not mutate or reenter the handler during a run. Ordinary input filenames never enter
the numerical task. The task may read its own library-generated checkpoint.

## QSGW iteration state and checkpoints

Checkpoint state, serialization and input comparison live in
`src/core/qsgw/checkpoint`. For filenames, switches and continuation examples,
see the [QSGW user guide](../user_guide/qsgw.md#checkpoint-continuation).
The source-only note `theory/qsgw/staticization.md` describes the self-energy map.

### Reference state and saved payload

Keep the original DFT reference immutable throughout a run and its
continuations. Current Hamiltonians are represented in that initial KS
basis; current eigenvalues and wavefunctions describe the updated states.
Restoring the saved wavefunctions directly preserves their phases and
degenerate-subspace frame at the restart boundary.

| Payload | Contents and purpose |
| --- | --- |
| Input reference | Grid reference eigenpairs, occupation weights, Fermi energy and DFT Vxc; grid/path k coordinates and grid weights. These bind the saved state to the initial basis and input. |
| Live state | Current reference-basis Hamiltonian, eigenvalues, complex wavefunctions, occupation weights and Fermi energy. Dimensions include spin, spinor, k point, band and AO counts. |
| Iteration summary | Completed iteration number, maximum eigenvalue change, Hamiltonian residual norms, Fermi energy, gap, electron count, effective mixing fraction and convergence flag. |
| Update settings | Cut mode, retained unoccupied count, shift, effective mixing fraction, frequency count and head mode. |
| Optional band channel | Band reference and Vxc, current band Hamiltonian and live eigenpairs/occupations. Its presence and original path must agree on continuation. |
| Live head/wing context | With active `current_hr`, native DFT H(R), S(R), r(R), lattice and mesh; with active `rotate_reference`, the immutable reference velocity matrices. |

Immutable reference operators are stored for comparison with freshly loaded
inputs. They do not replace the original input dataset. RI coefficients,
Coulomb matrices, frequency-dependent self-energy and screened interaction
are not serialized here; exchange and screening are rebuilt for the next
update. This format has no adaptive-mixer history.

The reader checks the numerical context above, including localized R support,
and compares real/imaginary entries with an absolute `1e-12` tolerance. It
does not check every producer file or numerical option: for example, Padé
controls and arbitrary RI/Coulomb changes are not covered by
`CheckpointContext`. Retain those inputs and controls unchanged as well.

### Completed updates and MPI restoration

The publication points have different meanings:

1. Optional HR export writes the actual mixed/post-cut Hamiltonian before
   diagonalization. An HR snapshot alone does not certify a completed update.
2. After diagonalization and occupation updates succeed, root writes the
   completed state to `<checkpoint>.tmp`, flushes and closes it, then renames
   it to the selected checkpoint filename. A later interrupted update leaves
   the preceding completed checkpoint available.
3. Text iteration/eigenvalue traces and optional band tables are written
   after the checkpoint. They are not part of the checkpoint's atomic rename;
   an output failure can leave a valid checkpoint without all corresponding
   text records.

Checkpoint file operations are root-only. On restore, reference checks,
payload reads and trailing-data checks finish before the temporary restored
state replaces the live state. `collective_root_stage` communicates failure
before any dependent state broadcast. On success, the numerical task broadcasts the
completed iteration number, live mean field and Hamiltonian, including the
band channel when present.

No MPI process partition is serialized. Continuation can change rank count
within the supported route, which requires complete replicated wavefunctions
on every rank. The next update rebuilds its compute objects from the restored
state. Work inside an unfinished self-energy calculation is not resumed.

### Binary format and numerical validation

This is an experimental native-binary format. Its magic string, version,
byte-order marker and velocity tags describe the serialization convention;
they are not content hashes or a portable interchange schema. Incompatible
byte order is rejected.

The writer currently selects versions as follows:

| Active velocity convention | No cut | Shifted cut mode 1 |
| --- | --- | --- |
| No live-velocity context, including legacy `reference` | 3 | 8 |
| `rotate_reference` | 5 | 9 |
| Native-baseline `current_hr` | 10 | 11 |

Legacy versions 1--3 are readable only when their update context still agrees.
Version 4 belongs to another equation/mixer branch; versions 6/7 use the
historical cut convention and are rejected. Versions 5/9 that encode the old
whole-H current-HR map cannot be continued as native-baseline versions 10/11.
Do not reinterpret a prior payload or silently change its physical map when
extending the format. Update both reader and writer, select an explicit
compatible version or reject the old convention, and document the boundary.

Use `regression_tests/verify_qsgw_restart.py` for continuous versus split-run
numerical comparisons. The restored boundary must agree exactly; subsequent
fixed-reference operators and eigenvalues use the existing numerical gates.
Later eigenvector phases or degenerate frames can differ between independent
runs. `regression_tests/verify_qsgw_headwing.py` therefore aligns the entering
eigenframes before comparing velocity and also compares head/wing feedback
and frequency nodes. Both scripts require NumPy, the matching MPI runtime
and a Python version with `tarfile` extraction filters, such as Python 3.12+.

For changes to native reconstruction, retain the finite-metric, native-tail
and SOC numerical controls in `test_rpa_headwing`. These checks establish
implementation consistency; they do not establish material convergence.
