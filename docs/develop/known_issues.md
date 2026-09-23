# Known Issues

## Caveats

### QSGW driver access to library internals

The experimental `qsgw` and `qsgw_band` tasks currently access and modify
internal dataset objects owned by a `LibrpaHandler` in `driver/tasks/qsgw.cpp`.
This uses the early-prototype exception to the driver/API separation guideline;
the iteration workflow is not yet exposed through the public C API.
Task-specific input and output remain in `driver/qsgw`, with numerical
operations in `src/core/qsgw`.

### `nan` in GW output

`nan` values in the `sigc` output indicate that the quasi-particle equation (QPE) solver was unable to find a stable solution.
This issue is occasionally encountered for unoccupied states far above the Fermi level (typically beyond 20 eV).

When such cases occur only sporadically and not at consecutive k-points, they usually have little effect on the visual appearance of a plotted band structure, and may be ignored for plotting purposes.
To improve convergence, try a different QPE solver or adjust the solver controls; see runtime options
{ref}`option_qpe_solver <runtime-parameter-option-qpe-solver>`,
{ref}`qpe_solver_thres <runtime-parameter-qpe-solver-thres>`,
{ref}`qpe_solver_n_iter_max <runtime-parameter-qpe-solver-n-iter-max>`,
{ref}`qpe_solver_damp_factor <runtime-parameter-qpe-solver-damp-factor>`, and
{ref}`use_qpe_adaptive_damp <runtime-parameter-use-qpe-adaptive-damp>`.

If `nan` values occur in larger contiguous regions, the GW calculation should be examined more carefully.
This may indicate that the calculation has become unreliable, for example because the filtering threshold is too large or because the analytic continuation produces multiple nearby solutions.

### FHI-aims inputs for high-symmetry rotated cells

We have encountered a case where FHI-aims input generated from a high-symmetry but rotated primitive cell led to a discrepancy between LibRPA calculations with symmetry enabled and disabled.
In that case, the underlying FHI-aims eigenvalue output also differed between the standard and rotated cells.

When investigating EXX or GW symmetry discrepancies with FHI-aims input, be careful with high-symmetry rotated cells and compare against LibRPA results obtained without symmetry.
A relevant issue can be found in [FHI-aims' issue tracker](https://aims-git.rz-berlin.mpg.de/aims/FHIaims/-/work_items/854) (unpublic).

## Planned Improvements

- [ ] Adapt RPA force work by Mohammad in the [backup branch](https://github.com/AESM-Group/LibRPA/tree/master-backup-240416)
