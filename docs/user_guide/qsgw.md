# Experimental fixed-basis QSGW

The `qsgw` and `qsgw_band` driver tasks construct self-energies in the immutable
initial KS basis and update the eigenvalues and wavefunctions each iteration.
The effective Hamiltonian replaces the initial DFT exchange-correlation
potential with exchange and a Hermitian correlation potential. Hartree updates
are not implemented. Only same-grid analytic head correction is supported;
independent-grid head updates and wing updates are rejected.

QSGW-specific input parsing and output live in `driver/qsgw`. Numerical
operations in `src/qsgw` accept matrices and mean-field objects, not filenames.
The existing generic ELSI matrix reader remains in `src/io`.

## Input metadata

Set `qsgw_input_contract` to a file relative to `input_dir`. The contract starts
with `# librpa-qsgw-input-contract-v2`; its table header is `role file`.
Metadata specify producer, internal units, reference basis/gauge, spin/band/AO
counts, grid/path sizes, and enabled updates. The bundled H2O and Si datasets
provide complete examples.

The `vxc_scf_manifest` and (for bands) `vxc_band_manifest` roles point to Vxc
manifests. These start with `# librpa-qsgw-vxc-manifest-v3`, followed by `kind`,
`producer`, `units`, `basis`, and `gauge`. Their table is:

```text
spin k_index kx ky kz rows columns file
```

Spin and k indices are one-based. Paths are relative to the manifest directory.
ABACUS matrices are in Ry and are converted to Ha; NAO matrices are projected
into the initial state basis, while state-basis matrices are already projected.
FHI-aims matrices use the Hartree state-basis ELSI representation.

The new formats omit SHA256 entirely. Convert old contract-v1 / manifest-v2
files by removing the hash column and changing the version header. Old formats
are rejected instead of silently accepting unchecked hashes. Dimensions, units,
k-point matching, finite matrix values and supported execution modes remain
validated by the driver, with errors for unsupported input.

## Numerical checks

The small default cases cover two updates of molecular H2O with linear mixing
and Si k333 with analytic head correction. A manual Si k666 band case is under
`regression_tests/manual/Si_k666_qsgw_band`. These are numerical workflow
regressions, not evidence of converged material predictions.
