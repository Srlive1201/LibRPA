# Runtime Parameters

## Overview

Runtime parameters and options in LibRPA are generally managed through the `LibrpaOptions` object.

- When using the **API**, you can directly assign values to the corresponding attributes to control the runtime behavior of LibRPA.

- When using the **driver**, parameters are read from a file named `librpa.in` in the working directory, using the syntax `key = value`.
  For example:
  ```
  task = rpa
  nfreq = 16
  cs_threshold = 1e-4
  input_dir = ../librpa_dataset
  ```

  If `librpa.in` is not found, the driver stops with an error.
  Unlike the API, the driver requires `task` to be explicitly specified.
  For all other parameters, default values are used when the corresponding keywords are not provided.

## Additional Parameters for Driver

| Parameter Name | Description                                              | Type   | Default Value (Options)                         |
|----------------|----------------------------------------------------------|--------|-------------------------------------------------|
| `task`         | Task type                                                | string | (required: rpa, g0w0, exx, g0w0_band, exx_band) |
| `input_dir`    | Input directory to read the AO dataset                   | string | `.`                                             |
| `cs_threshold` | Screening threshold when reading the RI coefficient data | double | 1e-6                                            |

## Common Parameter Settings for LibRPA

| Parameter Name          | Description                                        | Type   | Default Value (Options)                                     |
|-------------------------|----------------------------------------------------|--------|-------------------------------------------------------------|
| `output_dir`            | Output directory for results                       | string | `.`                                                         |
| `parallel_routing`      | Parallel scheme of LibRPA                          | string | auto (auto, atompair, rtau, libri)                          |
| `output_level`          | Verbosity level                                    | int    | 2 (0=silent, 1=critical, 2=info, 3=warn, 4=debug)           |
| `vq_threshold`          | Real-space Coulomb matrices screening threshold    | double | 0.0                                                         |
| `use_kpara_scf_eigvec`  | Flag for parallel distribution of SCF eigenvectors | bool   | `false`                                                     |
| `tfgrids_type`          | Type of time-frequency integration grid            | string | minimax (gl, gci, gcii, minimax, evenspaced, evenspaced_tf) |
| `nfreq`                 | Number of frequency integration grid points        | int    | 6                                                           |
| `tfgrids_freq_min`      | Minimum frequency for grid (Hartree)               | double | 0.0                                                         |
| `tfgrids_freq_interval` | Frequency interval for even-spaced grid (Hartree)  | double | 0.1                                                         |
| `tfgrids_freq_max`      | Maximum frequency for grid (Hartree)               | double | 10.0                                                        |
| `tfgrids_time_min`      | Minimum time for grid (Hartree^-1)                 | double | 0.0                                                         |
| `tfgrids_time_interval` | Time interval for even-spaced grid (Hartree^-1)    | double | 0.1                                                         |

## RPA-specific Parameters

| Parameter Name        | Description                                                           | Type   | Default Value (Options) |
|-----------------------|-----------------------------------------------------------------------|--------|-------------------------|
| `gf_threshold`        | Real-space Green's function screening threshold for response function | double | 1e-4                    |
| `use_scalapack_ecrpa` | Flag to use ScaLapack for calculating $E_\text{c}^{\text{RPA}}$       | bool   | `true`                  |

## GW-specific Parameters

| Parameter Name           | Description                                                                  | Type   | Default Value (Options)             |
|--------------------------|------------------------------------------------------------------------------|--------|-------------------------------------|
| `n_params_anacon`        | Number of parameters for analytic continuation                               | int    | 16                                  |
| `use_scalapack_gw_wc`    | Flag to use ScaLAPACK for computing Wc from chi0                             | bool   | `false`                             |
| `replace_w_head`         | Flag to replace head of dielectric matrix by macroscopic dielectric function | bool   | `false`                             |
| `option_dielect_func`    | Option for computing dielectric function on imaginary axis                   | int    | 0 (0=direct, 1=model fit, 2=spline) |
| `sqrt_coulomb_threshold` | Threshold for eigenvalues to perform square root of Coulomb matrices         | double | 1e-6                                |

## LibRI Parameters

| Parameter Name            | Description                                                           | Type   | Default Value (Options) |
|---------------------------|-----------------------------------------------------------------------|--------|-------------------------|
| `libri_chi0_threshold_C`  | Threshold of LRI triple coefficients for response function            | double | 1e-8                    |
| `libri_chi0_threshold_G`  | Threshold of Green's function for response function                   | double | 1e-8                    |
| `libri_exx_threshold_C`   | Threshold of LRI triple coefficients for exact exchange               | double | 1e-8                    |
| `libri_exx_threshold_D`   | Threshold of density matrices for exact exchange                      | double | 1e-8                    |
| `libri_exx_threshold_V`   | Threshold of Coulomb matrices for exact exchange                      | double | 1e-8                    |
| `libri_g0w0_threshold_C`  | Threshold of LRI triple coefficients for G0W0 correlation self-energy | double | 1e-8                    |
| `libri_g0w0_threshold_G`  | Threshold of Green's function for G0W0 correlation self-energy        | double | 1e-8                    |
| `libri_g0w0_threshold_Wc` | Threshold of screened Coulomb matrix for G0W0 correlation self-energy | double | 1e-8                    |

## Output Control Parameters

| Parameter Name          | Description                                                               | Type | Default Value (Options) |
|-------------------------|---------------------------------------------------------------------------|------|-------------------------|
| `output_gw_sigc_mat`    | Output correlation self-energy matrix (k-space, imaginary frequencies)    | bool | false                   |
| `output_gw_sigc_mat_rt` | Output correlation self-energy matrix (real-space, imaginary time)        | bool | false                   |
| `output_gw_sigc_mat_rf` | Output correlation self-energy matrix (real-space, imaginary frequencies) | bool | false                   |

For details on all parameters, you can visit the API documentation of struct {librpa}`LibrpaOptions`.
