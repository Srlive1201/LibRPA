# Runtime Parameters

## Overview

Runtime parameters and options in LibRPA are generally managed through the `LibrpaOptions` object.

- When using the **API**, you can directly assign values to the corresponding attributes to control the runtime behavior of LibRPA.
  API defaults are set by {librpa}`librpa_init_options`, the C++ {librpa}`librpa::Options` constructor, or `opts%init()` in the Fortran API.

- When using the **driver**, parameters are read from a file named `librpa.in` in the working directory, using the syntax `key = value`.
  For example:
  ```
  task = rpa
  input_dir = ../librpa_dataset
  nfreq = 16
  cs_threshold = 1e-4
  ```

  If `librpa.in` is not found, the driver stops with an error.
  Unlike the API, the driver requires `task` to be explicitly specified.
  For all other parameters, default values are used when the corresponding keywords are not provided.
  Please refer to [manual page of driver usage](driver_usage) for more information.

```{note}
Rows marked `Experimental` expose development, diagnostic, or less-tested code paths.
Use them as opt-in features and validate results against reference calculations before production use.
Blank status entries indicate the regular, documented path.
```

For API use, switch values are set with `LIBRPA_SWITCH_ON` and `LIBRPA_SWITCH_OFF` in C/C++.
The Fortran wrapper exposes these fields as logical values.

(additional-parameters-for-driver)=
## Additional Parameters for Driver

These parameters are available only through the standalone driver.

| Parameter Name        | Description                                                          | Type   | Default Value (Options)                         | Status       |
|-----------------------|----------------------------------------------------------------------|--------|-------------------------------------------------|--------------|
| `task`                | Task type                                                            | string | required (`rpa`, `g0w0`, `exx`, `g0w0_band`, `exx_band`) |              |
| `constants_choice`    | Physical constants source                                            | string | `internal` (`internal`, `aims`)                |              |
| `input_dir`           | Input directory to find and read the AO dataset                      | string | `./`                                           |              |
| `cs_threshold`        | Screening threshold when reading the RI coefficient data             | double | 1e-6                                           |              |
| `use_spinor_wfc`      | Read wavefunctions in spinor format                                  | bool   | `false`                                        | Experimental |
| `output_energy_qp`    | Output quasiparticle energies for external BSE workflows             | bool   | `false`                                        | Experimental |
| `output_hamgnn`       | Output GW energies for HamGNN machine-learning workflows             | bool   | `false`                                        | Experimental |
| `use_pyatb`           | Use PyATB mean-field data for dielectric head/wing calculations      | bool   | `false`                                        | Experimental |
| `output_gw_spec_func` | Output GW spectral-function data                                     | bool   | `false`                                        | Experimental |
| `sf_omega_start`      | Starting frequency for spectral-function output                      | double | 0.0                                            | Experimental |
| `sf_omega_end`        | Ending frequency for spectral-function output                        | double | 1.0                                            | Experimental |
| `sf_omega_step`       | Frequency step for spectral-function output                          | double | 0.1                                            | Experimental |
| `sf_gf_omega_shift`   | Broadening/shift used for Green's function in spectral output        | double | 0.01                                           | Experimental |
| `sf_sigc_omega_shift` | Broadening/shift used for correlation self-energy in spectral output | double | 0.01                                           | Experimental |
| `sf_state_start`      | First state index for spectral-function output                       | int    | 0                                              | Experimental |
| `sf_state_end`        | Last state index for spectral-function output                        | int    | 10000                                          | Experimental |

(common-parameter-settings-for-librpa)=
## Common Parameter Settings for LibRPA

The default values below are API defaults unless a driver-specific default is listed.

| Parameter Name          | Description                                            | Type   | Default Value (Options)                                      | Status       |
|-------------------------|--------------------------------------------------------|--------|--------------------------------------------------------------|--------------|
| `output_dir`            | Output directory for results                           | string | API: `.`; driver: `librpa.d/`                                |              |
| `parallel_routing`      | Parallel scheme of LibRPA                              | enum/string | API: `AUTO`; driver: `auto` (`auto`, `atompair`, `rtau`, `libri`) |              |
| `output_level`          | Verbosity level                                        | int/string | API: `LIBRPA_VERBOSE_INFO`; driver: `info` (`silent`, `critical`, `info`, `warn`, `debug`) |              |
| `vq_threshold`          | Real-space Coulomb matrix screening threshold          | double | 0.0                                                          |              |
| `use_kpara_scf_eigvec`  | Use k-point-parallel distribution of SCF eigenvectors  | bool   | `false`                                                      | Experimental |
| `tfgrids_type`          | Type of time-frequency integration grid                | enum/string | API: `TFGRID_UNSET`; driver: `minimax` (`GL`, `GC-I`, `GL-II`, `minimax`, `evenspaced`, `evenspaced_tf`) |              |
| `nfreq`                 | Number of frequency integration grid points            | int    | 6                                                            |              |
| `tfgrids_freq_min`      | Minimum frequency for grid (Hartree)                   | double | 0.005                                                        |              |
| `tfgrids_freq_interval` | Frequency interval for even-spaced grid (Hartree)      | double | 0.0                                                          |              |
| `tfgrids_freq_max`      | Maximum frequency for grid (Hartree)                   | double | 1000.0                                                       |              |
| `tfgrids_time_min`      | Minimum time for grid (Hartree^-1)                     | double | 0.005                                                        |              |
| `tfgrids_time_interval` | Time interval for even-spaced grid (Hartree^-1)        | double | 0.0                                                          |              |
| `minimax_emin`          | Minimum transition energy for minimax grid generation  | double | -1.0 (< 0 for automatic setup)                               | Experimental |
| `minimax_emax`          | Maximum transition energy for minimax grid generation  | double | -1.0 (< 0 for automatic setup)                               | Experimental |
| `minimax_regulation`    | Regulation parameter for minimax transfomration matrix | double | 0.0                                                          | Experimental |

When using the API, set `tfgrids_type` explicitly before computation.
The driver maps the unset value to `minimax` for backward compatibility.

## Coulomb and Band Summation Controls

| Parameter Name      | Description                                                      | Type | Default Value (Options) | Status       |
|---------------------|------------------------------------------------------------------|------|-------------------------|--------------|
| `use_fullcoul_eps`  | Use full Coulomb interaction in $\varepsilon = 1 - v \chi^0$     | bool | `true`                  | Experimental |
| `use_fullcoul_exx`  | Use full Coulomb interaction in the exact-exchange operator      | bool | `false`                 | Experimental |
| `use_fullcoul_wc`   | Use full Coulomb interaction in $W^c = (\varepsilon^{-1} - 1) v$ | bool | `false`                 | Experimental |
| `n_bands_chi0`      | Maximum number of bands for response-function construction       | int  | -1 (< 0 for all bands)  | Experimental |
| `n_bands_sigc`      | Maximum number of bands for correlation self-energy construction | int  | -1 (< 0 for all bands)  | Experimental |

## RPA-specific Parameters

| Parameter Name        | Description                                                           | Type   | Default Value (Options) | Status |
|-----------------------|-----------------------------------------------------------------------|--------|-------------------------|--------|
| `gf_threshold`        | Real-space Green's function screening threshold for response function | double | 0.0                     |        |
| `use_scalapack_ecrpa` | Use ScaLAPACK to calculate $E_\text{c}^{\text{RPA}}$                  | bool   | `true`                  |        |

## ABF Compression Parameters

| Parameter Name    | Description                                                                  | Type | Default Value (Options) | Status       |
|-------------------|------------------------------------------------------------------------------|------|-------------------------|--------------|
| `use_shrink_abfs` | Use a compressed auxiliary basis                                             | bool | `false`                 | Experimental |
| `use_shrink_chi`  | Build response matrices directly in the compressed auxiliary basis           | bool | `false`                 | Experimental |

## GW-specific Parameters

| Parameter Name           | Description                                                                  | Type   | Default Value (Options)               | Status       |
|--------------------------|------------------------------------------------------------------------------|--------|---------------------------------------|--------------|
| `n_params_anacon`        | Number of parameters for analytic continuation                               | int    | -1 (uses all `nfreq` data)            |              |
| `use_scalapack_gw_wc`    | Use ScaLAPACK for computing $W^c$ from $\chi^0$                              | bool   | `true`                                |              |
| `use_cholesky_gw_wc`     | Use Cholesky factorization for computing $W^c$ from $\chi^0$                 | bool   | `false`                               | Experimental |
| `replace_w_head`         | Replace the dielectric matrix head by the macroscopic dielectric function    | bool   | `false`                               | Experimental |
| `option_dielect_func`    | Option for dielectric function on the imaginary axis                         | int    | 0 (0=direct, 1=spline, 2=model fit)   | Experimental |
| `use_2d_dielectric`      | Use the 2D dielectric-function branch where supported                        | bool   | `false`                               | Experimental |
| `load_sigc_from_file`    | Load correlation self-energy matrix from file where supported                | bool   | `false`                               | Experimental |
| `sqrt_coulomb_threshold` | Threshold for eigenvalues when taking the square root of Coulomb matrices    | double | 0.0                                   |              |

## LibRI Parameters

These parameters are active for LibRI-enabled builds and LibRI routing.

| Parameter Name            | Description                                                           | Type   | Default Value (Options) | Status |
|---------------------------|-----------------------------------------------------------------------|--------|-------------------------|--------|
| `libri_chi0_threshold_C`  | Threshold of LRI triple coefficients for response function            | double | 0.0                     |        |
| `libri_chi0_threshold_G`  | Threshold of Green's function for response function                   | double | 0.0                     |        |
| `libri_exx_threshold_C`   | Threshold of LRI triple coefficients for exact exchange               | double | 0.0                     |        |
| `libri_exx_threshold_D`   | Threshold of density matrices for exact exchange                      | double | 0.0                     |        |
| `libri_exx_threshold_V`   | Threshold of Coulomb matrices for exact exchange                      | double | 0.0                     |        |
| `libri_g0w0_threshold_C`  | Threshold of LRI triple coefficients for G0W0 correlation self-energy | double | 0.0                     |        |
| `libri_g0w0_threshold_G`  | Threshold of Green's function for G0W0 correlation self-energy        | double | 0.0                     |        |
| `libri_g0w0_threshold_Wc` | Threshold of screened Coulomb matrix for G0W0 correlation self-energy | double | 0.0                     |        |

## Output Control Parameters

| Parameter Name            | Description                                                               | Type | Default Value (Options)                   | Status       |
|---------------------------|---------------------------------------------------------------------------|------|-------------------------------------------|--------------|
| `output_gw_sigc_mat`      | Output correlation self-energy matrix (k-space, imaginary frequencies)    | bool | `false`                                   | Experimental |
| `output_gw_sigc_mat_rt`   | Output correlation self-energy matrix (real-space, imaginary time)        | bool | `false`                                   | Experimental |
| `output_gw_sigc_mat_rf`   | Output correlation self-energy matrix (real-space, imaginary frequencies) | bool | `false`                                   | Experimental |
| `option_output_Wc_Rf_mat` | Output $W^c$ matrix in real space and imaginary frequency where supported | int  | 0 (0=off, 1=lowest frequency, 2=all frequencies) | Experimental |

## Driver Compatibility Aliases

The driver still accepts a few older input keys for backward compatibility,
and may be removed in the future:

| Alias Name       | Preferred Name       | Notes                                      |
|------------------|----------------------|--------------------------------------------|
| `debug`          | `output_level`       | `debug = true` maps to `output_level = debug` |
| `use_soc`        | `use_spinor_wfc`     | Old spelling for selecting spinor-format driver input |
| `tfgrid_type`    | `tfgrids_type`       | Old spelling                               |
| `gf_R_threshold` | `gf_threshold`       | Old spelling                               |

For details on all parameters, you can visit the API documentation of struct {librpa}`LibrpaOptions`.
