# Input Parameters



You can create a file named `librpa.in` in your working directory and input the parameter settings as follows:

```
task = rpa
nfreq = 16
gf_threshold = 1e-4
```
If `librpa.in` or the related keyword is not found, the default value will be used.

## Common Parameter Settings for LibRPA

| Parameter Name       | Description                                               | Type   | Default Value (Options)         |
|----------------------|-----------------------------------------------------------|--------|---------------------------------|
| `task`               | Task type                                                 | string | rpa (g0w0)                      |
| `tfgrid_type`        | Type of time-frequency integration grid                   | string | minimax (GL, etc.)   |
| `nfreq`              | Number of frequency integration grid points               | int    | 6                              |
| `gf_threshold`       | Green's function screening threshold                      | double | 1e-3                            |
| `Cs_threshold`       | Auxiliary basis coefficient tensor screening threshold    | double | 1e-4                            |
| `chi_parallel_routing` | Parallel scheme for calculating $\chi^0$                 | string | auto (atompair, rtau, libri)     |
| `use_libri_chi0`     | Whether to use LibRI for calculating $\chi^0$             | bool   | false                           |
| `use_scalapack_ecrpa`| Whether to use ScaLapack for calculating $E_\text{c}^{\text{RPA}}$ | bool   | true                           |
| `debug`              | Whether to enable debug mode                              | bool   | false                           |


For details on all parameters, you can visit {params}`params.h`.