# Input Parameters

For Driver Usage, you can create a file named `librpa.in` in your working directory and input the parameter settings as follows:

```
task = rpa
nfreq = 16
gf_R_threshold = 1e-4
```
If `librpa.in` or the related keyword is not found, the default value will be used.

## Common Parameter Settings for LibRPA

| Parameter Name         | Description                                                           | Type   | Default Value (Options)      |
|------------------------|-----------------------------------------------------------------------|--------|------------------------------|
| `task`                 | Task type                                                             | string | rpa (rpa, g0w0, exx)         |
| `tfgrid_type`          | Type of time-frequency integration grid                               | string | minimax (minimax)            |
| `nfreq`                | Number of frequency integration grid points                           | int    | 6                            |
| `gf_R_threshold`       | Real-space Green's function screening threshold for response function | double | 1e-3                         |
| `cs_threshold`         | Auxiliary basis coefficient tensor screening threshold                | double | 1e-4                         |
| `parallel_routing`     | Parallel scheme of LibRPA                                             | string | auto (atompair, rtau, libri) |
| `use_scalapack_ecrpa`  | Flag to use ScaLapack for calculating $E_\text{c}^{\text{RPA}}$       | bool   | true                         |
| `debug`                | Flag to enable debug mode                                             | bool   | false                        |

For details on all parameters, you can visit the API documentation of struct {librpa}`LibRPAParams`.
