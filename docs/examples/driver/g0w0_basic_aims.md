# One-shot GW with FHI-aims dataset

To run GW calculations, LibRPA needs to be compiled with LibRI switched on.

## Quasi-particle energy levels of water molecule

Folder structure
```
.
├── dataset
│   ├── control.in
│   └── geometry.in
└── librpa
    └── librpa.in
```

`control.in` for FHI-aims dataset
```{literalinclude} ../../../regression_tests/testcases/g0w0_aims_mole_H2O_libri/control.in
:lines: 1-16
```

`geometry.in` for FHI-aims dataset
```{literalinclude} ../../../regression_tests/testcases/g0w0_aims_mole_H2O_libri/geometry.in
```

`librpa.in` for the LibRPA driver
```
# Driver parameters
task = g0w0
input_dir = ../dataset

# API runtime
nfreq = 32
parallel_routing = libri
```

## Quasi-particle band structure of Si

Folder structure
```
.
├── dataset
│   ├── control.in
│   └── geometry.in
└── librpa
    └── librpa.in
```

`control.in` for FHI-aims dataset
```{literalinclude} ../../../regression_tests/testcases/g0w0_band_aims_Si_libri/control.in
:lines: 1-25
```

`geometry.in` for FHI-aims dataset
```{literalinclude} ../../../regression_tests/testcases/g0w0_band_aims_Si_libri/geometry.in
```

`librpa.in` for the LibRPA driver
```
# Driver parameters
task = g0w0_band
input_dir = ../dataset

# API runtime
nfreq = 6
option_dielect_func = 0
replace_w_head = t
parallel_routing = libri
```
