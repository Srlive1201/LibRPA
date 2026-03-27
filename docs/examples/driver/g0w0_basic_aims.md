# One-shot GW with FHI-aims dataset

This page provides input files for example one-shot GW calculations with LibRPA
based on an FHI-aims dataset.

Please refer to the [tutorial page](../../tutorial/g0w0/index.md) for the
corresponding workflow and usage instructions.

To run GW calculations for quasiparticle (QP) properties, LibRPA must be
compiled with LibRI enabled. See the
[installation](../../user_guide/install.md) and
[compile options](../../user_guide/compile_options.md) pages for details.

## QP energy levels of water molecule

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
:end-at: [insert species
```

`geometry.in` for FHI-aims dataset
```{literalinclude} ../../../regression_tests/testcases/g0w0_aims_mole_H2O_libri/geometry.in
```

`librpa.in` for the LibRPA driver
```{literalinclude} ../../../regression_tests/testcases/g0w0_aims_mole_H2O_libri/librpa/librpa.in
```

## QP energy levels of Lithium atom (spin-polarized)

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
```{literalinclude} ../../../regression_tests/testcases/g0w0_aims_Li_atom_libri/control.in
:end-at: [insert species
```

`geometry.in` for FHI-aims dataset
```{literalinclude} ../../../regression_tests/testcases/g0w0_aims_Li_atom_libri/geometry.in
```

`librpa.in` for the LibRPA driver
```{literalinclude} ../../../regression_tests/testcases/g0w0_aims_Li_atom_libri/librpa/librpa.in
```

## QP band structure of Si

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
:end-at: [insert species
```

`geometry.in` for FHI-aims dataset
```{literalinclude} ../../../regression_tests/testcases/g0w0_band_aims_Si_libri/geometry.in
```

`librpa.in` for the LibRPA driver
```{literalinclude} ../../../regression_tests/testcases/g0w0_band_aims_Si_libri/librpa/librpa.in
```
