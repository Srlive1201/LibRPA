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
:end-at: [insert species
```

`geometry.in` for FHI-aims dataset
```{literalinclude} ../../../regression_tests/testcases/g0w0_aims_mole_H2O_libri/geometry.in
```

`librpa.in` for the LibRPA driver
```{literalinclude} ../../../regression_tests/testcases/g0w0_aims_mole_H2O_libri/librpa/librpa.in
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
:end-at: [insert species
```

`geometry.in` for FHI-aims dataset
```{literalinclude} ../../../regression_tests/testcases/g0w0_band_aims_Si_libri/geometry.in
```

`librpa.in` for the LibRPA driver
```{literalinclude} ../../../regression_tests/testcases/g0w0_band_aims_Si_libri/librpa/librpa.in
```
