# Guidelines in LibRPA developments

## Fortran binding

Files are under `./binding/fortran`. Notes:

- `librpa_f03.f90`: main binding file containing all objects and API definitions.
- `librpa_f03_stubs.f90`:
  Stub file created from `librpa_f03.f90`. Do not edit manually. Instead, update via:

  ```bash
  ../../utilities/convert_fortran_module_to_stub.py librpa_f03.f90 librpa_f03_stubs.f90
  ```

  when under `binding/fortran`.
