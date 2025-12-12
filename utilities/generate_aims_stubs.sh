#!/usr/bin/env bash

LIBRPA_ROOT="$(dirname "$(dirname "$(realpath "${BASH_SOURCE[0]}")")")"
FILE_STUB="${LIBRPA_ROOT}/binding/fortran/librpa_f03_stubs.f90"
SED=gsed

if [[ ! -f "$FILE_STUB" ]]; then
  echo "Fortran stub file not found: $FILE_STUB"
  exit 2
else
  echo -n "Processing Fortran stub file $FILE_STUB to adapt to FHI-aims ... "
fi

FILE_AIMS_STUB="librpa_f03_stubs_aims.f90"

if "$SED" -e "s/write(\*,'(A)') info_str/call localorb_info(info_str,use_unit,'(A)')/" \
  -e '0,/implicit none/{/implicit none/i\
   use localorb_io, only: localorb_info, use_unit
}' \
  "${FILE_STUB}" >"${FILE_AIMS_STUB}"; then
  echo "done"
else
  echo "fail"
fi
