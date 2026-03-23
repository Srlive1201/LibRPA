#pragma once
#include "librpa_enums.h"

// Global functions that should work no matter the LibRPA environment is set or not,
// except for the test function `librpa_test`.

#ifdef __cplusplus
extern "C" {
#endif

const char* librpa_get_build_info(void);

int librpa_get_major_version(void);

int librpa_get_minor_version(void);

int librpa_get_patch_version(void);

void librpa_init_global(LibrpaSwitch switch_redirect_stdout = LIBRPA_SWITCH_OFF, const char *redirect_path = "stdout",
                        LibrpaSwitch switch_process_output = LIBRPA_SWITCH_ON);

void librpa_finalize_global(void);

void librpa_test(void);

void librpa_print_profile(void);

#ifdef __cplusplus
}
#endif
