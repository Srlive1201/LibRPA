#pragma once
#include "librpa_enums.h"

#ifdef __cplusplus
extern "C" {
#endif

int librpa_get_major_version(void);

int librpa_get_minor_version(void);

int librpa_get_micro_version(void);

void librpa_init_global(LibrpaSwitch switch_redirect_stdout, const char *redirect_path);

void librpa_finalize_global(void);

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
namespace librpa
{

void init_global(LibrpaSwitch switch_redirect_stdout, const char *redirect_path);

void finalize_global(void);

}
#endif
