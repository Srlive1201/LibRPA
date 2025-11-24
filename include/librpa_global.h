#pragma once
#include "librpa_enums.h"

#ifdef __cplusplus
extern "C" {
#endif

void librpa_init_global(LibrpaSwitch switch_redirect_stdout, const char *redirect_path);

void librpa_finalize_global();

#ifdef __cplusplus
}
#endif

#ifdef __cplusplus
namespace librpa
{

void init_global(LibrpaSwitch switch_redirect_stdout, const char *redirect_path);

void finalize_global();

}
#endif
