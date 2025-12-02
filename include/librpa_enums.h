#pragma once

// Define enumeration constants and enumerator types used in LibRPA

// C enums
#ifdef __cplusplus
extern "C" {
#endif

#define LIBRPA_UNSET -101
#define LIBRPA_AUTO -51

#define LIBRPA_SWITCH_OFF 0
#define LIBRPA_SWITCH_ON 1

#define LIBRPA_VERBOSE_DEBUG 4
#define LIBRPA_VERBOSE_WARN 3
#define LIBRPA_VERBOSE_INFO 2
#define LIBRPA_VERBOSE_CRITICAL 1
#define LIBRPA_VERBOSE_SILENT 0

// #define LIBRPA_KIND_AIMS 100
// #define LIBRPA_KIND_ABACUS 101
// #define LIBRPA_KIND_OPENMX 102
// #define LIBRPA_KIND_PYSCF 103

#define LIBRPA_ROUTING_COUNT 5
typedef enum
{
    ROUTING_UNSET = LIBRPA_UNSET,
    AUTO = LIBRPA_AUTO,
    RTAU = 0,
    ATOMPAIR = 1,
    LIBRI = 2,
} LibrpaParallelRouting;

#define LIBRPA_TFGRID_COUNT 7
typedef enum
{
    TFGRID_UNSET = LIBRPA_UNSET,
    GaussLegendre = 0,
    GaussChebyshevI = 1,
    GaussChebyshevII = 2,
    Minimax = 3,
    EvenSpaced = 4,
    EvenSpaced_TF = 5,
} LibrpaTimeFreqGrid;

typedef int LibrpaSwitch;

typedef int LibrpaKind;

typedef int LibrpaVerbose;

#ifdef __cplusplus
}
#endif
