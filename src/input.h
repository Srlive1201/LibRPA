#ifndef INPUT_H
#define INPUT_H

/* #include "atoms.h" */
#include <string>
#include <vector>
#include "vector3.h"
#include "vector3_order.h"
#include "matrix3.h"
extern double cs_threshold;
extern double vq_threshold;

extern int kv_nmp[3];
extern Matrix3 latvec;
extern Matrix3 G;
extern std::vector<Vector3_Order<double>> klist;
extern Vector3<double> *kvec_c;

void READ_AIMS_STRU(const std::string &file_path);
#endif
