#pragma once

#include <complex>

#ifdef __cplusplus
extern "C"
{
#endif /* __cplusplus */

////////////////////
/* BLACS bindings */
////////////////////
    void blacs_gridinit_(int *ictxt, const char *layout, const int *nprow, const int *npcol);
    void blacs_gridmap_(int *ictxt, const int *usermap, const int *ldup, const int *nprow0, const int *npcol0);
    void blacs_gridinfo_(const int *ictxt, int *nprow, int *npcol, int *myprow, int *mypcol);
    void blacs_setup_(int *pid, int *nprocs);
    void blacs_pinfo_(int *pid, int *nprocs);
    void blacs_get_(const int *ictxt, const int *what, int *val);
    void blacs_set_(const int *ictxt, const int *what, const int *val);
    int numroc_(const int *n, const int *nb, const int *iproc, const int *srcproc, const int *nprocs);
    void descinit_(int *desc,
                   const int *m, const int *n, const int *mb, const int *nb,
                   const int *irsrc, const int *icsrc, const int *ictxt, const int *lld,
                   int *info);
    void blacs_pcoord_(const int *ictxt, const int *pid, int *prow, int *pcol);
    int blacs_pnum_(const int *ictxt, const int *prow, const int *pcol);
    void blacs_gridexit_(const int *ictxt);

    int Csys2blacs_handle(int SysCtxt);
    void Cblacs_pinfo(int *myid, int *nprocs);
    void Cblacs_gridmap(int* icontxt, int *usermap, int ldumap, int nprow, int npcol);
    void Cblacs_gridinfo(int ictxt, int *nprow, int *npcol, int *prow, int *pcol);
    void Cblacs_gridinit(int *ictxt, char *layout, int nprow, int npcol);
    int  Cblacs_pnum(int icontxt, int prow, int pcol);
    void Cblacs_pcoord(int icontxt, int pnum, int *prow, int *pcol);

////////////////////
/* PBLAS bindings */
////////////////////
    void pdpotrf_(char *uplo, int *n, double *a, int *ia, int *ja, int *desca, int *info);
    void pzpotrf_(char *uplo, int *n, std::complex<double> *a, int *ia, int *ja, int *desca, int *info);
    void pdtran_(int *m, int *n, double *alpha,
                 double *a, int *ia, int *ja, int *desca,
                 double *beta, 
                 double *c, int *ic, int *jc, int *descc);

// matrix-vector
    void psgemv_(const char *transa, const int *M, const int *N, const float *alpha,
                 const float *A, const int *IA, const int *JA, const int *DESCA,
                 const float *B, const int *IB, const int *JB, const int *DESCB,
                 const int *K, const float *beta,
                 float *C, const int *IC, const int *JC, const int *DESCC,
                 const int *L);
    void pdgemv_(const char *transa, const int *M, const int *N, const double *alpha,
                 const double *A, const int *IA, const int *JA, const int *DESCA,
                 const double *B, const int *IB, const int *JB, const int *DESCB,
                 const int *K, const double *beta,
                 double *C, const int *IC, const int *JC, const int *DESCC,
                 const int *L);
    void pcgemv_(const char *transa, const int *M, const int *N, const std::complex<float> *alpha,
                 const std::complex<float> *A, const int *IA, const int *JA, const int *DESCA,
                 const std::complex<float> *B, const int *IB, const int *JB, const int *DESCB,
                 const int *K, const std::complex<float> *beta,
                 std::complex<float> *C, const int *IC, const int *JC, const int *DESCC,
                 const int *L);
    void pzgemv_(const char *transa, const int *M, const int *N, const std::complex<double> *alpha,
                 const std::complex<double> *A, const int *IA, const int *JA, const int *DESCA,
                 const std::complex<double> *B, const int *IB, const int *JB, const int *DESCB,
                 const int *K, const std::complex<double> *beta,
                 std::complex<double> *C, const int *IC, const int *JC, const int *DESCC,
                 const int *L);

// matrix-matrix
    void psgemm_(const char *transa, const char *transb,
                 const int *M, const int *N, const int *K,
                 const float *alpha,
                 const float *A, const int *IA, const int *JA, const int *DESCA,
                 const float *B, const int *IB, const int *JB, const int *DESCB,
                 const float *beta,
                 float *C, const int *IC, const int *JC, const int *DESCC);
    void pcgemm_(const char *transa, const char *transb,
                 const int *M, const int *N, const int *K,
                 const std::complex<float> *alpha,
                 const std::complex<float> *A, const int *IA, const int *JA, const int *DESCA,
                 const std::complex<float> *B, const int *IB, const int *JB, const int *DESCB,
                 const std::complex<float> *beta,
                 std::complex<float> *C, const int *IC, const int *JC, const int *DESCC);
    void pdgemm_(const char *transa, const char *transb,
                 const int *M, const int *N, const int *K,
                 const double *alpha,
                 const double *A, const int *IA, const int *JA, const int *DESCA,
                 const double *B, const int *IB, const int *JB, const int *DESCB,
                 const double *beta,
                 double *C, const int *IC, const int *JC, const int *DESCC);
    void pzgemm_(const char *transa, const char *transb,
                 const int *M, const int *N, const int *K,
                 const std::complex<double> *alpha,
                 const std::complex<double> *A, const int *IA, const int *JA, const int *DESCA,
                 const std::complex<double> *B, const int *IB, const int *JB, const int *DESCB,
                 const std::complex<double> *beta,
                 std::complex<double> *C, const int *IC, const int *JC, const int *DESCC);
    void pdsymm_(const char *side, const char *uplo,
                 const int *m, const int *n, const double *alpha,
                 const double *a, const int *ia, const int *ja, const int *desca,
                 const double *b, const int *ib, const int *jb, const int *descb,
                 const double *beta,
                 double *c, const int *ic, const int *jc, const int *descc);
    void pdtrmm_(const char *side, const char *uplo, const char *transa, const char *diag,
                 const int *m, const int *n,
                 const double *alpha,
                 const double *a, const int *ia, const int *ja, const int *desca,
                 double *b, const int *ib, const int *jb, const int *descb);
    void pztrmm_(const char *side, const char *uplo, const char *transa, const char *diag,
                 const int *m, const int *n,
                 const double *alpha, const std::complex<double> *a, const int *ia, const int *ja, const int *desca,
                 std::complex<double> *b, const int *ib, const int *jb, const int *descb);

// factoriazation
    void psgetrf_(const int *M, const int *N, 
                  float *A, const int *IA, const int *JA, const int *DESCA,
                  int *ipiv, int *info);
    void pdgetrf_(const int *M, const int *N, 
                  double *A, const int *IA, const int *JA, const int *DESCA,
                  int *ipiv, int *info);
    void pcgetrf_(const int *M, const int *N, 
                  std::complex<float> *A, const int *IA, const int *JA, const int *DESCA,
                  int *ipiv, int *info);
    void pzgetrf_(const int *M, const int *N, 
                  std::complex<double> *A, const int *IA, const int *JA, const int *DESCA,
                  int *ipiv, int *info);

// Eigenvalues of symmetric and hermitian matrices
    void pssygvx_(const int *itype, const char *jobz, const char *range, const char *uplo,
                  const int *n, float *A, const int *ia, const int *ja, const int*desca,
                  float *B, const int *ib, const int *jb, const int*descb,
                  const float *vl, const float *vu, const int *il, const int *iu,
                  const float *abstol, int *m, int *nz, float *w, const float *orfac,
                  float *Z, const int *iz, const int *jz, const int *descz,
                  float *work, int *lwork, int *iwork, int *liwork, int *ifail, int *iclustr,
                  float *gap, int *info);
    void pdsygvx_(const int *itype, const char *jobz, const char *range, const char *uplo,
                  const int *n, double *A, const int *ia, const int *ja, const int*desca,
                  double *B, const int *ib, const int *jb, const int*descb,
                  const double *vl, const double *vu, const int *il, const int *iu,
                  const double *abstol, int *m, int *nz, double *w, const double *orfac,
                  double *Z, const int *iz, const int *jz, const int *descz,
                  double *work, int *lwork, int *iwork, int *liwork, int *ifail, int *iclustr,
                  double *gap, int *info);
    void pchegvx_(const int *itype, const char *jobz, const char *range, const char *uplo,
                  const int *n, std::complex<float> *A, const int *ia, const int *ja, const int *desca,
                  std::complex<float> *B, const int *ib, const int *jb, const int *descb,
                  const float *vl, const float *vu, const int *il, const int *iu,
                  const float *abstol, int *m, int *nz, float *w, const float *orfac,
                  std::complex<float> *Z, const int *iz, const int *jz, const int*descz,
                  std::complex<float> *work, int *lwork, float *rwork, int *lrwork, int *iwork, int*liwork, int *ifail, int *iclustr,
                  float *gap, int *info);
    void pzhegvx_(const int *itype, const char *jobz, const char *range, const char *uplo,
                  const int *n, std::complex<double> *A, const int *ia, const int *ja, const int *desca,
                  std::complex<double> *B, const int *ib, const int *jb, const int *descb,
                  const double *vl, const double *vu, const int *il, const int *iu,
                  const double *abstol, int *m, int *nz, double *w, const double *orfac,
                  std::complex<double> *Z, const int *iz, const int *jz, const int*descz,
                  std::complex<double> *work, int *lwork, double *rwork, int *lrwork, int *iwork, int*liwork, int *ifail, int *iclustr,
                  double *gap, int *info);

// Matrix inversion
    void psgetri_(const int *n, 
                  float *A, const int *ia, const int *ja, const int *desca,
                  int *ipiv, const float *work, const int *lwork, const int *iwork, const int *liwork,
                  const int *info);
    void pdgetri_(const int *n, 
                  double *A, const int *ia, const int *ja, const int *desca,
                  int *ipiv, const double *work, const int *lwork, const int *iwork, const int *liwork,
                  const int *info);
    void pcgetri_(const int *n, 
                  std::complex<float> *A, const int *ia, const int *ja, const int *desca,
                  int *ipiv, const std::complex<float> *work, const int *lwork, const int *iwork, const int *liwork,
                  const int *info);
    void pzgetri_(const int *n, 
                  std::complex<double> *A, const int *ia, const int *ja, const int *desca,
                  int *ipiv, const std::complex<double> *work, const int *lwork, const int *iwork, const int *liwork,
                  const int *info);

// Matrix addition
    void psgeadd_(const char *transa, const int *m, const int *n,
                  float *alpha,
                  float *a, const int *ia, const int *ja, const int *desca,
                  float *beta,
                  float *c, const int *ic, const int *jc, const int *descc);
    void pdgeadd_(const char *transa, const int *m, const int *n,
                  double *alpha,
                  double *a, const int *ia, const int *ja, const int *desca,
                  double *beta,
                  double *c, const int *ic, const int *jc, const int *descc);
    void pcgeadd_(const char *transa, const int *m, const int *n,
                  const std::complex<float> *alpha,
                  const std::complex<float> *a, const int *ia, const int *ja, const int *desca,
                  const std::complex<float> *beta,
                  std::complex<float> *c, const int *ic, const int *jc, const int *descc);
    void pzgeadd_(const char *transa, const int *m, const int *n,
                  const std::complex<double> *alpha,
                  const std::complex<double> *a, const int *ia, const int *ja, const int *desca,
                  const std::complex<double> *beta,
                  std::complex<double> *c, const int *ic, const int *jc, const int *descc);

// submatrix redistribution
    void psgemr2d_(const int *m, const int *n,
                   const float *a, const int *ia, const int *ja, const int *desca,
                   float *b, const int *ib, const int *jb, const int *descb,
                   const int *ictxt);
    void pdgemr2d_(const int *m, const int *n,
                   const double *a, const int *ia, const int *ja, const int *desca,
                   double *b, const int *ib, const int *jb, const int *descb,
                   const int *ictxt);
    void pcgemr2d_(const int *m, const int *n,
                   const std::complex<float> *a, const int *ia, const int *ja, const int *desca,
                   std::complex<float>*b, const int *ib, const int *jb, const int *descb,
                   const int *ictxt);
    void pzgemr2d_(const int *m, const int *n,
                   const std::complex<double> *a, const int *ia, const int *ja, const int *desca,
                   std::complex<double> *b, const int *ib, const int *jb, const int *descb,
                   const int *ictxt);

#ifdef __cplusplus
}
#endif /* __cplusplus */
