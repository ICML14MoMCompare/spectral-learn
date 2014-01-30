#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

#include "tests.h"

void
test_syr2k (void) {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 2;
   float alpha = -0.3f;
   float beta = 0.1f;
   float A[] = { -0.635f, 0.805f };
   int lda = 2;
   float B[] = { 0.773f, 0.375f };
   int ldb = 2;
   float C[] = { 0.616f };
   int ldc = 1;
   float C_expected[] = { 0.174988f };
   cblas_ssyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssyr2k(case 1622)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 2;
   float alpha = -0.3f;
   float beta = 0.1f;
   float A[] = { -0.396f, -0.131f };
   int lda = 2;
   float B[] = { -0.603f, -0.288f };
   int ldb = 2;
   float C[] = { -0.434f };
   int ldc = 1;
   float C_expected[] = { -0.20931f };
   cblas_ssyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssyr2k(case 1623)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 2;
   float alpha = -0.3f;
   float beta = 0.1f;
   float A[] = { -0.085f, -0.444f };
   int lda = 1;
   float B[] = { 0.936f, 0.752f };
   int ldb = 1;
   float C[] = { -0.64f };
   int ldc = 1;
   float C_expected[] = { 0.184069f };
   cblas_ssyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssyr2k(case 1624)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 2;
   float alpha = -0.3f;
   float beta = 0.1f;
   float A[] = { 0.655f, 0.45f };
   int lda = 1;
   float B[] = { 0.16f, -0.747f };
   int ldb = 1;
   float C[] = { 0.576f };
   int ldc = 1;
   float C_expected[] = { 0.19641f };
   cblas_ssyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssyr2k(case 1625)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 112;
   int N = 1;
   int K = 2;
   float alpha = -0.3f;
   float beta = 0.0f;
   float A[] = { 0.259f, -0.334f };
   int lda = 1;
   float B[] = { -0.911f, -0.426f };
   int ldb = 1;
   float C[] = { 0.432f };
   int ldc = 1;
   float C_expected[] = { 0.056199f };
   cblas_ssyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssyr2k(case 1626)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 112;
   int N = 1;
   int K = 2;
   float alpha = -0.3f;
   float beta = 0.0f;
   float A[] = { -0.765f, 0.7f };
   int lda = 1;
   float B[] = { 0.487f, 0.768f };
   int ldb = 1;
   float C[] = { 0.836f };
   int ldc = 1;
   float C_expected[] = { -0.099027f };
   cblas_ssyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssyr2k(case 1627)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 112;
   int N = 1;
   int K = 2;
   float alpha = -0.3f;
   float beta = 0.0f;
   float A[] = { -0.584f, 0.056f };
   int lda = 2;
   float B[] = { 0.928f, -0.101f };
   int ldb = 2;
   float C[] = { -0.529f };
   int ldc = 1;
   float C_expected[] = { 0.328565f };
   cblas_ssyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssyr2k(case 1628)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 112;
   int N = 1;
   int K = 2;
   float alpha = -0.3f;
   float beta = 0.0f;
   float A[] = { 0.25f, 0.8f };
   int lda = 2;
   float B[] = { 0.489f, -0.642f };
   int ldb = 2;
   float C[] = { 0.322f };
   int ldc = 1;
   float C_expected[] = { 0.23481f };
   cblas_ssyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], flteps, "ssyr2k(case 1629)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 2;
   double alpha = 0;
   double beta = 0;
   double A[] = { 0.591, 0.21 };
   int lda = 2;
   double B[] = { -0.718, -0.579 };
   int ldb = 2;
   double C[] = { -0.856 };
   int ldc = 1;
   double C_expected[] = { -0.0 };
   cblas_dsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsyr2k(case 1630)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 2;
   double alpha = 0;
   double beta = 0;
   double A[] = { -0.971, -0.824 };
   int lda = 2;
   double B[] = { -0.227, 0.457 };
   int ldb = 2;
   double C[] = { 0.521 };
   int ldc = 1;
   double C_expected[] = { 0.0 };
   cblas_dsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsyr2k(case 1631)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 2;
   double alpha = 0;
   double beta = 0;
   double A[] = { -0.274, 0.583 };
   int lda = 1;
   double B[] = { 0.668, -0.83 };
   int ldb = 1;
   double C[] = { 0.907 };
   int ldc = 1;
   double C_expected[] = { 0.0 };
   cblas_dsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsyr2k(case 1632)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 2;
   double alpha = 0;
   double beta = 0;
   double A[] = { -0.512, -0.436 };
   int lda = 1;
   double B[] = { -0.443, -0.259 };
   int ldb = 1;
   double C[] = { -0.667 };
   int ldc = 1;
   double C_expected[] = { 0.0 };
   cblas_dsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsyr2k(case 1633)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 112;
   int N = 1;
   int K = 2;
   double alpha = -0.3;
   double beta = 0;
   double A[] = { 0.741, -0.341 };
   int lda = 1;
   double B[] = { 0.743, -0.315 };
   int ldb = 1;
   double C[] = { -0.776 };
   int ldc = 1;
   double C_expected[] = { -0.3947868 };
   cblas_dsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsyr2k(case 1634)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 112;
   int N = 1;
   int K = 2;
   double alpha = -0.3;
   double beta = 0;
   double A[] = { 0.03, 0.175 };
   int lda = 1;
   double B[] = { -0.832, 0.291 };
   int ldb = 1;
   double C[] = { 0.281 };
   int ldc = 1;
   double C_expected[] = { -0.015579 };
   cblas_dsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsyr2k(case 1635)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 112;
   int N = 1;
   int K = 2;
   double alpha = -0.3;
   double beta = 0;
   double A[] = { 0.476, 0.496 };
   int lda = 2;
   double B[] = { -0.626, -0.159 };
   int ldb = 2;
   double C[] = { -0.964 };
   int ldc = 1;
   double C_expected[] = { 0.226104 };
   cblas_dsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsyr2k(case 1636)");
     }
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 112;
   int N = 1;
   int K = 2;
   double alpha = -0.3;
   double beta = 0;
   double A[] = { -0.489, 0.611 };
   int lda = 2;
   double B[] = { -0.285, -0.673 };
   int ldb = 2;
   double C[] = { -0.11 };
   int ldc = 1;
   double C_expected[] = { 0.1631028 };
   cblas_dsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[i], C_expected[i], dbleps, "dsyr2k(case 1637)");
     }
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 2;
   float alpha[2] = {1.0f, 0.0f};
   float beta[2] = {0.0f, 0.1f};
   float A[] = { 0.796f, 0.872f, -0.919f, 0.748f };
   int lda = 2;
   float B[] = { -0.945f, 0.915f, -0.252f, -0.276f };
   int ldb = 2;
   float C[] = { 0.07f, -0.957f };
   int ldc = 1;
   float C_expected[] = { -2.12843f, -0.054104f };
   cblas_csyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyr2k(case 1638) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyr2k(case 1638) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 2;
   float alpha[2] = {1.0f, 0.0f};
   float beta[2] = {0.0f, 0.1f};
   float A[] = { 0.984f, 0.526f, 0.284f, 0.806f };
   int lda = 2;
   float B[] = { -0.509f, -0.178f, 0.188f, -0.221f };
   int ldb = 2;
   float C[] = { -0.388f, 0.795f };
   int ldc = 1;
   float C_expected[] = { -0.43092f, -0.747044f };
   cblas_csyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyr2k(case 1639) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyr2k(case 1639) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 2;
   float alpha[2] = {1.0f, 0.0f};
   float beta[2] = {0.0f, 0.1f};
   float A[] = { -0.16f, 0.628f, -0.06f, -0.645f };
   int lda = 1;
   float B[] = { 0.846f, 0.545f, 0.032f, 0.493f };
   int ldb = 1;
   float C[] = { -0.041f, -0.621f };
   int ldc = 1;
   float C_expected[] = { -0.26101f, 0.783636f };
   cblas_csyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyr2k(case 1640) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyr2k(case 1640) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 2;
   float alpha[2] = {1.0f, 0.0f};
   float beta[2] = {0.0f, 0.1f};
   float A[] = { -0.478f, -0.556f, 0.519f, 0.177f };
   int lda = 1;
   float B[] = { -0.946f, 0.423f, -0.859f, 0.736f };
   int ldb = 1;
   float C[] = { -0.54f, -0.035f };
   int ldc = 1;
   float C_expected[] = { 0.226066f, 1.05345f };
   cblas_csyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyr2k(case 1641) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyr2k(case 1641) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 112;
   int N = 1;
   int K = 2;
   float alpha[2] = {0.0f, 0.0f};
   float beta[2] = {1.0f, 0.0f};
   float A[] = { -0.582f, 0.09f, -0.176f, 0.784f };
   int lda = 1;
   float B[] = { 0.687f, -0.859f, 0.945f, 0.756f };
   int ldb = 1;
   float C[] = { -0.663f, -0.186f };
   int ldc = 1;
   float C_expected[] = { -0.663f, -0.186f };
   cblas_csyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyr2k(case 1642) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyr2k(case 1642) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 112;
   int N = 1;
   int K = 2;
   float alpha[2] = {0.0f, 0.0f};
   float beta[2] = {1.0f, 0.0f};
   float A[] = { 0.231f, -0.452f, -0.112f, -0.837f };
   int lda = 1;
   float B[] = { -0.258f, 0.464f, -0.224f, 0.893f };
   int ldb = 1;
   float C[] = { -0.448f, 0.046f };
   int ldc = 1;
   float C_expected[] = { -0.448f, 0.046f };
   cblas_csyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyr2k(case 1643) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyr2k(case 1643) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 112;
   int N = 1;
   int K = 2;
   float alpha[2] = {0.0f, 0.0f};
   float beta[2] = {1.0f, 0.0f};
   float A[] = { 0.115f, 0.178f, -0.193f, -0.491f };
   int lda = 2;
   float B[] = { 0.545f, -0.665f, 0.979f, -0.4f };
   int ldb = 2;
   float C[] = { 0.522f, 0.712f };
   int ldc = 1;
   float C_expected[] = { 0.522f, 0.712f };
   cblas_csyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyr2k(case 1644) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyr2k(case 1644) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 112;
   int N = 1;
   int K = 2;
   float alpha[2] = {0.0f, 0.0f};
   float beta[2] = {1.0f, 0.0f};
   float A[] = { -0.725f, -0.808f, -0.244f, 0.145f };
   int lda = 2;
   float B[] = { 0.447f, -0.413f, -0.226f, -0.585f };
   int ldb = 2;
   float C[] = { -0.531f, 0.227f };
   int ldc = 1;
   float C_expected[] = { -0.531f, 0.227f };
   cblas_csyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "csyr2k(case 1645) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "csyr2k(case 1645) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 2;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {-0.3, 0.1};
   double A[] = { 0.337, -0.737, -0.993, 0.69 };
   int lda = 2;
   double B[] = { -0.39, -0.836, -0.32, 0.368 };
   int ldb = 2;
   double C[] = { 0.844, -0.763 };
   int ldc = 1;
   double C_expected[] = { 0.3494384, 0.5248712 };
   cblas_zsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyr2k(case 1646) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyr2k(case 1646) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 2;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {-0.3, 0.1};
   double A[] = { 0.386, -0.465, 0.719, -0.378 };
   int lda = 2;
   double B[] = { 0.099, -0.879, 0.864, 0.141 };
   int ldb = 2;
   double C[] = { -0.599, -0.47 };
   int ldc = 1;
   double C_expected[] = { 0.1664126, 0.5082238 };
   cblas_zsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyr2k(case 1647) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyr2k(case 1647) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 2;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {-0.3, 0.1};
   double A[] = { 0.914, 0.128, -0.262, -0.26 };
   int lda = 1;
   double B[] = { 0.431, 0.276, 0.75, 0.904 };
   int ldb = 1;
   double C[] = { 0.287, 0.537 };
   int ldc = 1;
   double C_expected[] = { -0.3532044, 0.0216788 };
   cblas_zsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyr2k(case 1648) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyr2k(case 1648) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 2;
   double alpha[2] = {-0.3, 0.1};
   double beta[2] = {-0.3, 0.1};
   double A[] = { 0.618, 0.72, 0.392, -0.737 };
   int lda = 1;
   double B[] = { 0.783, 0.531, 0.375, 0.203 };
   int ldb = 1;
   double C[] = { 0.058, -0.116 };
   int ldc = 1;
   double C_expected[] = { -0.3837348, -0.2968344 };
   cblas_zsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyr2k(case 1649) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyr2k(case 1649) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 112;
   int N = 1;
   int K = 2;
   double alpha[2] = {0, 1};
   double beta[2] = {1, 0};
   double A[] = { -0.372, -0.735, -0.711, 0.051 };
   int lda = 1;
   double B[] = { 0.257, 0.097, 0.338, -0.484 };
   int ldb = 1;
   double C[] = { -0.142, -0.197 };
   int ldc = 1;
   double C_expected[] = { -0.414766, -0.676886 };
   cblas_zsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyr2k(case 1650) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyr2k(case 1650) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 112;
   int N = 1;
   int K = 2;
   double alpha[2] = {0, 1};
   double beta[2] = {1, 0};
   double A[] = { 0.1, -0.878, 0.28, -0.381 };
   int lda = 1;
   double B[] = { -0.208, 0.309, -0.276, 0.123 };
   int ldb = 1;
   double C[] = { 0.483, -0.541 };
   int ldc = 1;
   double C_expected[] = { -0.22324, -0.10083 };
   cblas_zsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyr2k(case 1651) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyr2k(case 1651) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 112;
   int N = 1;
   int K = 2;
   double alpha[2] = {0, 1};
   double beta[2] = {1, 0};
   double A[] = { -0.918, 0.515, -0.985, 0.067 };
   int lda = 2;
   double B[] = { -0.034, 0.408, 0.66, -0.945 };
   int ldb = 2;
   double C[] = { -0.063, -0.018 };
   int ldc = 1;
   double C_expected[] = { -1.228982, -1.549386 };
   cblas_zsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyr2k(case 1652) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyr2k(case 1652) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 112;
   int N = 1;
   int K = 2;
   double alpha[2] = {0, 1};
   double beta[2] = {1, 0};
   double A[] = { 0.443, -0.009, -0.245, -0.008 };
   int lda = 2;
   double B[] = { 0.495, -0.239, -0.973, -0.032 };
   int ldb = 2;
   double C[] = { -0.85, -0.799 };
   int ldc = 1;
   double C_expected[] = { -0.660584, 0.111526 };
   cblas_zsyr2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zsyr2k(case 1653) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zsyr2k(case 1653) imag");
     };
   };
  };


}
