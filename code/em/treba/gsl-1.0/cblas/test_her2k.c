#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

#include "tests.h"

void
test_her2k (void) {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 2;
   float alpha[2] = {0.0f, 0.1f};
   float beta = 0.1f;
   float A[] = { 0.531f, 0.721f, -0.848f, 0.826f };
   int lda = 2;
   float B[] = { -0.711f, -0.2f, -0.92f, -0.676f };
   int ldb = 2;
   float C[] = { -0.447f, 0.701f };
   int ldc = 1;
   float C_expected[] = { 0.30322f, 0.0f };
   cblas_cher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cher2k(case 1654) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cher2k(case 1654) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 2;
   float alpha[2] = {0.0f, 0.1f};
   float beta = 0.1f;
   float A[] = { 0.68f, 0.079f, 0.837f, -0.814f };
   int lda = 2;
   float B[] = { -0.986f, 0.024f, 0.584f, -0.248f };
   int ldb = 2;
   float C[] = { 0.477f, -0.551f };
   int ldc = 1;
   float C_expected[] = { 0.120103f, 0.0f };
   cblas_cher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cher2k(case 1655) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cher2k(case 1655) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 111;
   int N = 1;
   int K = 2;
   float alpha[2] = {0.0f, 0.1f};
   float beta = 0.1f;
   float A[] = { 0.354f, -0.63f, -0.85f, 0.426f };
   int lda = 1;
   float B[] = { 0.787f, -0.228f, -0.568f, 0.83f };
   int ldb = 1;
   float C[] = { 0.428f, -0.388f };
   int ldc = 1;
   float C_expected[] = { 0.0331132f, 0.0f };
   cblas_cher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cher2k(case 1656) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cher2k(case 1656) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 111;
   int N = 1;
   int K = 2;
   float alpha[2] = {0.0f, 0.1f};
   float beta = 0.1f;
   float A[] = { -0.49f, 0.224f, -0.606f, 0.46f };
   int lda = 1;
   float B[] = { -0.191f, -0.815f, 0.464f, 0.066f };
   int ldb = 1;
   float C[] = { 0.302f, 0.023f };
   int ldc = 1;
   float C_expected[] = { 0.0679396f, 0.0f };
   cblas_cher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cher2k(case 1657) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cher2k(case 1657) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 113;
   int N = 1;
   int K = 2;
   float alpha[2] = {-1.0f, 0.0f};
   float beta = 0.0f;
   float A[] = { 0.943f, 0.075f, 0.15f, -0.141f };
   int lda = 1;
   float B[] = { -0.962f, 0.422f, -0.592f, -0.789f };
   int ldb = 1;
   float C[] = { 0.728f, 0.601f };
   int ldc = 1;
   float C_expected[] = { 1.70613f, 0.0f };
   cblas_cher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cher2k(case 1658) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cher2k(case 1658) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 113;
   int N = 1;
   int K = 2;
   float alpha[2] = {-1.0f, 0.0f};
   float beta = 0.0f;
   float A[] = { -0.93f, -0.386f, 0.565f, 0.141f };
   int lda = 1;
   float B[] = { -0.801f, 0.022f, 0.558f, -0.932f };
   int ldb = 1;
   float C[] = { 0.068f, 0.501f };
   int ldc = 1;
   float C_expected[] = { -1.84059f, 0.0f };
   cblas_cher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cher2k(case 1659) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cher2k(case 1659) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 113;
   int N = 1;
   int K = 2;
   float alpha[2] = {-1.0f, 0.0f};
   float beta = 0.0f;
   float A[] = { -0.383f, 0.124f, 0.458f, -0.221f };
   int lda = 2;
   float B[] = { -0.107f, 0.199f, 0.18f, 0.122f };
   int ldb = 2;
   float C[] = { 0.896f, -0.874f };
   int ldc = 1;
   float C_expected[] = { -0.24227f, 0.0f };
   cblas_cher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cher2k(case 1660) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cher2k(case 1660) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 113;
   int N = 1;
   int K = 2;
   float alpha[2] = {-1.0f, 0.0f};
   float beta = 0.0f;
   float A[] = { 0.131f, 0.692f, 0.533f, -0.672f };
   int lda = 2;
   float B[] = { -0.435f, -0.453f, 0.195f, -0.579f };
   int ldb = 2;
   float C[] = { -0.547f, 0.736f };
   int ldc = 1;
   float C_expected[] = { -0.245124f, 0.0f };
   cblas_cher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cher2k(case 1661) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cher2k(case 1661) imag");
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
   double beta = 0.1;
   double A[] = { 0.972, -0.353, 0.712, -0.53 };
   int lda = 2;
   double B[] = { 0.787, -0.379, 0.889, 0.901 };
   int ldb = 2;
   double C[] = { 0.002, 0.266 };
   int ldc = 1;
   double C_expected[] = { -0.4278924, 0.0 };
   cblas_zher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zher2k(case 1662) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zher2k(case 1662) imag");
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
   double beta = 0.1;
   double A[] = { -0.36, 0.192, 0.539, 0.198 };
   int lda = 2;
   double B[] = { -0.673, 0.781, 0.792, 0.335 };
   int ldb = 2;
   double C[] = { 0.719, -0.339 };
   int ldc = 1;
   double C_expected[] = { -0.485009, 0.0 };
   cblas_zher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zher2k(case 1663) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zher2k(case 1663) imag");
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
   double beta = 0.1;
   double A[] = { -0.143, 0.456, 0.677, -0.522 };
   int lda = 1;
   double B[] = { 0.851, 0.196, 0.586, 0.64 };
   int ldb = 1;
   double C[] = { 0.617, 0.118 };
   int ldc = 1;
   double C_expected[] = { 0.1081226, 0.0 };
   cblas_zher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zher2k(case 1664) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zher2k(case 1664) imag");
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
   double beta = 0.1;
   double A[] = { 0.801, 0.91, 0.376, -0.006 };
   int lda = 1;
   double B[] = { -0.613, -0.758, -0.966, 0.194 };
   int ldb = 1;
   double C[] = { -0.723, -0.765 };
   int ldc = 1;
   double C_expected[] = { 0.8583678, 0.0 };
   cblas_zher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zher2k(case 1665) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zher2k(case 1665) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 113;
   int N = 1;
   int K = 2;
   double alpha[2] = {-0.3, 0.1};
   double beta = 0.1;
   double A[] = { -0.359, 0.913, 0.738, -0.227 };
   int lda = 1;
   double B[] = { 0.787, 0.745, 0.036, -0.606 };
   int ldb = 1;
   double C[] = { -0.652, -0.281 };
   int ldc = 1;
   double C_expected[] = { -0.1172608, 0.0 };
   cblas_zher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zher2k(case 1666) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zher2k(case 1666) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 113;
   int N = 1;
   int K = 2;
   double alpha[2] = {-0.3, 0.1};
   double beta = 0.1;
   double A[] = { -0.933, 0.598, 0.952, 0.25 };
   int lda = 1;
   double B[] = { -0.508, -0.461, -0.727, 0.162 };
   int ldb = 1;
   double C[] = { 0.215, 0.943 };
   int ldc = 1;
   double C_expected[] = { 0.0795166, 0.0 };
   cblas_zher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zher2k(case 1667) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zher2k(case 1667) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 113;
   int N = 1;
   int K = 2;
   double alpha[2] = {-0.3, 0.1};
   double beta = 0.1;
   double A[] = { -0.735, 0.372, -0.251, -0.168 };
   int lda = 2;
   double B[] = { 0.217, 0.863, -0.179, -0.057 };
   int ldb = 2;
   double C[] = { 0.579, -0.305 };
   int ldc = 1;
   double C_expected[] = { 0.0744312, 0.0 };
   cblas_zher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zher2k(case 1668) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zher2k(case 1668) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 113;
   int N = 1;
   int K = 2;
   double alpha[2] = {-0.3, 0.1};
   double beta = 0.1;
   double A[] = { -0.587, -0.994, -0.625, 0.681 };
   int lda = 2;
   double B[] = { -0.577, -0.014, -0.434, 0.204 };
   int ldb = 2;
   double C[] = { 0.256, 0.093 };
   int ldc = 1;
   double C_expected[] = { -0.3526202, 0.0 };
   cblas_zher2k(order, uplo, trans, N, K, alpha, A, lda, B, ldb, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 1; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zher2k(case 1669) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zher2k(case 1669) imag");
     };
   };
  };


}
