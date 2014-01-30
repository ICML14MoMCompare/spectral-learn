#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

#include "tests.h"

void
test_herk (void) {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int uplo = 121;
   int trans = 111;
   int N = 2;
   int K = 1;
   float alpha = 1.0f;
   float beta = 1.0f;
   float A[] = { 0.934f, 0.664f, 0.426f, 0.263f };
   int lda = 1;
   float C[] = { 0.251f, -0.97f, 0.76f, -0.349f, 0.152f, -0.899f, -0.17f, 0.707f };
   int ldc = 2;
   float C_expected[] = { 1.56425f, 0.0f, 1.33252f, -0.311778f, 0.152f, -0.899f, 0.080645f, 0.0f };
   cblas_cherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cherk(case 1606) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cherk(case 1606) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 111;
   int N = 2;
   int K = 1;
   float alpha = 1.0f;
   float beta = 1.0f;
   float A[] = { 0.16f, 0.464f, -0.623f, 0.776f };
   int lda = 2;
   float C[] = { 0.771f, -0.449f, 0.776f, 0.112f, -0.134f, 0.317f, 0.547f, -0.551f };
   int ldc = 2;
   float C_expected[] = { 1.0119f, 0.0f, 0.776f, 0.112f, 0.126384f, -0.096232f, 1.5373f, 0.0f };
   cblas_cherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cherk(case 1607) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cherk(case 1607) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 113;
   int N = 2;
   int K = 1;
   float alpha = 0.1f;
   float beta = 1.0f;
   float A[] = { 0.787f, 0.057f, -0.49f, 0.47f };
   int lda = 2;
   float C[] = { -0.758f, 0.912f, 0.992f, -0.356f, 0.584f, 0.806f, 0.965f, 0.674f };
   int ldc = 2;
   float C_expected[] = { -0.695738f, 0.0f, 0.956116f, -0.316218f, 0.584f, 0.806f, 1.0111f, 0.0f };
   cblas_cherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cherk(case 1608) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cherk(case 1608) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 113;
   int N = 2;
   int K = 1;
   float alpha = 0.1f;
   float beta = 1.0f;
   float A[] = { 0.961f, -0.384f, 0.165f, 0.395f };
   int lda = 1;
   float C[] = { -0.186f, 0.404f, -0.873f, 0.09f, -0.451f, -0.972f, -0.203f, -0.304f };
   int ldc = 2;
   float C_expected[] = { -0.0789023f, 0.0f, -0.873f, 0.09f, -0.450312f, -0.927704f, -0.184675f, 0.0f };
   cblas_cherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cherk(case 1609) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cherk(case 1609) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 111;
   int N = 2;
   int K = 1;
   float alpha = 0.0f;
   float beta = -0.3f;
   float A[] = { 0.04f, 0.608f, 0.21f, -0.44f };
   int lda = 1;
   float C[] = { 0.285f, -0.943f, 0.581f, -0.56f, 0.112f, 0.529f, 0.16f, -0.913f };
   int ldc = 2;
   float C_expected[] = { -0.0855f, 0.0f, 0.581f, -0.56f, -0.0336f, -0.1587f, -0.048f, 0.0f };
   cblas_cherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cherk(case 1610) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cherk(case 1610) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 111;
   int N = 2;
   int K = 1;
   float alpha = 0.0f;
   float beta = -0.3f;
   float A[] = { -0.984f, -0.398f, -0.379f, 0.919f };
   int lda = 2;
   float C[] = { -0.44f, -0.087f, 0.156f, -0.945f, -0.943f, -0.355f, 0.577f, 0.053f };
   int ldc = 2;
   float C_expected[] = { 0.132f, 0.0f, -0.0468f, 0.2835f, -0.943f, -0.355f, -0.1731f, 0.0f };
   cblas_cherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cherk(case 1611) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cherk(case 1611) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 113;
   int N = 2;
   int K = 1;
   float alpha = 1.0f;
   float beta = -1.0f;
   float A[] = { 0.269f, -0.428f, -0.029f, 0.964f };
   int lda = 2;
   float C[] = { 0.473f, -0.932f, -0.689f, -0.072f, -0.952f, -0.862f, 0.001f, 0.282f };
   int ldc = 2;
   float C_expected[] = { -0.217455f, 0.0f, -0.689f, -0.072f, 0.531607f, 0.615096f, 0.929137f, 0.0f };
   cblas_cherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cherk(case 1612) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cherk(case 1612) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 113;
   int N = 2;
   int K = 1;
   float alpha = 1.0f;
   float beta = -1.0f;
   float A[] = { -0.303f, -0.037f, -0.411f, -0.243f };
   int lda = 1;
   float C[] = { 0.652f, -0.227f, -0.849f, 0.87f, -0.051f, -0.535f, 0.418f, -0.681f };
   int ldc = 2;
   float C_expected[] = { -0.558822f, 0.0f, 0.982524f, -0.928422f, -0.051f, -0.535f, -0.19003f, 0.0f };
   cblas_cherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], flteps, "cherk(case 1613) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], flteps, "cherk(case 1613) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 111;
   int N = 2;
   int K = 1;
   double alpha = 0.1;
   double beta = 1;
   double A[] = { -0.384, -0.851, 0.518, 0.492 };
   int lda = 1;
   double C[] = { -0.117, -0.194, -0.915, 0.069, 0.445, 0.089, 0.213, -0.889 };
   int ldc = 2;
   double C_expected[] = { -0.0298343, 0.0, -0.9767604, 0.043811, 0.445, 0.089, 0.2640388, 0.0 };
   cblas_zherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zherk(case 1614) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zherk(case 1614) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 111;
   int N = 2;
   int K = 1;
   double alpha = 0.1;
   double beta = 1;
   double A[] = { 0.13, 0.236, 0.788, 0.629 };
   int lda = 2;
   double C[] = { 0.021, -0.376, -0.804, 0.689, -0.912, 0.21, -0.581, 0.406 };
   int ldc = 2;
   double C_expected[] = { 0.0282596, 0.0, -0.804, 0.689, -0.8869116, 0.2204198, -0.4793415, 0.0 };
   cblas_zherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zherk(case 1615) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zherk(case 1615) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 121;
   int trans = 113;
   int N = 2;
   int K = 1;
   double alpha = 0;
   double beta = 1;
   double A[] = { 0.593, 0.846, -0.144, 0.128 };
   int lda = 2;
   double C[] = { 0.02, 0.313, 0.222, 0.301, 0.412, -0.645, -0.411, -0.02 };
   int ldc = 2;
   double C_expected[] = { 0.02, 0.313, 0.222, 0.301, 0.412, -0.645, -0.411, -0.02 };
   cblas_zherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zherk(case 1616) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zherk(case 1616) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 121;
   int trans = 113;
   int N = 2;
   int K = 1;
   double alpha = 0;
   double beta = 1;
   double A[] = { 0.857, 0.994, -0.933, 0.069 };
   int lda = 1;
   double C[] = { 0.253, -0.521, 0.937, -0.73, 0.24, 0.177, -0.27, -0.225 };
   int ldc = 2;
   double C_expected[] = { 0.253, -0.521, 0.937, -0.73, 0.24, 0.177, -0.27, -0.225 };
   cblas_zherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zherk(case 1617) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zherk(case 1617) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 111;
   int N = 2;
   int K = 1;
   double alpha = 0;
   double beta = 1;
   double A[] = { -0.343, -0.433, -0.381, -0.087 };
   int lda = 1;
   double C[] = { -0.695, 0.911, 0.719, -0.074, -0.621, -0.256, 0.216, -0.889 };
   int ldc = 2;
   double C_expected[] = { -0.695, 0.911, 0.719, -0.074, -0.621, -0.256, 0.216, -0.889 };
   cblas_zherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zherk(case 1618) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zherk(case 1618) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 111;
   int N = 2;
   int K = 1;
   double alpha = 0;
   double beta = 1;
   double A[] = { -0.887, 0.557, -0.43, 0.912 };
   int lda = 2;
   double C[] = { -0.083, 0.219, 0.417, 0.817, -0.294, -0.683, -0.633, 0.831 };
   int ldc = 2;
   double C_expected[] = { -0.083, 0.219, 0.417, 0.817, -0.294, -0.683, -0.633, 0.831 };
   cblas_zherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zherk(case 1619) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zherk(case 1619) imag");
     };
   };
  };


  {
   int order = 101;
   int uplo = 122;
   int trans = 113;
   int N = 2;
   int K = 1;
   double alpha = -0.3;
   double beta = -1;
   double A[] = { -0.531, 0.187, -0.777, -0.329 };
   int lda = 2;
   double C[] = { -0.173, 0.833, 0.155, -0.52, -0.99, 0.28, 0.455, 0.481 };
   int ldc = 2;
   double C_expected[] = { 0.077921, 0.0, 0.155, -0.52, 0.8846808, -0.1840006, -0.668591, 0.0 };
   cblas_zherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zherk(case 1620) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zherk(case 1620) imag");
     };
   };
  };


  {
   int order = 102;
   int uplo = 122;
   int trans = 113;
   int N = 2;
   int K = 1;
   double alpha = -0.3;
   double beta = -1;
   double A[] = { -0.287, 0.068, 0.917, -0.449 };
   int lda = 1;
   double C[] = { -0.248, -0.608, -0.124, -0.718, -0.037, -0.115, 0.998, -0.551 };
   int ldc = 2;
   double C_expected[] = { 0.2219021, 0.0, 0.2121133, 0.7379521, -0.037, -0.115, -1.310747, 0.0 };
   cblas_zherk(order, uplo, trans, N, K, alpha, A, lda, beta, C, ldc);
   {
     int i;
     for (i = 0; i < 4; i++) {
       gsl_test_rel(C[2*i], C_expected[2*i], dbleps, "zherk(case 1621) real");
       gsl_test_rel(C[2*i+1], C_expected[2*i+1], dbleps, "zherk(case 1621) imag");
     };
   };
  };


}
