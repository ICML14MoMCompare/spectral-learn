#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

#include "tests.h"

void
test_trsm (void) {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { 0.773f, 0.069f, 0.45f, 0.189f };
   int lda = 2;
   float B[] = { -0.037f, 0.788f, 0.015f, 0.028f, -0.804f, -0.357f };
   int ldb = 3;
   float B_expected[] = { 0.0183269f, -0.419738f, -0.0564036f, -0.0444444f, 1.27619f, 0.566667f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1830)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { -0.13f, -0.832f, 0.426f, 0.195f };
   int lda = 2;
   float B[] = { 0.504f, 0.996f, 0.872f, -0.35f, 0.518f, -0.8f };
   int ldb = 3;
   float B_expected[] = { -0.06384f, -0.428093f, -0.06192f, 0.105f, -0.1554f, 0.24f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1831)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { 0.755f, -0.053f, -0.132f, -0.515f };
   int lda = 2;
   float B[] = { -0.735f, 0.494f, 0.072f, -0.882f, -0.112f, 0.904f };
   int ldb = 3;
   float B_expected[] = { 0.292053f, -0.196291f, -0.0286093f, -0.588643f, -0.0149311f, 0.533935f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1832)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { -0.88f, -0.555f, 0.642f, 0.751f };
   int lda = 2;
   float B[] = { -0.411f, 0.134f, 0.657f, 0.072f, -0.007f, -0.34f };
   int ldb = 3;
   float B_expected[] = { 0.1233f, -0.0402f, -0.1971f, -0.100759f, 0.0279084f, 0.228538f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1833)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { -0.478f, 0.938f, -0.731f, 0.25f };
   int lda = 2;
   float B[] = { -0.859f, -0.409f, -0.154f, -0.54f, 0.146f, -0.106f };
   int ldb = 2;
   float B_expected[] = { -1.2897f, 0.4908f, -1.08763f, 0.648f, -0.102894f, 0.1272f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1834)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { 0.953f, 0.249f, -0.451f, -0.781f };
   int lda = 2;
   float B[] = { -0.4f, -0.546f, 0.839f, 0.392f, -0.445f, -0.818f };
   int ldb = 2;
   float B_expected[] = { 0.193874f, 0.1638f, -0.304738f, -0.1176f, 0.244175f, 0.2454f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1835)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { 0.831f, -0.997f, -0.366f, 0.307f };
   int lda = 2;
   float B[] = { 0.157f, -0.02f, 0.57f, 0.309f, -0.159f, 0.266f };
   int ldb = 2;
   float B_expected[] = { -0.0566787f, -0.164523f, -0.205776f, -0.970224f, 0.0574007f, -0.0735227f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1836)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { -0.842f, 0.674f, 0.03f, 0.628f };
   int lda = 2;
   float B[] = { -0.426f, 0.806f, 0.299f, 0.626f, -0.471f, 0.208f };
   int ldb = 2;
   float B_expected[] = { 0.1278f, -0.327937f, -0.0897f, -0.127342f, 0.1413f, -0.157636f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1837)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { -0.095f, 0.301f, 0.168f, 0.934f, 0.107f, 0.068f, 0.384f, -0.201f, 0.116f };
   int lda = 3;
   float B[] = { 0.534f, 0.773f, -0.304f, -0.402f, 0.642f, -0.102f };
   int ldb = 3;
   float B_expected[] = { 1.68632f, -6.91104f, 2.39525f, -1.26947f, 1.77114f, 1.06409f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1838)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { -0.738f, -0.353f, -0.616f, 0.304f, 0.403f, 0.739f, 0.996f, 0.329f, 0.273f };
   int lda = 3;
   float B[] = { -0.436f, 0.074f, 0.273f, -0.609f, 0.858f, 0.993f };
   int ldb = 3;
   float B_expected[] = { 0.1308f, 0.0239724f, -0.0190428f, 0.1827f, -0.192907f, -0.0427986f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1839)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { -0.956f, 0.878f, 0.156f, 0.217f, 0.082f, -0.869f, 0.595f, 0.845f, 0.064f };
   int lda = 3;
   float B[] = { -0.744f, 0.662f, -0.31f, 0.811f, 0.257f, 0.98f };
   int ldb = 3;
   float B_expected[] = { -3.27779f, -17.3962f, 1.45312f, 7.92713f, 46.3978f, -4.59375f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1840)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { 0.313f, -0.316f, 0.836f, 0.359f, -0.415f, 0.154f, -0.948f, -0.596f, -0.799f };
   int lda = 3;
   float B[] = { 0.29f, -0.291f, 0.652f, 0.614f, 0.922f, -0.063f };
   int ldb = 3;
   float B_expected[] = { -0.261918f, -0.0292776f, -0.1956f, -0.0710273f, -0.265336f, 0.0189f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1841)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { -0.634f, 0.561f, 0.883f, -0.136f, 0.203f, -0.531f, 0.733f, -0.332f, 0.705f };
   int lda = 3;
   float B[] = { 0.133f, -0.843f, -0.179f, 0.94f, -0.656f, 0.645f };
   int ldb = 2;
   float B_expected[] = { 0.0629338f, -0.398896f, 0.306695f, -1.6564f, 0.358145f, -0.639766f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1842)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { 0.742f, -0.438f, 0.991f, 0.614f, 0.108f, -0.125f, 0.736f, -0.383f, 0.0f };
   int lda = 3;
   float B[] = { -0.792f, -0.033f, -0.723f, 0.885f, 0.336f, 0.584f };
   int ldb = 2;
   float B_expected[] = { 0.2376f, 0.0099f, 0.0710136f, -0.271579f, -0.248475f, -0.286501f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1843)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { 0.761f, 0.466f, 0.907f, -0.85f, -0.342f, -0.058f, -0.379f, -0.416f, 0.599f };
   int lda = 3;
   float B[] = { -0.238f, 0.013f, 0.473f, -0.626f, 0.912f, -0.003f };
   int ldb = 2;
   float B_expected[] = { 0.336709f, 0.329497f, 0.492375f, -0.549378f, -0.456761f, 0.0015025f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1844)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { 0.567f, -0.532f, -0.817f, 0.85f, -0.135f, 0.797f, 0.981f, -0.75f, 0.856f };
   int lda = 3;
   float B[] = { -0.705f, 0.326f, 0.184f, 0.079f, -0.173f, 0.125f };
   int ldb = 2;
   float B_expected[] = { 0.20253f, -0.125146f, -0.0965643f, 0.0061875f, 0.0519f, -0.0375f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1845)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { -0.859f, 0.563f, -0.61f, 0.2f };
   int lda = 2;
   float B[] = { -0.241f, -0.357f, -0.683f, -0.718f, 0.69f, -0.486f };
   int ldb = 3;
   float B_expected[] = { -0.0841676f, -0.12468f, -0.238533f, 1.31393f, -0.684026f, 1.40047f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1846)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { 0.157f, -0.741f, 0.844f, 0.206f };
   int lda = 2;
   float B[] = { 0.816f, -0.692f, 0.765f, -0.408f, 0.404f, 0.764f };
   int ldb = 3;
   float B_expected[] = { -0.2448f, 0.2076f, -0.2295f, -0.0589968f, 0.0326316f, -0.399259f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1847)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { 0.187f, 0.354f, -0.931f, 0.18f };
   int lda = 2;
   float B[] = { -0.215f, -0.645f, 0.847f, 0.014f, 0.83f, 0.761f };
   int ldb = 3;
   float B_expected[] = { 0.228752f, -5.85232f, -7.67336f, -0.0233333f, -1.38333f, -1.26833f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1848)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { -0.923f, 0.27f, -0.319f, -0.856f };
   int lda = 2;
   float B[] = { 0.391f, 0.01f, 0.429f, 0.685f, 0.332f, -0.643f };
   int ldb = 3;
   float B_expected[] = { -0.182855f, -0.0347724f, -0.0671649f, -0.2055f, -0.0996f, 0.1929f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1849)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { 0.724f, 0.201f, 0.87f, -0.638f };
   int lda = 2;
   float B[] = { -0.533f, 0.183f, 0.569f, 0.85f, 0.642f, -0.051f };
   int ldb = 2;
   float B_expected[] = { 0.220856f, 0.387218f, -0.235773f, 0.0781772f, -0.266022f, -0.386739f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1850)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { 0.291f, 0.244f, 0.931f, 0.857f };
   int lda = 2;
   float B[] = { 0.008f, -0.478f, -0.252f, -0.155f, 0.419f, -0.192f };
   int ldb = 2;
   float B_expected[] = { -0.0024f, 0.145634f, 0.0756f, -0.0238836f, -0.1257f, 0.174627f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1851)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { -0.634f, -0.529f, -0.344f, 0.375f };
   int lda = 2;
   float B[] = { -0.295f, 0.551f, 0.832f, 0.744f, -0.326f, 0.111f };
   int ldb = 2;
   float B_expected[] = { 0.228207f, -0.4408f, 0.890317f, -0.5952f, -0.0801653f, -0.0888f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1852)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { 0.641f, 0.989f, 0.998f, -0.005f };
   int lda = 2;
   float B[] = { -0.168f, 0.465f, 0.36f, 0.356f, -0.858f, 0.879f };
   int ldb = 2;
   float B_expected[] = { 0.188365f, -0.1395f, -0.0023748f, -0.1068f, 0.518199f, -0.2637f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1853)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { -0.638f, 0.389f, 0.997f, 0.909f, -0.598f, -0.43f, -0.345f, -0.897f, 0.119f };
   int lda = 3;
   float B[] = { 0.64f, 0.779f, -0.129f, 0.016f, 0.599f, -0.668f };
   int ldb = 3;
   float B_expected[] = { 0.904844f, 0.156956f, 0.32521f, 2.08405f, -0.910426f, 1.68403f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1854)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { 0.289f, 0.641f, -0.876f, -0.503f, -0.062f, -0.987f, 0.1f, -0.105f, 0.757f };
   int lda = 3;
   float B[] = { -0.285f, 0.285f, 0.219f, -0.986f, -0.0f, -0.605f };
   int ldb = 3;
   float B_expected[] = { 0.124319f, -0.150346f, -0.0657f, 0.339965f, 0.17914f, 0.1815f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1855)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { 0.524f, 0.018f, 0.292f, -0.573f, 0.866f, 0.749f, 0.99f, 0.101f, 0.871f };
   int lda = 3;
   float B[] = { 0.522f, -0.269f, -0.142f, -0.266f, -0.505f, -0.55f };
   int ldb = 3;
   float B_expected[] = { -0.298855f, -0.104554f, 0.400719f, 0.15229f, 0.275707f, -0.0156298f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1856)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { 0.283f, 0.62f, -0.387f, -0.739f, -0.599f, 0.114f, 0.552f, 0.083f, -0.976f };
   int lda = 3;
   float B[] = { 0.202f, 0.169f, 0.7f, 0.473f, 0.86f, -0.557f };
   int ldb = 3;
   float B_expected[] = { -0.0606f, -0.0954834f, -0.168624f, -0.1419f, -0.362864f, 0.275547f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1857)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { -0.185f, 0.178f, -0.22f, -0.645f, -0.585f, -0.342f, -0.594f, -0.141f, 0.944f };
   int lda = 3;
   float B[] = { 0.22f, -0.895f, -0.301f, -0.683f, -0.009f, -0.451f };
   int ldb = 2;
   float B_expected[] = { 0.888147f, -0.569939f, -0.155048f, -0.384802f, 0.00286017f, 0.143326f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1858)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { -0.145f, 0.746f, 0.541f, 0.584f, -0.394f, 0.371f, -0.172f, -0.601f, 0.542f };
   int lda = 3;
   float B[] = { 0.529f, 0.636f, 0.668f, 0.848f, -0.816f, -0.925f };
   int ldb = 2;
   float B_expected[] = { -0.0854817f, -0.0918985f, -0.0532752f, -0.0876225f, 0.2448f, 0.2775f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1859)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { 0.416f, -0.526f, -0.486f, -0.716f, 0.361f, 0.365f, -0.492f, 0.544f, 0.721f };
   int lda = 3;
   float B[] = { 0.25f, 0.746f, 0.55f, 0.836f, -0.024f, 0.226f };
   int ldb = 2;
   float B_expected[] = { -0.180288f, -0.537981f, -0.719755f, -1.47861f, 0.25283f, 0.291864f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1860)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha = -0.3f;
   float A[] = { -0.735f, -0.606f, -0.124f, 0.641f, -0.074f, -0.053f, -0.734f, 0.907f, 0.558f };
   int lda = 3;
   float B[] = { 0.623f, 0.392f, -0.808f, -0.022f, -0.665f, -0.616f };
   int ldb = 2;
   float B_expected[] = { -0.1869f, -0.1176f, 0.129139f, -0.0646656f, 0.183169f, 0.16679f };
   cblas_strsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strsm(case 1861)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { -0.584, -0.058, -0.964, -0.214 };
   int lda = 2;
   double B[] = { 0.073, -0.734, -0.058, -0.115, 0.513, 0.503 };
   int ldb = 3;
   double B_expected[] = { -0.0178370247087, 0.149492702599, 0.0332751888363, 0.053738317757, -0.239719626168, -0.235046728972 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1862)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { 0.251, -0.8, 0.365, 0.809 };
   int lda = 2;
   double B[] = { -0.632, -0.611, 0.9, 0.063, -0.652, -0.841 };
   int ldb = 3;
   double B_expected[] = { -0.05816, -0.11326, 0.02272, 0.0063, -0.0652, -0.0841 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1863)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { -0.833, 0.934, -0.608, 0.49 };
   int lda = 2;
   double B[] = { 0.336, -0.541, -0.729, -0.382, 0.741, 0.546 };
   int ldb = 3;
   double B_expected[] = { -0.0403361344538, 0.0649459783914, 0.0875150060024, -0.128008917853, 0.231810520126, 0.220018619693 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1864)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { 0.824, 0.907, 0.632, -0.348 };
   int lda = 2;
   double B[] = { 0.351, -0.301, 0.602, 0.873, 0.031, -0.2 };
   int ldb = 3;
   double B_expected[] = { 0.0351, -0.0301, 0.0602, 0.0651168, 0.0221232, -0.0580464 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1865)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { 0.427, 0.193, -0.959, -0.679 };
   int lda = 2;
   double B[] = { -0.646, 0.741, -0.339, 0.049, 0.734, -0.182 };
   int ldb = 2;
   double B_expected[] = { -0.3963857167, -0.10913107511, -0.0955986383061, -0.00721649484536, 0.232096380888, 0.0268041237113 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1866)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { 0.282, 0.766, -0.422, -0.518 };
   int lda = 2;
   double B[] = { 0.269, 0.211, -0.911, -0.685, -0.777, -0.919 };
   int ldb = 2;
   double B_expected[] = { 0.0358042, 0.0211, -0.120007, -0.0685, -0.1164818, -0.0919 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1867)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { -0.877, -0.818, 0.191, 0.468 };
   int lda = 2;
   double B[] = { 0.517, 0.669, 0.337, -0.579, 0.885, -0.677 };
   int ldb = 2;
   double B_expected[] = { -0.0589509692132, 0.039910485435, -0.0384264538198, -0.190882135095, -0.100912200684, -0.321038846495 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1868)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { 0.469, 0.115, 0.284, 0.139 };
   int lda = 2;
   double B[] = { 0.889, -0.002, -0.686, -0.256, 0.028, 0.371 };
   int ldb = 2;
   double B_expected[] = { 0.0889, -0.0104235, -0.0686, -0.017711, 0.0028, 0.036778 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1869)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { -0.218, -0.819, -0.523, 0.042, 0.545, -0.292, 0.283, 0.224, 0.247 };
   int lda = 3;
   double B[] = { 0.677, 0.153, -0.272, -0.226, 0.987, -0.216 };
   int ldb = 3;
   double B_expected[] = { -0.310550458716, -0.438607019611, -1.28619894589, 0.103669724771, 0.336890834105, 0.530329512606 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1870)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { 0.241, 0.561, 0.164, 0.486, 0.891, -0.508, -0.596, -0.074, 0.576 };
   int lda = 3;
   double B[] = { -0.325, 0.382, 0.368, 0.761, -0.349, 0.324 };
   int ldb = 3;
   double B_expected[] = { -0.0325, 0.0564325, 0.07079771, 0.0761, -0.0775921, -0.0194971868 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1871)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { 0.76, 0.58, -0.203, 0.053, 0.792, 0.355, -0.685, 0.449, -0.367 };
   int lda = 3;
   double B[] = { 0.861, -0.44, 0.842, -0.019, -0.382, -0.579 };
   int ldb = 3;
   double B_expected[] = { -0.0986936127734, 0.0745114634079, -0.229427792916, 0.149297547123, -0.137672708006, 0.157765667575 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1872)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { 0.802, 0.298, 0.159, 0.333, 0.515, 0.715, -0.32, -0.217, 0.301 };
   int lda = 3;
   double B[] = { -0.268, 0.1, -0.631, 0.472, 0.796, 0.278 };
   int ldb = 3;
   double B_expected[] = { -0.0457623309, -0.0036927, -0.0631, 0.0275803442, 0.0856326, 0.0278 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1873)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { -0.028, 0.186, -0.435, -0.747, 0.212, 0.257, 0.804, -0.595, 0.64 };
   int lda = 3;
   double B[] = { 0.729, -0.847, -0.577, 0.056, -0.493, 0.619 };
   int ldb = 2;
   double B_expected[] = { -2.60357142857, 3.025, -9.44607479784, 10.685259434, -5.58819230648, 6.23051463001 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1874)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { -0.74, -0.091, 0.484, 0.769, 0.91, 0.817, -0.26, 0.579, 0.393 };
   int lda = 3;
   double B[] = { 0.109, 0.969, -0.668, 0.544, 0.753, 0.796 };
   int ldb = 2;
   double B_expected[] = { 0.0109, 0.0969, -0.0751821, -0.0201161, 0.1216644359, 0.1164412219 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1875)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { 0.123, -0.328, -0.482, 0.083, -0.125, -0.712, -0.757, -0.009, 0.237 };
   int lda = 3;
   double B[] = { -0.18, 0.358, 0.839, -0.725, 0.73, -0.095 };
   int ldb = 2;
   double B_expected[] = { -5.40775366883, 2.28950005146, -2.42566413502, 0.808320675105, 0.308016877637, -0.0400843881857 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1876)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = 0.1;
   double A[] = { 0.255, -0.069, -0.137, -0.45, -0.24, 0.221, -0.509, -0.484, -0.131 };
   int lda = 3;
   double B[] = { -0.563, 0.993, 0.508, 0.771, 0.745, 0.233 };
   int ldb = 2;
   double B_expected[] = { -0.0437243505, 0.1074566983, 0.0343355, 0.0719507, 0.0745, 0.0233 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1877)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { -0.772, 0.079, -0.227, 0.998 };
   int lda = 2;
   double B[] = { -0.095, 0.012, -0.988, -0.722, 0.738, 0.05 };
   int ldb = 3;
   double B_expected[] = { -0.123056994819, 0.0155440414508, -1.27979274611, 0.733187878347, -0.740709398071, 0.051206039021 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1878)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { -0.045, 0.059, -0.61, -0.328 };
   int lda = 2;
   double B[] = { 0.302, -0.099, 0.521, 0.487, -0.961, 0.903 };
   int ldb = 3;
   double B_expected[] = { -0.302, 0.099, -0.521, -0.469182, 0.955159, -0.872261 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1879)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { -0.319, 0.642, 0.511, 0.762 };
   int lda = 2;
   double B[] = { 0.883, 0.987, 0.436, -0.783, 0.175, -0.973 };
   int ldb = 3;
   double B_expected[] = { 4.41405227952, 2.72615785879, 3.41221747752, 1.02755905512, -0.229658792651, 1.27690288714 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1880)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { 0.676, 0.038, 0.543, 0.296 };
   int lda = 2;
   double B[] = { 0.804, -0.28, -0.318, 0.382, -0.165, -0.007 };
   int ldb = 3;
   double B_expected[] = { -0.596574, 0.190405, 0.314199, -0.382, 0.165, 0.007 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1881)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { -0.722, -0.355, -0.14, -0.146 };
   int lda = 2;
   double B[] = { -0.44, 0.751, -0.995, 0.625, 0.16, -0.127 };
   int ldb = 2;
   double B_expected[] = { -0.609418282548, 5.72820931203, -1.37811634349, 5.60230334307, 0.221606648199, -1.08236253937 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1882)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { 0.817, -0.619, 0.548, 0.064 };
   int lda = 2;
   double B[] = { -0.756, -0.169, 0.429, -0.789, 0.79, 0.479 };
   int ldb = 2;
   double B_expected[] = { 0.756, -0.245288, -0.429, 1.024092, -0.79, -0.04608 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1883)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { 0.496, -0.734, -0.679, -0.697 };
   int lda = 2;
   double B[] = { -0.483, -0.508, -0.819, 0.237, 0.852, -0.512 };
   int ldb = 2;
   double B_expected[] = { -0.104772180312, -0.728837876614, 2.1543973018, 0.340028694405, -2.80479705651, -0.734576757532 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1884)");
     }
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { 0.049, 0.079, -0.8, -0.762 };
   int lda = 2;
   double B[] = { 0.426, 0.094, 0.794, -0.098, 0.442, -0.991 };
   int ldb = 2;
   double B_expected[] = { -0.418574, -0.094, -0.801742, 0.098, -0.520289, 0.991 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1885)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { -0.974, 0.848, -0.765, 0.528, -0.693, 0.252, -0.135, -0.507, 0.954 };
   int lda = 3;
   double B[] = { 0.395, 0.791, -0.787, 0.636, 0.271, -0.905 };
   int ldb = 3;
   double B_expected[] = { 1.01254427581, 1.4413950829, 0.824947589099, 0.548697105717, 0.736012415258, 0.948637316562 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1886)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { 0.919, -0.513, -0.38, 0.587, -0.862, 0.598, 0.714, 0.726, 0.491 };
   int lda = 3;
   double B[] = { -0.056, -0.802, -0.962, 0.656, -0.195, -0.679 };
   int ldb = 3;
   double B_expected[] = { 0.537869412, 0.226724, 0.962, -0.506244546, -0.211042, 0.679 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1887)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { -0.287, -0.009, -0.989, -0.062, 0.714, -0.293, -0.875, 0.371, 0.728 };
   int lda = 3;
   double B[] = { -0.14, -0.969, 0.702, -0.317, -0.739, -0.518 };
   int ldb = 3;
   double B_expected[] = { -0.487804878049, 1.31478445037, -2.22062403761, -1.10452961672, 0.939102470256, -1.09460224052 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1888)");
     }
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { 0.236, 0.605, 0.338, -0.926, 0.362, 0.562, -0.554, 0.076, 0.85 };
   int lda = 3;
   double B[] = { 0.113, 0.604, 0.859, 0.216, -0.6, -0.048 };
   int ldb = 3;
   double B_expected[] = { -0.113, -0.708638, -0.867745512, -0.216, 0.399984, -0.102062784 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1889)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { -0.476, 0.428, -0.214, -0.889, -0.526, -0.704, 0.458, -0.479, 0.077 };
   int lda = 3;
   double B[] = { 0.124, -0.007, 0.452, 0.966, 0.42, 0.369 };
   int ldb = 2;
   double B_expected[] = { -15.8695808776, -16.2060574143, 5.8264777048, 6.20050861686, -5.45454545455, -4.79220779221 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1890)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { 0.748, 0.242, -0.964, 0.422, -0.78, -0.595, -0.926, -0.474, 0.947 };
   int lda = 3;
   double B[] = { 0.242, -0.553, -0.899, -0.714, -0.084, -0.609 };
   int ldb = 2;
   double B_expected[] = { -0.560396352, 0.693808948, 0.938816, 1.002666, 0.084, 0.609 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1891)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { 0.808, -0.074, 0.359, -0.172, -0.934, -0.67, 0.92, -0.617, -0.383 };
   int lda = 3;
   double B[] = { 0.079, 0.978, 0.82, 0.444, -0.597, -0.64 };
   int ldb = 2;
   double B_expected[] = { -0.0977722772277, -1.2103960396, 0.885690737168, 0.571273347892, -3.19977295412, -3.80492250994 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1892)");
     }
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha = -1;
   double A[] = { 0.786, 0.922, -0.763, 0.498, -0.082, 0.538, 0.742, -0.391, -0.255 };
   int lda = 3;
   double B[] = { 0.911, 0.066, 0.895, 0.255, -0.547, -0.805 };
   int ldb = 2;
   double B_expected[] = { -0.911, -0.066, -0.055058, -0.194148, -0.118471796, 0.859093624 };
   cblas_dtrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrsm(case 1893)");
     }
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0.0f, 1.0f};
   float A[] = { -0.362f, -0.457f, -0.347f, -0.203f, -0.517f, 0.462f, 0.572f, 0.521f };
   int lda = 2;
   float B[] = { 0.118f, -0.593f, 0.773f, 0.053f, -0.419f, -0.096f, 0.846f, -0.311f, -0.364f, 0.161f, -0.496f, -0.393f };
   int ldb = 3;
   float B_expected[] = { -1.58885f, 0.58489f, -0.628497f, -0.878921f, 0.701485f, 1.08554f, 1.03347f, 0.537701f, -0.470639f, -0.207688f, -0.056162f, -0.815978f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1894) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1894) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0.0f, 1.0f};
   float A[] = { 0.845f, -0.654f, -0.43f, -0.834f, 0.206f, 0.414f, 0.761f, 0.961f };
   int lda = 2;
   float B[] = { 0.069f, 0.005f, -0.419f, 0.806f, 0.857f, 0.669f, 0.942f, 0.657f, 0.52f, 0.19f, -0.609f, -0.305f };
   int ldb = 3;
   float B_expected[] = { -1.07314f, -0.073878f, -1.32138f, -0.35386f, -0.029944f, 0.8495f, -0.657f, 0.942f, -0.19f, 0.52f, 0.305f, -0.609f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1895) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1895) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0.0f, 1.0f};
   float A[] = { -0.87f, 0.218f, 0.813f, 0.575f, -0.848f, 0.7f, -0.311f, 0.374f };
   int lda = 2;
   float B[] = { 0.117f, 0.758f, -0.189f, -0.768f, 0.857f, -0.269f, 0.796f, -0.592f, -0.499f, 0.977f, 0.643f, 0.282f };
   int ldb = 3;
   float B_expected[] = { 0.851499f, 0.0788813f, -0.881826f, -0.00372192f, -0.0586805f, -0.999761f, -1.37808f, -2.51525f, 2.45259f, 2.57925f, 1.0972f, 1.8459f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1896) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1896) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0.0f, 1.0f};
   float A[] = { -0.697f, -0.668f, 0.746f, -0.818f, 0.651f, 0.275f, -0.702f, -0.615f };
   int lda = 2;
   float B[] = { -0.876f, 0.842f, -0.848f, 0.901f, 0.75f, 0.361f, -0.702f, 0.039f, -0.41f, 0.541f, 0.489f, 0.025f };
   int ldb = 3;
   float B_expected[] = { -0.842f, -0.876f, -0.901f, -0.848f, -0.361f, 0.75f, 0.268242f, 0.099826f, -0.187649f, 0.389823f, 0.416261f, 0.100025f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1897) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1897) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0.0f, 1.0f};
   float A[] = { -0.394f, -0.65f, -0.082f, -0.632f, -0.53f, 0.483f, 0.149f, -0.192f };
   int lda = 2;
   float B[] = { -0.691f, 0.732f, 0.976f, 0.073f, 0.607f, 0.918f, -0.918f, 0.67f, 0.37f, -0.344f, -0.114f, -0.62f };
   int ldb = 2;
   float B_expected[] = { -1.39367f, -3.05481f, -3.35679f, 2.2248f, 4.33836f, -1.06673f, 1.29393f, -4.49373f, -1.89826f, 2.23995f, 1.93461f, 1.72783f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1898) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1898) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0.0f, 1.0f};
   float A[] = { 0.45f, 0.972f, -0.051f, 0.71f, -0.127f, -0.274f, 0.152f, 0.789f };
   int lda = 2;
   float B[] = { 0.683f, -0.915f, -0.773f, 0.088f, -0.28f, 0.17f, 0.818f, 0.293f, -0.551f, 0.365f, 0.899f, 0.257f };
   int ldb = 2;
   float B_expected[] = { 1.11563f, 0.560717f, -0.088f, -0.773f, -0.431343f, -0.256396f, -0.293f, 0.818f, -0.643965f, -0.507245f, -0.257f, 0.899f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1899) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1899) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0.0f, 1.0f};
   float A[] = { -0.993f, -0.028f, -0.547f, -0.251f, 0.781f, -0.315f, 0.865f, 0.229f };
   int lda = 2;
   float B[] = { 0.578f, 0.73f, -0.931f, 0.288f, 0.048f, 0.508f, -0.168f, 0.655f, 0.92f, -0.26f, 0.485f, 0.05f };
   int ldb = 2;
   float B_expected[] = { 0.718162f, -0.602325f, -0.0323644f, -1.24023f, 0.509813f, -0.0627138f, -0.410611f, 0.0227614f, -0.287729f, -0.918372f, -0.000635986f, -0.10338f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1900) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1900) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0.0f, 1.0f};
   float A[] = { 0.131f, -0.494f, 0.615f, -0.089f, -0.22f, -0.874f, 0.677f, 0.074f };
   int lda = 2;
   float B[] = { -0.276f, 0.539f, 0.647f, 0.986f, -0.34f, 0.983f, -0.819f, 0.144f, 0.361f, 0.561f, 0.178f, -0.433f };
   int ldb = 2;
   float B_expected[] = { -0.539f, -0.276f, -0.629951f, 0.768769f, -0.983f, -0.34f, 0.490805f, -0.697387f, -0.561f, 0.361f, 0.745886f, -0.093944f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1901) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1901) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0.0f, 1.0f};
   float A[] = { -0.123f, -0.726f, -0.011f, 0.245f, -0.205f, 0.77f, -0.81f, -0.973f, 0.354f, -0.835f, 0.552f, 0.396f, -0.524f, -0.204f, -0.814f, 0.284f, -0.976f, -0.835f };
   int lda = 3;
   float B[] = { -0.42f, 0.976f, -0.845f, 0.651f, -0.44f, -0.862f, 0.137f, 0.066f, -0.63f, 0.482f, -0.187f, 0.724f };
   int ldb = 3;
   float B_expected[] = { 0.783777f, -1.21156f, 0.66205f, -1.40548f, 0.886226f, 0.0391664f, -0.168468f, -0.119451f, 0.378144f, -0.774828f, 0.708857f, -0.807468f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1902) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1902) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0.0f, 1.0f};
   float A[] = { 0.921f, 0.167f, -0.41f, 0.578f, -0.372f, 0.106f, 0.551f, 0.668f, 0.295f, 0.855f, -0.167f, 0.976f, -0.782f, -0.777f, 0.278f, -0.98f, 0.038f, -0.832f };
   int lda = 3;
   float B[] = { 0.459f, 0.06f, 0.387f, 0.871f, -0.366f, 0.926f, 0.236f, -0.889f, 0.619f, 0.319f, -0.709f, 0.884f };
   int ldb = 3;
   float B_expected[] = { -0.06f, 0.459f, -0.630298f, 0.60987f, -0.409693f, 0.528127f, 0.889f, 0.236f, 0.181898f, 0.201918f, -0.300827f, -0.859254f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1903) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1903) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0.0f, 1.0f};
   float A[] = { -0.691f, -0.817f, 0.954f, -0.969f, -0.574f, -0.026f, 0.992f, 0.529f, 0.135f, -0.413f, -0.314f, -0.859f, -0.284f, -0.849f, 0.781f, 0.534f, -0.018f, 0.282f };
   int lda = 3;
   float B[] = { -0.028f, -0.429f, 0.066f, -0.854f, -0.316f, 0.514f, -0.465f, -0.857f, 0.286f, 0.415f, -0.486f, 0.538f };
   int ldb = 3;
   float B_expected[] = { 6.83575f, 2.7232f, 3.79999f, 5.15624f, -1.00015f, 1.88653f, 5.42614f, 2.69261f, 2.30584f, 3.85628f, -1.59513f, 2.00962f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1904) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1904) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0.0f, 1.0f};
   float A[] = { -0.839f, -0.318f, 0.175f, 0.72f, -0.683f, 0.395f, -0.279f, 0.151f, -0.71f, 0.445f, 0.533f, -0.38f, -0.749f, -0.833f, 0.871f, -0.426f, 0.195f, 0.889f };
   int lda = 3;
   float B[] = { 0.804f, -0.346f, 0.234f, 0.782f, 0.033f, 0.581f, 0.981f, -0.68f, 0.919f, -0.758f, 0.152f, -0.503f };
   int ldb = 3;
   float B_expected[] = { -0.20395f, 0.376748f, -0.290007f, -0.042249f, -0.581f, 0.033f, 1.15245f, 1.75457f, 0.255135f, 1.00089f, 0.503f, 0.152f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1905) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1905) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0.0f, 1.0f};
   float A[] = { -0.365f, -0.662f, 0.188f, -0.571f, 0.082f, 0.192f, -0.833f, -0.958f, 0.159f, -0.203f, 0.481f, 0.08f, -0.954f, 0.681f, -0.015f, 0.146f, -0.352f, -0.068f };
   int lda = 3;
   float B[] = { 0.779f, -0.691f, -0.516f, 0.148f, 0.721f, 0.217f, -0.976f, -0.963f, 0.532f, -0.366f, 0.176f, 0.4f };
   int ldb = 2;
   float B_expected[] = { -1.34375f, 0.302916f, 0.692272f, 0.158126f, -2.93098f, -5.71682f, 3.87247f, 3.8052f, 3.25028f, -6.53201f, -2.34332f, 2.30748f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1906) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1906) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0.0f, 1.0f};
   float A[] = { 0.779f, 0.065f, -0.616f, -0.245f, 0.823f, 0.689f, 0.06f, -0.164f, 0.768f, -0.727f, 0.897f, -0.556f, -0.875f, 0.862f, 0.863f, -0.085f, 0.171f, 0.063f };
   int lda = 3;
   float B[] = { -0.621f, 0.428f, 0.096f, 0.711f, 0.416f, -0.684f, 0.806f, 0.491f, 0.037f, -0.776f, -0.312f, 0.391f };
   int ldb = 2;
   float B_expected[] = { -0.428f, -0.621f, -0.711f, 0.096f, 0.811524f, 0.383068f, -0.464084f, 0.683636f, -0.866708f, -0.399047f, -0.587978f, -0.244543f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1907) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1907) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {0.0f, 1.0f};
   float A[] = { -0.349f, 0.224f, -0.323f, 0.457f, 0.081f, 0.443f, 0.809f, 0.037f, -0.543f, 0.554f, 0.779f, 0.632f, -0.852f, -0.148f, -0.649f, -0.78f, 0.469f, -0.515f };
   int lda = 3;
   float B[] = { 0.162f, 0.754f, -0.978f, -0.097f, 0.986f, 0.943f, 0.676f, 0.718f, 0.204f, 0.264f, -0.124f, -0.73f };
   int ldb = 2;
   float B_expected[] = { 0.0811068f, 1.92921f, -3.74716f, 1.18561f, 1.80842f, -0.638944f, 0.528341f, 1.20828f, -0.471728f, -0.083028f, 0.837267f, 0.654994f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1908) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1908) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {0.0f, 1.0f};
   float A[] = { -0.469f, -0.164f, -0.792f, -0.454f, 0.206f, 0.785f, 0.504f, -0.561f, 0.205f, 0.463f, -0.8f, 0.803f, 0.283f, 0.131f, 0.576f, -0.431f, 0.297f, -0.415f };
   int lda = 3;
   float B[] = { -0.364f, 0.853f, 0.056f, -0.78f, 0.05f, 0.223f, -0.166f, -0.097f, 0.24f, 0.721f, 0.023f, 0.508f };
   int ldb = 2;
   float B_expected[] = { -1.3696f, 0.527133f, 0.554099f, 0.524136f, -0.60708f, 0.820963f, -0.290931f, 0.260324f, -0.721f, 0.24f, -0.508f, 0.023f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1909) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1909) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {1.0f, 0.0f};
   float A[] = { 0.917f, 0.367f, -0.115f, -0.321f, -0.811f, -0.563f, 0.78f, -0.742f };
   int lda = 2;
   float B[] = { 0.797f, 0.166f, 0.737f, -0.685f, 0.677f, -0.04f, -0.652f, 0.327f, 0.094f, -0.656f, 0.496f, -0.646f };
   int ldb = 3;
   float B_expected[] = { 0.811592f, -0.143789f, 0.435059f, -0.92112f, 0.621302f, -0.292277f, -0.710488f, 0.0561583f, 0.694329f, -0.137285f, 0.752465f, 0.100199f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1910) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1910) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {1.0f, 0.0f};
   float A[] = { -0.384f, 0.783f, -0.086f, -0.649f, -0.574f, 0.216f, -0.809f, -0.608f };
   int lda = 2;
   float B[] = { 0.067f, -0.183f, -0.524f, 0.77f, 0.169f, 0.769f, -0.982f, -0.522f, -0.051f, -0.129f, 0.595f, 0.56f };
   int ldb = 3;
   float B_expected[] = { 0.067f, -0.183f, -0.524f, 0.77f, 0.169f, 0.769f, -0.857471f, -0.494255f, -0.595794f, -0.402856f, 0.110453f, 0.735815f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1911) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1911) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {1.0f, 0.0f};
   float A[] = { 0.228f, -0.644f, 0.731f, 0.458f, 0.051f, -0.725f, 0.731f, 0.537f };
   int lda = 2;
   float B[] = { -0.588f, 0.01f, -0.009f, -0.374f, 0.422f, 0.758f, -0.428f, 0.263f, 0.659f, 0.171f, -0.239f, 0.968f };
   int ldb = 3;
   float B_expected[] = { -0.232749f, -1.39168f, -0.124158f, 0.287962f, -1.55821f, 0.0298572f, -0.208619f, 0.513035f, 0.697138f, -0.278198f, 0.419466f, 1.01607f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1912) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1912) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {1.0f, 0.0f};
   float A[] = { -0.321f, 0.761f, 0.809f, -0.017f, -0.009f, -0.975f, 0.057f, 0.396f };
   int lda = 2;
   float B[] = { 0.377f, 0.776f, -0.686f, -0.561f, 0.29f, 0.601f, 0.755f, 0.518f, 0.313f, -0.394f, 0.945f, 0.395f };
   int ldb = 3;
   float B_expected[] = { -0.121255f, 1.51679f, -0.299033f, -0.259371f, -0.08662f, 1.52593f, 0.755f, 0.518f, 0.313f, -0.394f, 0.945f, 0.395f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1913) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1913) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {1.0f, 0.0f};
   float A[] = { 0.186f, 0.818f, -0.142f, -0.376f, 0.332f, 0.746f, 0.413f, -0.151f };
   int lda = 2;
   float B[] = { -0.374f, -0.787f, 0.221f, -0.104f, 0.74f, -0.548f, 0.88f, -0.66f, 0.65f, 0.046f, -0.839f, -0.783f };
   int ldb = 2;
   float B_expected[] = { -1.01366f, 0.226724f, 1.10152f, 1.79962f, -0.441403f, -1.00501f, 0.588898f, 0.222456f, 0.225271f, -0.743398f, -2.5862f, -2.65075f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1914) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1914) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {1.0f, 0.0f};
   float A[] = { -0.574f, 0.018f, -0.584f, -0.184f, 0.41f, 0.075f, 0.92f, 0.022f };
   int lda = 2;
   float B[] = { 0.524f, -0.234f, 0.198f, 0.079f, -0.449f, -0.433f, -0.14f, -0.201f, -0.242f, -0.368f, -0.298f, 0.693f };
   int ldb = 2;
   float B_expected[] = { 0.524f, -0.234f, -0.03439f, 0.13564f, -0.449f, -0.433f, 0.011615f, 0.010205f, -0.242f, -0.368f, -0.22638f, 0.86203f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1915) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1915) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {1.0f, 0.0f};
   float A[] = { 0.422f, -0.38f, 0.919f, -0.229f, -0.849f, -0.19f, 0.02f, -0.181f };
   int lda = 2;
   float B[] = { 0.971f, -0.339f, 0.203f, 0.083f, 0.461f, -0.623f, 0.334f, 0.653f, 0.694f, 0.42f, 0.239f, -0.061f };
   int ldb = 2;
   float B_expected[] = { 3.06394f, -0.745692f, -0.330599f, 1.15808f, 8.0252f, -0.902398f, -3.36278f, 2.21688f, 0.70369f, -0.872941f, 0.477097f, 1.26772f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1916) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1916) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {1.0f, 0.0f};
   float A[] = { 0.564f, 0.483f, -0.635f, 0.84f, 0.238f, 0.35f, 0.96f, 0.397f };
   int lda = 2;
   float B[] = { 0.963f, -0.513f, 0.989f, 0.404f, -0.352f, 0.924f, 0.052f, -0.059f, -0.771f, 0.341f, -0.566f, -0.844f };
   int ldb = 2;
   float B_expected[] = { 1.93037f, -1.08722f, 0.989f, 0.404f, -0.36854f, 0.842855f, 0.052f, -0.059f, -1.83937f, 0.2805f, -0.566f, -0.844f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1917) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1917) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {1.0f, 0.0f};
   float A[] = { 0.182f, 0.699f, 0.303f, -0.273f, -0.363f, 0.02f, 0.991f, -0.206f, -0.347f, 0.269f, -0.384f, 0.797f, 0.392f, -0.966f, 0.347f, 0.87f, 0.016f, -0.097f };
   int lda = 3;
   float B[] = { 0.587f, 0.875f, -0.848f, 0.154f, -0.887f, -0.709f, 0.824f, -0.895f, 0.159f, 0.933f, -0.011f, -0.393f };
   int ldb = 3;
   float B_expected[] = { -14.9753f, 2.31554f, 0.613295f, 24.1527f, 5.64728f, -10.0758f, -3.79783f, -5.34545f, -5.38045f, 2.99977f, 3.92602f, -0.760993f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1918) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1918) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {1.0f, 0.0f};
   float A[] = { 0.637f, -0.57f, 0.322f, -0.303f, 0.618f, 0.261f, 0.654f, -0.238f, 0.66f, -0.485f, 0.223f, -0.196f, -0.252f, 0.929f, -0.012f, 0.965f, 0.783f, 0.489f };
   int lda = 3;
   float B[] = { 0.894f, 0.93f, 0.648f, 0.914f, 0.7f, -0.138f, 0.63f, -0.173f, -0.671f, -0.327f, -0.922f, 0.816f };
   int ldb = 3;
   float B_expected[] = { -0.0695574f, 0.64143f, 0.518948f, 1.08197f, 0.7f, -0.138f, 1.8231f, -0.404044f, -0.62533f, -0.68968f, -0.922f, 0.816f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1919) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1919) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {1.0f, 0.0f};
   float A[] = { 0.274f, 0.721f, -0.445f, 0.14f, 0.023f, -0.945f, 0.859f, -0.522f, -0.227f, 0.722f, 0.165f, 0.969f, -0.212f, -0.816f, 0.908f, -0.652f, -0.208f, -0.229f };
   int lda = 3;
   float B[] = { 0.011f, -0.818f, 0.067f, -0.191f, -0.911f, 0.84f, -0.162f, -0.951f, -0.502f, -0.21f, 0.492f, 0.767f };
   int ldb = 3;
   float B_expected[] = { -0.986296f, -0.390076f, -0.910328f, -1.26205f, -3.05033f, 0.930902f, -1.22716f, -0.241667f, -1.07925f, -0.600129f, -2.84941f, 5.27338f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1920) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1920) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {1.0f, 0.0f};
   float A[] = { -0.186f, 0.118f, -0.545f, 0.784f, 0.057f, 0.39f, 0.77f, -0.518f, -0.97f, 0.271f, 0.488f, 0.637f, -0.482f, -0.993f, -0.797f, -0.945f, 0.257f, 0.3f };
   int lda = 3;
   float B[] = { -0.783f, 0.649f, 0.698f, 0.046f, -0.153f, 0.473f, -0.996f, -0.211f, 0.84f, 0.201f, -0.457f, 0.918f };
   int ldb = 3;
   float B_expected[] = { -0.783f, 0.649f, 0.964728f, -0.859324f, 0.406086f, 0.235086f, -0.996f, -0.211f, 1.71622f, -0.152458f, 0.78435f, 1.32759f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1921) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1921) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {1.0f, 0.0f};
   float A[] = { -0.681f, -0.342f, -0.195f, -0.053f, 0.016f, -0.191f, 0.989f, -0.718f, -0.59f, 0.646f, -0.41f, -0.809f, -0.359f, -0.783f, -0.902f, 0.917f, -0.703f, 0.795f };
   int lda = 3;
   float B[] = { -0.27f, 0.037f, 0.349f, 0.36f, -0.293f, 0.128f, -0.481f, -0.834f, -0.815f, -0.6f, 0.728f, 0.122f };
   int ldb = 2;
   float B_expected[] = { -0.69977f, -2.39368f, 2.17354f, 1.74016f, 0.260417f, -1.25151f, 0.175881f, 1.93577f, 0.085191f, 0.949825f, -0.368302f, -0.590043f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1922) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1922) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {1.0f, 0.0f};
   float A[] = { -0.132f, 0.33f, 0.357f, 0.32f, 0.833f, -0.111f, -0.192f, -0.643f, -0.622f, -0.663f, -0.58f, 0.423f, -0.874f, 0.86f, -0.281f, -0.992f, 0.055f, 0.137f };
   int lda = 3;
   float B[] = { 0.104f, -0.906f, -0.712f, 0.103f, -0.474f, -0.591f, 0.073f, -0.906f, -0.261f, -0.391f, 0.881f, -0.345f };
   int ldb = 2;
   float B_expected[] = { 0.126148f, -1.31009f, -0.0285057f, -0.554776f, -0.159469f, -0.959783f, 0.662801f, -0.128993f, -0.261f, -0.391f, 0.881f, -0.345f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1923) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1923) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {1.0f, 0.0f};
   float A[] = { -0.353f, -0.581f, -0.648f, 0.894f, 0.825f, -0.23f, -0.529f, 0.213f, 0.568f, 0.296f, 0.372f, 0.442f, 0.515f, -0.409f, 0.222f, -0.246f, -0.524f, 0.318f };
   int lda = 3;
   float B[] = { -0.467f, 0.632f, 0.672f, 0.777f, -0.609f, 0.511f, -0.991f, 0.311f, -0.617f, -0.732f, -0.585f, 0.152f };
   int ldb = 2;
   float B_expected[] = { -0.437806f, -1.06979f, -1.49004f, 0.251317f, -2.40924f, 1.62379f, -1.09482f, 3.75003f, -1.80514f, -2.07012f, -4.8059f, -0.418185f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1924) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1924) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {1.0f, 0.0f};
   float A[] = { -0.201f, -0.918f, -0.514f, -0.889f, -0.948f, 0.34f, 0.818f, 0.557f, 0.341f, 0.484f, 0.235f, 0.561f, 0.874f, -0.342f, -0.411f, -0.975f, -0.85f, -0.621f };
   int lda = 3;
   float B[] = { -0.389f, -0.252f, 0.322f, -0.763f, -0.839f, -0.744f, -0.946f, -0.312f, 0.051f, -0.686f, -0.626f, -0.043f };
   int ldb = 2;
   float B_expected[] = { -0.389f, -0.252f, 0.322f, -0.763f, -0.814918f, -1.21935f, -0.102185f, -0.417924f, -0.896001f, -0.04892f, -0.790606f, -0.720266f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1925) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1925) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3f, 0.1f};
   float A[] = { 0.571f, 0.12f, 0.665f, 0.425f, -0.977f, -0.772f, -0.944f, -0.154f };
   int lda = 2;
   float B[] = { -0.357f, -0.213f, 0.57f, 0.134f, 0.089f, 0.046f, 0.027f, 0.825f, -0.127f, 0.658f, -0.332f, 0.247f };
   int ldb = 3;
   float B_expected[] = { 0.205417f, 0.092557f, -0.315204f, -0.0368205f, -0.0507703f, -0.0192512f, 0.238158f, 0.270895f, -0.257649f, 0.296502f, -0.140106f, 0.100105f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1926) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1926) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3f, 0.1f};
   float A[] = { 0.051f, 0.966f, 0.04f, -0.765f, 0.276f, -0.798f, 0.766f, -0.37f };
   int lda = 2;
   float B[] = { 0.532f, 0.59f, 0.305f, 0.443f, 0.036f, 0.655f, -0.145f, -0.864f, -0.483f, -0.45f, -0.327f, -0.365f };
   int ldb = 3;
   float B_expected[] = { -0.2186f, -0.1238f, -0.1358f, -0.1024f, -0.0763f, -0.1929f, 0.043937f, 0.416881f, 0.116996f, 0.194683f, -0.0099165f, 0.142885f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1927) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1927) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3f, 0.1f};
   float A[] = { -0.163f, -0.238f, -0.032f, 0.494f, 0.863f, 0.96f, 0.669f, 0.415f };
   int lda = 2;
   float B[] = { -0.724f, -0.682f, 0.034f, 0.352f, 0.42f, 0.253f, 0.186f, -0.061f, 0.278f, -0.764f, -0.484f, 0.051f };
   int ldb = 3;
   float B_expected[] = { -0.532386f, -1.09223f, -1.1606f, 1.43429f, 1.04476f, 0.724237f, -0.0783541f, 0.00655162f, -0.179639f, 0.27272f, 0.193877f, 0.0250509f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1928) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1928) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3f, 0.1f};
   float A[] = { 0.533f, 0.575f, 0.808f, -0.631f, 0.185f, 0.296f, -0.757f, -0.279f };
   int lda = 2;
   float B[] = { -0.744f, -0.881f, -0.594f, 0.629f, -0.924f, 0.017f, -0.089f, -0.052f, 0.959f, -0.486f, 0.39f, -0.378f };
   int ldb = 3;
   float B_expected[] = { 0.303415f, 0.198103f, 0.0879903f, -0.363588f, 0.245042f, -0.149137f, 0.0319f, 0.0067f, -0.2391f, 0.2417f, -0.0792f, 0.1524f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1929) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1929) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3f, 0.1f};
   float A[] = { -0.001f, -0.948f, -0.97f, -0.285f, -0.664f, -0.977f, -0.746f, 0.192f };
   int lda = 2;
   float B[] = { 0.997f, -0.852f, 0.87f, -0.955f, 0.007f, -0.071f, -0.263f, -0.077f, -0.856f, 0.228f, -0.81f, 0.476f };
   int ldb = 2;
   float B_expected[] = { 0.375027f, 0.225237f, -0.432345f, -0.0987217f, 0.0232012f, -0.00529874f, -0.112225f, 0.0682749f, -0.162707f, -0.246664f, 0.267117f, 0.237712f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1930) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1930) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3f, 0.1f};
   float A[] = { 0.804f, 0.476f, -0.898f, -0.966f, 0.51f, -0.346f, 0.622f, -0.749f };
   int lda = 2;
   float B[] = { -0.964f, 0.453f, 0.799f, -0.949f, -0.055f, 0.803f, 0.99f, -0.162f, 0.913f, -0.081f, -0.057f, 0.014f };
   int ldb = 2;
   float B_expected[] = { 0.2439f, -0.2323f, -0.349565f, 0.398684f, -0.0638f, -0.2464f, -0.333516f, 0.295339f, -0.2658f, 0.1156f, 0.191256f, 0.0231108f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1931) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1931) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3f, 0.1f};
   float A[] = { -0.144f, -0.611f, -0.054f, 0.618f, 0.213f, 0.49f, -0.465f, -0.488f };
   int lda = 2;
   float B[] = { -0.225f, -0.663f, 0.073f, -0.379f, -0.297f, 0.822f, -0.038f, -0.935f, -0.81f, 0.885f, -0.065f, 0.412f };
   int ldb = 2;
   float B_expected[] = { 0.287563f, -0.439427f, 0.113582f, -0.141015f, -0.375321f, -0.339988f, 0.189826f, -0.395838f, -0.655583f, 0.0702722f, -0.117522f, 0.15645f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1932) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1932) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3f, 0.1f};
   float A[] = { -0.367f, 0.502f, -0.309f, 0.404f, 0.531f, -0.188f, 0.181f, 0.583f };
   int lda = 2;
   float B[] = { 0.861f, -0.648f, 0.906f, -0.402f, 0.455f, 0.412f, 0.34f, -0.248f, 0.107f, 0.507f, 0.088f, -0.593f };
   int ldb = 2;
   float B_expected[] = { -0.350389f, 0.252194f, -0.2316f, 0.2112f, -0.245348f, -0.0757932f, -0.0772f, 0.1084f, -0.148061f, -0.0704181f, 0.0329f, 0.1867f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1933) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1933) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3f, 0.1f};
   float A[] = { -0.476f, 0.202f, -0.66f, 0.774f, -0.943f, -0.99f, -0.035f, 0.901f, -0.742f, -0.085f, -0.335f, -0.591f, 0.799f, 0.515f, 0.753f, 0.76f, -0.042f, -0.011f };
   int lda = 3;
   float B[] = { 0.025f, -0.976f, -0.44f, 0.741f, -0.126f, 0.527f, 0.743f, 0.216f, 0.661f, -0.071f, 0.564f, -0.093f };
   int ldb = 3;
   float B_expected[] = { -6.73789f, 0.501263f, -2.62173f, -2.22684f, -0.664138f, 3.89034f, 4.11106f, 5.79368f, -1.20958f, 3.39994f, 4.05469f, -0.945199f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1934) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1934) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3f, 0.1f};
   float A[] = { -0.999f, 0.418f, 0.687f, 0.6f, 0.106f, -0.737f, -0.165f, 0.263f, 0.998f, -0.092f, 0.555f, -0.671f, -0.162f, -0.814f, 0.317f, 0.582f, 0.302f, -0.48f };
   int lda = 3;
   float B[] = { 0.699f, 0.128f, 0.296f, -0.021f, 0.654f, 0.14f, 0.008f, 0.94f, -0.963f, 0.333f, -0.481f, -0.917f };
   int ldb = 3;
   float B_expected[] = { -0.312717f, 0.0986958f, 0.0456624f, 0.163957f, -0.2102f, 0.0234f, 0.143952f, 0.0170999f, 0.276937f, -0.480541f, 0.236f, 0.227f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1935) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1935) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3f, 0.1f};
   float A[] = { -0.932f, 0.532f, -0.763f, -0.029f, -0.524f, -0.938f, 0.007f, -0.445f, -0.659f, 0.709f, -0.581f, 0.825f, -0.904f, -0.453f, 0.119f, 0.964f, -0.649f, 0.48f };
   int lda = 3;
   float B[] = { -0.571f, 0.138f, 0.038f, -0.175f, 0.737f, 0.567f, -0.569f, 0.062f, 0.522f, -0.625f, 0.156f, 0.799f };
   int ldb = 3;
   float B_expected[] = { -0.0819591f, 0.15247f, -0.121808f, -0.00810757f, 0.287388f, -0.154159f, -0.0982488f, 0.13709f, -0.190946f, -0.223188f, 0.0729118f, 0.274542f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1936) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1936) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3f, 0.1f};
   float A[] = { 0.527f, 0.272f, 0.431f, 0.642f, -0.239f, -0.254f, -0.231f, 0.766f, 0.85f, -0.09f, 0.679f, -0.898f, 0.192f, -0.651f, -0.869f, 0.859f, 0.68f, 0.03f };
   int lda = 3;
   float B[] = { 0.867f, 0.816f, -0.643f, 0.509f, -0.594f, -0.833f, -0.174f, 0.51f, 0.676f, 0.115f, 0.261f, -0.409f };
   int ldb = 3;
   float B_expected[] = { -0.3417f, -0.1581f, 0.184172f, -0.515263f, 0.82684f, 0.153742f, 0.0012f, -0.1704f, -0.0834964f, -0.0053432f, -0.216529f, 0.104369f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1937) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1937) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3f, 0.1f};
   float A[] = { 0.606f, -0.757f, 0.503f, -0.649f, -0.269f, -0.484f, 0.626f, -0.107f, -0.867f, -0.047f, -0.779f, 0.675f, 0.249f, 0.645f, -0.755f, 0.242f, 0.941f, 0.189f };
   int lda = 3;
   float B[] = { -0.402f, 0.252f, -0.214f, 0.745f, 0.342f, -0.98f, -0.096f, 0.38f, -0.543f, 0.605f, 0.63f, -0.059f };
   int ldb = 2;
   float B_expected[] = { 0.349049f, -0.0955741f, -0.472341f, -0.259287f, -0.176304f, -0.239347f, 0.191174f, 0.170679f, 0.152979f, -0.219859f, -0.203592f, 0.0448683f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1938) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1938) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3f, 0.1f};
   float A[] = { -0.797f, -0.288f, 0.943f, -0.821f, -0.565f, 0.73f, -0.146f, -0.967f, 0.473f, -0.095f, 0.877f, 0.178f, -0.159f, 0.021f, -0.988f, 0.296f, 0.279f, -0.513f };
   int lda = 3;
   float B[] = { -0.455f, 0.859f, -0.21f, 0.702f, -0.591f, -0.235f, 0.519f, 0.279f, -0.444f, 0.816f, -0.507f, 0.893f };
   int ldb = 2;
   float B_expected[] = { -0.136371f, -0.712172f, -0.311667f, -0.302476f, 0.337384f, -0.259056f, -0.027248f, -0.327988f, 0.0516f, -0.2892f, 0.0628f, -0.3186f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1939) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1939) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3f, 0.1f};
   float A[] = { 0.791f, 0.19f, -0.549f, 0.994f, -0.822f, 0.679f, -0.586f, 0.042f, -0.159f, -0.86f, 0.065f, 0.943f, -0.545f, 0.403f, 0.199f, 0.76f, 0.159f, 0.715f };
   int lda = 3;
   float B[] = { -0.336f, 0.317f, 0.502f, 0.543f, 0.027f, 0.802f, 0.391f, 0.716f, -0.154f, 0.436f, 0.738f, -0.029f };
   int ldb = 2;
   float B_expected[] = { 0.119543f, -0.133991f, -0.212552f, -0.193533f, -0.239565f, -0.0842153f, -0.531028f, 0.229828f, 0.61223f, 0.265016f, 0.850081f, -0.810046f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1940) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1940) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   float alpha[2] = {-0.3f, 0.1f};
   float A[] = { -0.182f, -0.821f, -0.756f, -0.479f, -0.191f, -0.989f, -0.466f, 0.018f, 0.85f, 0.516f, -0.826f, 0.209f, -0.321f, -0.988f, -0.936f, -0.745f, -0.57f, -0.362f };
   int lda = 3;
   float B[] = { -0.501f, 0.915f, -0.928f, 0.722f, -0.542f, -0.828f, -0.875f, -0.981f, 0.425f, 0.347f, -0.929f, -0.596f };
   int ldb = 2;
   float B_expected[] = { 0.0588f, -0.3246f, 0.2062f, -0.3094f, 0.134369f, -0.0793628f, 0.368285f, -0.125876f, -0.344423f, -0.219222f, 0.402199f, -0.204129f };
   cblas_ctrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrsm(case 1941) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrsm(case 1941) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { 0.976, -0.41, -0.313, -0.779, -0.164, 0.571, 0.056, -0.526 };
   int lda = 2;
   double B[] = { -0.177, 0.837, 0.391, -0.853, -0.633, 0.693, -0.392, -0.356, -0.708, 0.926, -0.093, -0.337 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1942) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1942) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { -0.383, 0.141, 0.889, -0.007, -0.148, -0.068, 0.481, 0.675 };
   int lda = 2;
   double B[] = { 0.469, 0.735, -0.47, -0.164, 0.994, -0.483, -0.354, 0.357, 0.51, 0.523, 0.934, -0.592 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1943) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1943) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { -0.089, -0.391, -0.317, -0.349, 0.618, -0.541, -0.84, 0.31 };
   int lda = 2;
   double B[] = { 0.931, -0.257, -0.048, 0.633, -0.32, -0.576, -0.682, 0.953, -0.412, 0.408, -0.809, 0.092 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1944) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1944) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { 0.599, -0.01, -0.045, 0.567, 0.827, -0.969, -0.729, 0.538 };
   int lda = 2;
   double B[] = { 0.971, -0.626, -0.77, -0.882, 0.434, 0.269, -0.456, 0.497, 0.289, 0.957, 0.447, -0.921 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1945) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1945) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { 0.441, -0.501, 0.607, -0.4, -0.976, -0.523, -0.136, -0.492 };
   int lda = 2;
   double B[] = { 0.639, 0.872, -0.436, 0.518, 0.164, -0.04, 0.489, 0.201, 0.723, -0.958, 0.934, -0.549 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1946) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1946) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { -0.603, -0.475, 0.598, -0.666, -0.733, 0.04, 0.491, -0.592 };
   int lda = 2;
   double B[] = { 0.71, -0.827, 0.947, -0.364, 0.235, 0.294, 0.298, -0.401, -0.193, -0.008, 0.122, -0.47 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1947) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1947) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { -0.73, -0.823, 0.636, -0.965, 0.886, -0.236, 0.501, -0.301 };
   int lda = 2;
   double B[] = { 0.259, 0.701, -0.033, 0.616, -0.646, -0.177, -0.886, 0.589, -0.736, -0.303, -0.995, 0.982 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1948) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1948) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { 0.829, -0.889, 0.382, 0.083, 0.006, -0.76, -0.338, -0.601 };
   int lda = 2;
   double B[] = { 0.006, 0.381, 0.241, 0.096, -0.672, 0.664, 0.952, -0.376, -0.803, 0.344, -0.09, -0.175 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1949) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1949) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { 0.879, -0.511, -0.814, -0.94, 0.91, 0.761, 0.223, 0.03, -0.689, -0.739, -0.814, 0.463, 0.389, 0.615, -0.175, 0.129, -0.904, 0.102 };
   int lda = 3;
   double B[] = { 0.383, 0.328, 0.589, -0.29, 0.912, 0.327, 0.629, 0.883, -0.578, -0.708, 0.168, -0.982 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1950) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1950) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { -0.91, -0.182, 0.333, 0.193, 0.14, 0.538, 0.161, -0.034, -0.614, -0.154, 0.881, 0.842, 0.183, -0.229, 0.099, 0.062, -0.121, 0.179 };
   int lda = 3;
   double B[] = { -0.138, 0.109, -0.87, -0.161, 0.917, 0.443, 0.798, 0.677, -0.574, 0.327, -0.626, 0.446 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1951) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1951) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { -0.491, -0.021, -0.833, 0.921, -0.71, 0.282, 0.638, 0.223, -0.434, 0.921, -0.949, 0.457, -0.665, -0.844, -0.633, -0.874, -0.73, 0.637 };
   int lda = 3;
   double B[] = { -0.047, 0.714, 0.678, 0.756, 0.003, 0.359, 0.507, -0.197, -0.726, 0.873, -0.118, -0.996 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1952) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1952) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { 0.372, 0.354, -0.537, 0.948, -0.348, 0.808, 0.573, -0.797, 0.818, 0.701, -0.749, -0.801, -0.959, -0.781, 0.727, -0.189, 0.244, 0.414 };
   int lda = 3;
   double B[] = { 0.852, -0.714, 0.455, 0.171, -0.128, 0.554, 0.342, -0.203, 0.669, 0.619, -0.76, 0.759 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1953) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1953) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { 0.1, -0.975, 0.885, -0.608, -0.303, 0.87, -0.763, 0.409, 0.501, 0.522, -0.176, 0.679, -0.681, -0.815, -0.878, 0.86, 0.348, -0.65 };
   int lda = 3;
   double B[] = { -0.245, 0.954, -0.465, -0.931, 0.327, 0.288, -0.067, 0.252, 0.124, -0.073, -0.731, 0.176 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1954) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1954) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { 0.572, 0.045, -0.465, 0.113, 0.996, -0.597, 0.712, 0.945, 0.053, -0.436, 0.36, 0.035, -0.489, -0.012, 0.23, 0.22, 0.068, -0.586 };
   int lda = 3;
   double B[] = { -0.543, -0.809, -0.641, -0.744, 0.507, -0.742, -0.279, -0.835, -0.097, -0.968, 0.984, -0.813 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1955) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1955) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { 0.837, 0.576, -0.396, 0.013, -0.567, 0.59, 0.513, 0.824, 0.045, 0.486, 0.386, 0.766, 0.222, 0.042, 0.091, -0.008, 0.43, 0.102 };
   int lda = 3;
   double B[] = { 0.16, -0.958, -0.125, 0.833, 0.344, 0.213, 0.2, -0.689, 0.81, 0.415, -0.198, 0.001 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1956) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1956) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 111;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {0, 0};
   double A[] = { -0.351, 0.7, -0.495, 0.448, -0.229, 0.925, -0.269, 0.251, -0.783, -0.223, 0.582, 0.373, -0.095, -0.383, -0.087, -0.043, -0.315, -0.999 };
   int lda = 3;
   double B[] = { -0.067, -0.104, 0.92, -0.333, 0.367, 0.995, 0.86, 0.425, 0.12, -0.756, 0.441, -0.214 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1957) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1957) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.572, -0.073, 0.878, -0.688, -0.615, -0.213, -0.643, 0.809 };
   int lda = 2;
   double B[] = { -0.973, -0.481, 0.071, -0.71, -0.669, 0.717, -0.09, -0.304, -0.427, 0.625, 0.539, -0.565 };
   int ldb = 3;
   double B_expected[] = { 0.574560994608, 0.155494672389, 0.0371747871512, 0.389534544514, 0.283820482207, -0.45678514825, 0.591891359193, 0.214411302729, -0.27258111691, 0.507180331171, 0.645135319443, -0.46315922005 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1958) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1958) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.022, 0.475, 0.444, 0.252, -0.871, 0.867, -0.093, 0.264 };
   int lda = 2;
   double B[] = { 0.696, 0.259, 0.494, 0.162, -0.9, 0.143, 0.436, 0.487, -0.733, 0.138, -0.618, 0.572 };
   int ldb = 3;
   double B_expected[] = { -0.2347, -0.0081, -0.1644, 0.0008, 0.2557, -0.1329, -0.0773344, -0.0397592, 0.2792952, -0.0736264, -0.0188216, -0.2388288 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1959) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1959) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.918, -0.459, 0.547, 0.887, 0.4, -0.497, 0.49, -0.313 };
   int lda = 2;
   double B[] = { 0.028, 0.482, -0.59, -0.533, -0.594, 0.544, -0.717, -0.524, 0.07, -0.839, 0.538, -0.548 };
   int ldb = 3;
   double B_expected[] = { -0.258092239243, -0.278373561582, 0.128448307703, -0.0949352940165, 0.35005709854, -0.355276452021, 0.308556833073, 0.371588344391, -0.148348709879, 0.433197660833, -0.356526626221, 0.217565644883 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1960) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1960) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.86, -0.532, -0.396, -0.116, 0.766, -0.818, -0.335, 0.271 };
   int lda = 2;
   double B[] = { -0.029, -0.754, -0.566, -0.108, 0.904, -0.038, 0.07, -0.476, -0.48, 0.961, 0.864, -0.593 };
   int ldb = 3;
   double B_expected[] = { -0.058812, 0.130312, 0.419002, 0.272588, -0.330474, -0.264172, 0.0266, 0.1498, 0.0479, -0.3363, -0.1999, 0.2643 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1961) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1961) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.043, 0.25, -0.831, 0.609, -0.896, 0.886, 0.653, 0.065 };
   int lda = 2;
   double B[] = { 0.548, 0.076, 0.429, 0.873, -0.559, -0.329, -0.326, -0.174, 0.633, 0.489, 0.317, -0.896 };
   int ldb = 2;
   double B_expected[] = { 0.239257797324, 0.64684765886, 0.889006221152, 0.139062311692, 0.0322336011438, -0.807944179397, -0.977615509726, -1.02501063893, -0.164440783851, 0.983483814822, 1.28991055447, 1.90436729944 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1962) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1962) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.268, -0.062, -0.017, 0.326, 0.561, -0.203, -0.665, 0.338 };
   int lda = 2;
   double B[] = { -0.46, 0.954, 0.823, 0.945, -0.825, 0.882, -0.214, -0.095, -0.935, -0.245, 0.902, 0.904 };
   int ldb = 2;
   double B_expected[] = { 0.0426, -0.3322, -0.297862, -0.006188, 0.1593, -0.3471, 0.054794, 0.234161, 0.305, -0.02, -0.528045, -0.107865 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1963) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1963) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.075, 0.178, -0.321, -0.056, -0.124, -0.483, 0.685, -0.052 };
   int lda = 2;
   double B[] = { -0.47, -0.363, 0.766, -0.961, -0.391, -0.691, 0.42, -0.339, 0.45, -0.975, 0.991, -0.198 };
   int ldb = 2;
   double B_expected[] = { 0.874038948948, -0.779868445448, -0.234271045009, 0.514916650598, 0.810533012472, -1.05664738101, -0.149515922946, 0.198430908039, 2.17245126703, 0.115946317124, -0.420252834642, 0.199484456348 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1964) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1964) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.206, -0.461, -0.681, 0.358, 0.21, -0.318, 0.082, -0.097 };
   int lda = 2;
   double B[] = { 0.576, -0.249, 0.718, 0.424, 0.728, -0.464, 0.774, 0.541, -0.112, 0.803, 0.275, -0.638 };
   int ldb = 2;
   double B_expected[] = { -0.343295, 0.186865, -0.2578, -0.0554, -0.3973645, 0.2566785, -0.2863, -0.0849, 0.0189315, -0.0963345, -0.0187, 0.2189 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1965) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1965) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.117, 0.983, -0.929, -0.69, -0.144, 0.28, 0.658, 0.304, -0.657, 0.543, -0.051, -0.98, -0.846, -0.484, 0.052, 0.691, 0.613, -0.178 };
   int lda = 3;
   double B[] = { -0.688, 0.453, -0.63, 0.067, 0.193, 0.359, -0.792, 0.307, -0.501, -0.616, -0.595, 0.817 };
   int ldb = 3;
   double B_expected[] = { -0.566587593051, 0.340892661842, -0.458137993587, -0.0857620879204, -0.102500656517, -0.173972458173, -1.32599192297, -0.284341349955, -0.284178293736, -0.823318590512, 0.278700120014, -0.415972885216 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1966) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1966) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.773, -0.614, 0.782, -0.728, -0.727, -0.715, 0.858, -0.065, 0.922, 0.178, 0.588, 0.215, -0.92, -0.443, -0.583, -0.244, 0.996, -0.539 };
   int lda = 3;
   double B[] = { 0.159, 0.669, -0.692, 0.808, -0.146, 0.489, -0.385, -0.646, 0.704, -0.968, 0.551, -0.281 };
   int ldb = 3;
   double B_expected[] = { 0.0796383322, -0.0678193334, 0.0951193, -0.2156591, -0.0051, -0.1613, -0.2408434996, -0.0853028168, -0.0037554, 0.3083308, -0.1372, 0.1394 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1967) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1967) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.864, -0.382, -0.715, 0.227, -0.973, -0.709, -0.247, -0.601, 0.467, -0.133, 0.988, 0.937, -0.272, -0.334, 0.719, 0.992, 0.203, -0.646 };
   int lda = 3;
   double B[] = { 0.285, -0.409, -0.347, -0.925, -0.616, 0.422, 0.631, -0.954, -0.053, -0.255, -0.749, -0.979 };
   int ldb = 3;
   double B_expected[] = { -0.0215414266825, -0.165475896999, 0.469240391843, 0.538308411392, 1.71185240759, 0.063655952267, -0.0586080545035, -0.378370049976, 0.536158413721, 0.02961076215, 0.67769157898, -0.0939027988826 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1968) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1968) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.373, -0.335, -0.106, 0.542, -0.504, 0.574, -0.318, 0.043, -0.801, -0.331, 0.699, 0.776, -0.56, 0.131, 0.742, -0.692, -0.614, -0.874 };
   int lda = 3;
   double B[] = { -0.823, 0.929, -0.55, 0.172, -0.44, 0.067, 0.99, -0.013, 0.513, -0.438, -0.591, -0.302 };
   int ldb = 3;
   double B_expected[] = { 0.154, -0.361, 0.181249, -0.22802, 0.187552082, 0.008181148, -0.2957, 0.1029, -0.1997079, 0.2281373, 0.0457001502, -0.1796150434 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1969) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1969) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.912, 0.523, 0.314, -0.205, -0.895, 0.033, 0.157, -0.936, -0.582, 0.104, -0.868, 0.851, -0.131, 0.836, 0.993, 0.319, -0.684, -0.035 };
   int lda = 3;
   double B[] = { 0.07, -0.556, 0.018, -0.245, -0.405, 0.77, 0.888, 0.01, -0.81, -0.42, 0.66, -0.387 };
   int ldb = 2;
   double B_expected[] = { -0.132542904863, 0.151203976135, 0.45996395874, -0.700981460432, -0.771115355304, 0.0234040392321, 1.04091400336, -0.314874142966, -0.418936175202, -0.0443526810935, 0.218699329114, -0.27741882532 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1970) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1970) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.904, 0.983, -0.777, 0.503, 0.061, -0.442, 0.797, 0.415, -0.49, -0.466, 0.386, -0.147, -0.793, -0.381, -0.481, 0.33, 0.69, 0.35 };
   int lda = 3;
   double B[] = { 0.152, 0.832, 0.687, -0.287, -0.571, -0.187, -0.456, 0.631, 0.976, 0.833, -0.527, -0.188 };
   int ldb = 2;
   double B_expected[] = { -0.3155234788, -0.5211211034, -0.2870272698, 0.3910522396, -0.0411631, 0.0498567, 0.1600099, -0.2914973, -0.3761, -0.1523, 0.1769, 0.0037 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1971) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1971) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.527, 0.434, 0.025, 0.505, 0.724, 0.961, -0.071, 0.675, -0.334, 0.259, 0.167, 0.898, 0.116, 0.723, 0.086, 0.042, -0.483, -0.862 };
   int lda = 3;
   double B[] = { -0.874, 0.252, 0.924, 0.251, 0.559, -0.619, -0.131, -0.286, 0.09, -0.111, 0.062, -0.973 };
   int ldb = 2;
   double B_expected[] = { 0.116195543731, -0.404988360492, -0.325886265381, 0.300824742268, 0.86553022636, 0.0931927221532, -0.0931167995431, -0.760087414797, 0.774460770553, -0.204189465459, -0.501996021978, -0.354684266966 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1972) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1972) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 112;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.383, -0.184, 0.14, 0.131, -0.494, -0.025, -0.396, -0.183, 0.519, 0.806, -0.737, 0.764, -0.03, 0.622, -0.826, 0.605, 0.638, 0.935 };
   int lda = 3;
   double B[] = { 0.975, -0.816, -0.996, -0.038, -0.316, -0.31, -0.003, -0.974, 0.364, -0.217, 0.909, -0.656 };
   int ldb = 2;
   double B_expected[] = { -0.2109, 0.3423, 0.3026, -0.0882, 0.2001673, 0.0411059, 0.0443818, 0.2646074, -0.0213138923, 0.1426909311, 0.1794588402, 0.4128021586 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1973) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1973) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.433, -0.405, -0.008, 0.13, 0.377, -0.664, 0.421, -0.779 };
   int lda = 2;
   double B[] = { 0.022, -0.326, -0.905, 0.323, -0.722, 0.282, -0.877, -0.793, -0.906, -0.999, -0.607, -0.979 };
   int ldb = 3;
   double B_expected[] = { 0.0831887207906, -0.153137570623, -0.510564586332, -0.0447544052299, -0.412732352054, -0.0239182507667, 0.35364638809, -0.274824473121, 0.341954849059, -0.294570686181, 0.328230337479, -0.181800438645 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1974) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1974) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.007, 0.289, 0.434, 0.931, 0.776, -0.861, 0.83, -0.753 };
   int lda = 2;
   double B[] = { 0.775, -0.299, -0.45, 0.923, 0.251, 0.934, 0.388, -0.958, -0.732, 0.263, -0.5, 0.097 };
   int ldb = 3;
   double B_expected[] = { -0.2026, 0.1672, 0.0427, -0.3219, -0.1687, -0.2551, -0.0883348, 0.0650146, 0.4744571, 0.0273583, 0.4510139, -0.1254463 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1975) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1975) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.768, -0.738, 0.236, 0.721, 0.691, -0.963, -0.36, -0.376 };
   int lda = 2;
   double B[] = { -0.822, 0.174, 0.799, 0.8, -0.985, -0.169, 0.652, -0.529, -0.51, -0.506, -0.542, -0.786 };
   int ldb = 3;
   double B_expected[] = { -0.212429545832, 0.508667487335, 0.591670151369, 0.238559438419, 0.40264717438, -0.154881488703, 0.500259801606, -0.0994508738781, -0.130621162022, -0.416426547, -0.0684577231932, -0.575944733113 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1976) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1976) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 141;
   int uplo = 122;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.294, 0.843, 0.52, 0.53, 0.392, 0.293, 0.209, 0.497 };
   int lda = 2;
   double B[] = { 0.765, -0.547, 0.451, -0.581, 0.166, 0.834, -0.541, 0.278, -0.832, 0.66, -0.718, -0.664 };
   int ldb = 3;
   double B_expected[] = { -0.1872365, 0.3339085, -0.0667796, 0.3834252, -0.2809938, -0.2009734, 0.1345, -0.1375, 0.1836, -0.2812, 0.2818, 0.1274 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1977) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1977) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.642, 0.513, 0.428, 0.273, -0.612, 0.531, -0.664, 0.801 };
   int lda = 2;
   double B[] = { 0.429, -0.049, -0.661, 0.36, -0.247, 0.523, -0.227, 0.459, -0.902, 0.328, 0.37, -0.225 };
   int ldb = 2;
   double B_expected[] = { -0.161443909893, -0.0392846195877, 0.158306491417, 0.236544282705, 0.158671944063, -0.1560767799, 0.00300493937503, 0.254905467713, 0.369328020399, 0.00134777953987, -0.306971508873, -0.0836654236493 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1978) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1978) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 121;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.229, -0.461, -0.279, 0.674, -0.797, -0.286, 0.397, 0.329 };
   int lda = 2;
   double B[] = { 0.402, 0.728, 0.824, -0.691, -0.362, 0.437, 0.192, 0.788, -0.259, 0.599, 0.79, 0.076 };
   int ldb = 2;
   double B_expected[] = { -0.1934, -0.1782, -0.383205, 0.202987, 0.0649, -0.1673, -0.1325225, -0.3690995, 0.0178, -0.2056, -0.289215, -0.112754 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1979) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1979) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.942, -0.832, 0.595, -0.092, 0.01, 0.001, 0.944, 0.256 };
   int lda = 2;
   double B[] = { 0.73, 0.488, -0.363, -0.01, -0.112, 0.169, -0.268, -0.13, -0.657, 0.573, 0.91, 0.632 };
   int ldb = 2;
   double B_expected[] = { 0.158268746344, 0.226988691038, 0.117355164571, -0.00345029435376, -0.0289643723553, 0.0722018494696, 0.0888981803586, 0.0370317099277, -0.233113998714, -0.101765761072, -0.305361921327, -0.187259165106 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1980) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1980) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 141;
   int uplo = 122;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.508, 0.053, -0.516, 0.785, -0.451, -0.53, 0.551, 0.235 };
   int lda = 2;
   double B[] = { -0.09, 0.46, 0.948, 0.918, -0.337, 0.012, -0.786, -0.676, 0.906, -0.38, -0.566, 0.645 };
   int ldb = 2;
   double B_expected[] = { -0.0713482, -0.5355066, -0.3762, -0.1806, 0.1589574, 0.2649562, 0.3034, 0.1242, 0.0168633, 0.1582089, 0.1053, -0.2501 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1981) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1981) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.499, -0.268, 0.234, 0.032, -0.158, 0.684, -0.878, 0.613, 0.968, 0.812, 0.013, 0.34, -0.485, -0.565, 0.316, 0.286, -0.459, 0.637 };
   int lda = 3;
   double B[] = { -0.964, 0.804, 0.197, 0.141, 0.942, 0.474, 0.741, -0.441, -0.738, -0.703, -0.27, 0.98 };
   int ldb = 3;
   double B_expected[] = { 0.561582612433, -0.70128258354, -0.0253749021391, 0.0631927226609, 0.295313488523, -0.305260767297, -0.0937671252683, 0.884164549696, 0.000683977216651, 0.260184505619, 0.344358828778, 0.221445372699 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1982) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1982) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 121;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.453, 0.917, 0.131, 0.361, 0.087, 0.441, -0.439, 0.439, 0.777, 0.131, 0.535, 0.646, 0.508, 0.746, -0.347, -0.911, -0.874, -0.525 };
   int lda = 3;
   double B[] = { -0.739, -0.776, -0.049, 0.548, -0.39, -0.856, -0.757, 0.307, -0.533, -0.342, 0.431, 0.618 };
   int ldb = 3;
   double B_expected[] = { 0.2794424312, 0.1451980676, -0.2891898, -0.1549434, 0.2026, 0.2178, 0.2242026328, -0.0997909546, 0.3882643, 0.0019799, -0.1911, -0.1423 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1983) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1983) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.835, -0.775, -0.384, -0.128, -0.41, -0.511, -0.282, -0.341, -0.856, -0.662, 0.721, -0.939, 0.175, -0.899, 0.832, -0.519, 0.652, -0.318 };
   int lda = 3;
   double B[] = { -0.654, 0.105, -0.39, 0.645, 0.867, 0.045, -0.842, -0.896, -0.249, 0.419, 0.575, 0.561 };
   int ldb = 3;
   double B_expected[] = { -0.177337134492, -0.0485464421929, -0.0947130836909, 0.143712701441, -0.0502556531648, 0.286334558029, -0.109929498786, -0.323108217437, -0.0362323282558, 0.21056630482, -0.514117706819, 0.0792536824901 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1984) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1984) imag");
     };
   };
  };


  {
   int order = 101;
   int side = 142;
   int uplo = 122;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.136, 0.272, 0.676, 0.673, -0.659, 0.668, 0.991, -0.569, -0.489, 0.581, -0.232, -0.249, -0.396, -0.832, 0.763, -0.092, 0.117, 0.108 };
   int lda = 3;
   double B[] = { 0.721, -0.141, -0.604, 0.318, 0.387, 0.73, -0.549, 0.302, 0.101, 0.721, -0.064, 0.673 };
   int ldb = 3;
   double B_expected[] = { -0.2022, 0.1144, 0.4148738, -0.1541186, -0.5047180206, 0.1126569022, 0.1345, -0.1455, -0.318479, -0.13854, 0.114359797, -0.242815912 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1985) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1985) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.578, 0.601, -0.43, -0.187, -0.934, 0.635, 0.157, -0.561, -0.964, 0.025, 0.435, -0.674, -0.575, 0.275, 0.609, 0.228, -0.202, -0.267 };
   int lda = 3;
   double B[] = { 0.505, -0.347, 0.213, -0.392, -0.465, -0.918, -0.737, -0.974, -0.051, 0.97, 0.066, 0.604 };
   int ldb = 2;
   double B_expected[] = { -0.206165616299, -0.811510964363, -0.328765954464, -0.593889594613, -0.410790365608, 0.365230809488, -0.377900693873, 0.166778025696, -0.558066070138, 0.728199798382, -0.271362172482, 0.505674752215 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1986) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1986) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 121;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.311, -0.737, -0.738, -0.214, -0.387, -0.043, -0.168, 0.563, 0.165, 0.007, -0.121, 0.408, -0.75, -0.641, -0.997, -0.347, 0.523, -0.922 };
   int lda = 3;
   double B[] = { 0.46, 0.376, -0.623, -0.092, 0.233, 0.981, -0.435, -0.493, 0.405, 0.855, -0.391, 0.572 };
   int ldb = 2;
   double B_expected[] = { -0.311417159, -0.418726217, 0.2053384662, -0.1587052684, -0.449331, -0.414523, 0.1666068, -0.1265226, -0.207, -0.216, 0.0601, -0.2107 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1987) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1987) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 113;
   int diag = 131;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { -0.288, -0.421, 0.451, 0.234, 0.67, -0.483, 0.273, 0.131, 0.005, 0.091, -0.706, -0.191, 0.285, -0.434, 0.648, -0.556, -0.886, 0.798 };
   int lda = 3;
   double B[] = { 0.359, -0.682, -0.618, 0.479, 0.463, 0.468, -0.43, 0.058, -0.361, -0.058, -0.028, -0.729 };
   int ldb = 2;
   double B_expected[] = { 0.432870841901, -0.202296442916, -0.484714722217, 0.00498299287046, -1.27917612947, -3.59551100448, 2.13407463306, 3.62604336509, 2.50059207751, 0.44116664838, -3.08374361183, -0.156015309482 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1988) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1988) imag");
     };
   };
  };


  {
   int order = 102;
   int side = 142;
   int uplo = 122;
   int trans = 113;
   int diag = 132;
   int M = 2;
   int N = 3;
   double alpha[2] = {-0.3, 0.1};
   double A[] = { 0.288, 0.961, -0.571, -0.341, -0.443, 0.116, -0.928, 0.157, 0.035, 0.822, 0.733, -0.15, 0.851, -0.634, -0.769, -0.709, 0.346, -0.943 };
   int lda = 3;
   double B[] = { -0.708, 0.945, -0.144, 0.505, 0.827, -0.467, 0.883, 0.194, -0.607, -0.332, 0.716, -0.117 };
   int ldb = 2;
   double B_expected[] = { 0.1179, -0.3543, -0.0073, -0.1659, -0.2548954, -0.0197092, -0.3450402, -0.0621396, 0.4925104482, -0.0516973464, 0.0565040266, 0.1296638568 };
   cblas_ztrsm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrsm(case 1989) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrsm(case 1989) imag");
     };
   };
  };


}
