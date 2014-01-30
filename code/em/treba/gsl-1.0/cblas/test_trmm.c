#include <gsl/gsl_test.h>
#include <gsl/gsl_ieee_utils.h>
#include <gsl/gsl_math.h>
#include "gsl_cblas.h"

#include "tests.h"

void
test_trmm (void) {
const double flteps = 1e-4, dbleps = 1e-6;
  {
   int order = 101;
   int side = 141;
   int uplo = 121;
   int trans = 111;
   int diag = 131;
   int M = 2;
   int N = 3;
   float alpha = 1.0f;
   float A[] = { 0.565f, 0.967f, -0.969f, 0.184f };
   int lda = 2;
   float B[] = { 0.842f, -0.918f, -0.748f, -0.859f, -0.463f, 0.292f };
   int ldb = 3;
   float B_expected[] = { -0.354923f, -0.966391f, -0.140256f, -0.158056f, -0.085192f, 0.053728f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1670)");
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
   float alpha = 1.0f;
   float A[] = { -0.748f, 0.548f, 0.245f, 0.761f };
   int lda = 2;
   float B[] = { 0.349f, -0.552f, -0.682f, -0.71f, 0.475f, -0.59f };
   int ldb = 3;
   float B_expected[] = { -0.04008f, -0.2917f, -1.00532f, -0.71f, 0.475f, -0.59f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1671)");
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
   float alpha = 1.0f;
   float A[] = { 0.788f, 0.617f, -0.998f, -0.97f };
   int lda = 2;
   float B[] = { -0.4f, 0.773f, 0.074f, -0.388f, 0.825f, -0.608f };
   int ldb = 3;
   float B_expected[] = { -0.3152f, 0.609124f, 0.058312f, 0.77556f, -1.5717f, 0.515908f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1672)");
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
   float alpha = 1.0f;
   float A[] = { 0.01f, 0.387f, -0.953f, -0.374f };
   int lda = 2;
   float B[] = { 0.364f, 0.09f, 0.588f, -0.263f, 0.584f, 0.463f };
   int ldb = 3;
   float B_expected[] = { 0.364f, 0.09f, 0.588f, -0.609892f, 0.49823f, -0.097364f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1673)");
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
   float alpha = 1.0f;
   float A[] = { -0.586f, -0.426f, 0.765f, -0.239f };
   int lda = 2;
   float B[] = { -0.673f, -0.724f, 0.217f, -0.672f, -0.378f, -0.005f };
   int ldb = 2;
   float B_expected[] = { -0.159482f, 0.173036f, -0.641242f, 0.160608f, 0.217683f, 0.001195f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1674)");
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
   float alpha = 1.0f;
   float A[] = { -0.668f, 0.962f, 0.515f, 0.292f };
   int lda = 2;
   float B[] = { -0.145f, -0.337f, 0.718f, -0.866f, -0.454f, -0.439f };
   int ldb = 2;
   float B_expected[] = { -0.318555f, -0.337f, 0.27201f, -0.866f, -0.680085f, -0.439f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1675)");
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
   float alpha = 1.0f;
   float A[] = { -0.125f, -0.676f, 0.181f, 0.741f };
   int lda = 2;
   float B[] = { 0.354f, -0.366f, 0.455f, 0.134f, -0.564f, -0.303f };
   int ldb = 2;
   float B_expected[] = { -0.04425f, -0.51051f, -0.056875f, -0.208286f, 0.0705f, 0.156741f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1676)");
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
   float alpha = 1.0f;
   float A[] = { -0.162f, 0.542f, -0.839f, -0.935f };
   int lda = 2;
   float B[] = { 0.216f, 0.766f, -0.228f, -0.097f, 0.205f, 0.875f };
   int ldb = 2;
   float B_expected[] = { 0.216f, 0.883072f, -0.228f, -0.220576f, 0.205f, 0.98611f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1677)");
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
   float alpha = 1.0f;
   float A[] = { -0.353f, -0.854f, -0.502f, 0.591f, -0.934f, -0.729f, 0.063f, 0.352f, 0.126f };
   int lda = 3;
   float B[] = { 0.2f, -0.626f, -0.694f, -0.889f, -0.251f, -0.42f };
   int ldb = 3;
   float B_expected[] = { -0.0706f, 0.413884f, 0.26851f, 0.313817f, 0.99364f, 0.576337f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1678)");
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
   float alpha = 1.0f;
   float A[] = { -0.864f, -0.046f, -0.755f, 0.12f, 0.525f, 0.917f, 0.571f, -0.098f, -0.226f };
   int lda = 3;
   float B[] = { -0.905f, -0.296f, -0.927f, -0.813f, 0.624f, -0.366f };
   int ldb = 3;
   float B_expected[] = { -0.905f, -0.25437f, -0.515157f, -0.813f, 0.661398f, 0.820023f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1679)");
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
   float alpha = 1.0f;
   float A[] = { -0.69f, -0.927f, -0.281f, -0.918f, -0.527f, -0.652f, -0.393f, -0.954f, 0.651f };
   int lda = 3;
   float B[] = { -0.587f, 0.788f, -0.629f, -0.444f, 0.515f, 0.081f };
   int ldb = 3;
   float B_expected[] = { -0.071157f, 0.18479f, -0.409479f, -0.198243f, -0.348679f, 0.052731f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1680)");
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
   float alpha = 1.0f;
   float A[] = { -0.082f, -0.077f, 0.811f, 0.852f, 0.224f, 0.443f, -0.509f, 0.171f, 0.986f };
   int lda = 3;
   float B[] = { -0.982f, 0.388f, -0.493f, -0.497f, -0.605f, 0.433f };
   int ldb = 3;
   float B_expected[] = { -0.400487f, 0.303697f, -0.493f, -1.23286f, -0.530957f, 0.433f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1681)");
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
   float alpha = 1.0f;
   float A[] = { 0.97f, -0.666f, 0.066f, -0.176f, 0.402f, 0.286f, -0.703f, 0.962f, 0.912f };
   int lda = 3;
   float B[] = { -0.644f, -0.97f, 0.814f, -0.777f, 0.812f, 0.254f };
   int ldb = 2;
   float B_expected[] = { -0.62468f, -0.9409f, 0.440572f, -0.141634f, 1.97634f, 0.166084f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1682)");
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
   float alpha = 1.0f;
   float A[] = { 0.714f, 0.468f, 0.859f, -0.547f, 0.076f, 0.542f, 0.512f, -0.987f, -0.167f };
   int lda = 3;
   float B[] = { -0.238f, -0.336f, 0.402f, 0.945f, -0.242f, -0.062f };
   int ldb = 2;
   float B_expected[] = { -0.238f, -0.336f, 0.532186f, 1.12879f, -0.76063f, -1.16675f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1683)");
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
   float alpha = 1.0f;
   float A[] = { -0.723f, 0.041f, 0.333f, -0.682f, 0.193f, 0.581f, 0.963f, -0.757f, 0.396f };
   int lda = 3;
   float B[] = { 0.047f, -0.701f, -0.25f, -0.779f, 0.435f, 0.612f };
   int ldb = 2;
   float B_expected[] = { 0.100624f, 0.67868f, 0.204485f, 0.205225f, 0.17226f, 0.242352f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1684)");
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
   float alpha = 1.0f;
   float A[] = { -0.13f, 0.511f, -0.544f, 0.938f, -0.126f, -0.873f, 0.118f, -0.75f, 0.674f };
   int lda = 3;
   float B[] = { -0.927f, -0.558f, -0.289f, -0.66f, 0.83f, 0.363f };
   int ldb = 2;
   float B_expected[] = { -1.5262f, -1.09273f, -1.01359f, -0.976899f, 0.83f, 0.363f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1685)");
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
   float alpha = 0.1f;
   float A[] = { -0.625f, -0.123f, -0.48f, -0.088f };
   int lda = 2;
   float B[] = { 0.376f, -0.46f, -0.813f, 0.419f, 0.792f, 0.226f };
   int ldb = 3;
   float B_expected[] = { -0.0235f, 0.02875f, 0.0508125f, -0.008312f, -0.0013116f, 0.0080111f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1686)");
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
   float alpha = 0.1f;
   float A[] = { -0.038f, -0.105f, -0.946f, 0.474f };
   int lda = 2;
   float B[] = { -0.757f, 0.974f, -0.045f, -0.809f, 0.654f, 0.611f };
   int ldb = 3;
   float B_expected[] = { -0.0757f, 0.0974f, -0.0045f, -0.0729515f, 0.055173f, 0.0615725f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1687)");
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
   float alpha = 0.1f;
   float A[] = { -0.328f, 0.713f, 0.781f, 0.084f };
   int lda = 2;
   float B[] = { -0.097f, 0.442f, -0.563f, 0.065f, -0.18f, 0.63f };
   int ldb = 3;
   float B_expected[] = { 0.0082581f, -0.0285556f, 0.0676694f, 5.46e-04f, -0.001512f, 0.005292f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1688)");
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
   float alpha = 0.1f;
   float A[] = { 0.261f, -0.659f, -0.536f, 0.694f };
   int lda = 2;
   float B[] = { -0.498f, 0.692f, 0.125f, 0.706f, -0.118f, -0.907f };
   int ldb = 3;
   float B_expected[] = { -0.0876416f, 0.0755248f, 0.0611152f, 0.0706f, -0.0118f, -0.0907f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1689)");
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
   float alpha = 0.1f;
   float A[] = { -0.669f, 0.416f, 0.761f, -0.359f };
   int lda = 2;
   float B[] = { -0.305f, -0.675f, -0.442f, 0.566f, 0.064f, 0.962f };
   int ldb = 2;
   float B_expected[] = { 0.0204045f, 0.001022f, 0.0295698f, -0.0539556f, -0.0042816f, -0.0296654f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1690)");
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
   float alpha = 0.1f;
   float A[] = { 0.565f, 0.386f, 0.643f, -0.028f };
   int lda = 2;
   float B[] = { 0.863f, -0.241f, 0.766f, 0.656f, -0.977f, 0.274f };
   int ldb = 2;
   float B_expected[] = { 0.0863f, 0.0313909f, 0.0766f, 0.114854f, -0.0977f, -0.0354211f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1691)");
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
   float alpha = 0.1f;
   float A[] = { 0.116f, 0.534f, 0.043f, 0.73f };
   int lda = 2;
   float B[] = { -0.758f, -0.63f, -0.043f, 0.666f, -0.088f, 0.382f };
   int ldb = 2;
   float B_expected[] = { -0.0424348f, -0.04599f, 0.0350656f, 0.048618f, 0.019378f, 0.027886f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1692)");
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
   float alpha = 0.1f;
   float A[] = { 0.48f, -0.63f, -0.786f, -0.437f };
   int lda = 2;
   float B[] = { 0.945f, 0.528f, -0.855f, -0.587f, 0.062f, 0.372f };
   int ldb = 2;
   float B_expected[] = { 0.061236f, 0.0528f, -0.048519f, -0.0587f, -0.017236f, 0.0372f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1693)");
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
   float alpha = 0.1f;
   float A[] = { -0.822f, -0.068f, 0.119f, -0.244f, -0.05f, 0.685f, 0.752f, -0.059f, -0.935f };
   int lda = 3;
   float B[] = { -0.431f, -0.753f, -0.319f, 0.164f, 0.979f, 0.885f };
   int ldb = 3;
   float B_expected[] = { 0.0367525f, -0.0180865f, 0.0298265f, -0.0096065f, 0.0557275f, -0.0827475f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1694)");
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
   float alpha = 0.1f;
   float A[] = { 0.97f, -0.408f, 0.174f, -0.308f, 0.997f, -0.484f, 0.322f, -0.183f, 0.849f };
   int lda = 3;
   float B[] = { -0.571f, 0.696f, -0.256f, -0.178f, 0.098f, 0.004f };
   int ldb = 3;
   float B_expected[] = { -0.0899512f, 0.0819904f, -0.0256f, -0.0217288f, 0.0096064f, 4.0e-04f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1695)");
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
   float alpha = 0.1f;
   float A[] = { -0.831f, 0.73f, 0.407f, 0.721f, 0.086f, -0.294f, 0.941f, -0.656f, -0.066f };
   int lda = 3;
   float B[] = { -0.051f, -0.343f, -0.98f, 0.722f, -0.372f, 0.466f };
   int ldb = 3;
   float B_expected[] = { 0.0042381f, -0.0066269f, 0.0241697f, -0.0599982f, 0.048857f, 0.0892678f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1696)");
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
   float alpha = 0.1f;
   float A[] = { 0.472f, 0.137f, -0.341f, 0.386f, -0.578f, 0.863f, -0.415f, -0.547f, -0.023f };
   int lda = 3;
   float B[] = { 0.582f, 0.141f, -0.306f, -0.047f, -0.162f, -0.784f };
   int ldb = 3;
   float B_expected[] = { 0.0582f, 0.0365652f, -0.0624657f, -0.0047f, -0.0180142f, -0.0675881f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1697)");
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
   float alpha = 0.1f;
   float A[] = { -0.775f, 0.762f, -0.038f, -0.8f, 0.626f, -0.701f, 0.639f, 0.239f, 0.34f };
   int lda = 3;
   float B[] = { 0.42f, 0.917f, 0.485f, 0.844f, -0.832f, 0.179f };
   int ldb = 2;
   float B_expected[] = { -0.124515f, -0.127149f, 0.0104762f, 0.0571125f, -0.028288f, 0.006086f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1698)");
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
   float alpha = 0.1f;
   float A[] = { -0.675f, 0.283f, 0.785f, -0.0f, -0.592f, -0.661f, 0.149f, -0.129f, 0.149f };
   int lda = 3;
   float B[] = { 0.964f, -0.575f, -0.215f, 0.953f, 0.527f, -0.418f };
   int ldb = 2;
   float B_expected[] = { 0.104252f, -0.0637282f, -0.0282983f, 0.100692f, 0.0527f, -0.0418f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1699)");
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
   float alpha = 0.1f;
   float A[] = { -0.225f, -0.943f, 0.839f, 0.759f, 0.752f, 0.807f, 0.288f, -0.276f, 0.434f };
   int lda = 3;
   float B[] = { -0.234f, 0.275f, 0.658f, -0.423f, -0.807f, -0.683f };
   int ldb = 2;
   float B_expected[] = { 0.005265f, -0.0061875f, 0.0715478f, -0.0577421f, -0.0015558f, -0.0407058f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1700)");
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
   float alpha = 0.1f;
   float A[] = { -0.043f, -0.983f, 0.479f, -0.136f, 0.048f, 0.745f, -0.408f, -0.731f, -0.953f };
   int lda = 3;
   float B[] = { 0.917f, 0.682f, -0.32f, 0.557f, -0.302f, 0.989f };
   int ldb = 2;
   float B_expected[] = { 0.0917f, 0.0682f, -0.122141f, -0.0113406f, -0.0101157f, 0.173064f };
   cblas_strmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], flteps, "strmm(case 1701)");
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
   double alpha = 0;
   double A[] = { -0.561, -0.114, -0.148, 0.488 };
   int lda = 2;
   double B[] = { 0.684, 0.38, 0.419, -0.361, 0.378, -0.423 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1702)");
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
   double alpha = 0;
   double A[] = { -0.378, 0.607, 0.41, 0.418 };
   int lda = 2;
   double B[] = { 0.146, -0.688, -0.953, -0.983, 0.237, 0.128 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1703)");
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
   double alpha = 0;
   double A[] = { -0.31, 0.277, -0.587, 0.885 };
   int lda = 2;
   double B[] = { -0.221, -0.831, -0.319, -0.547, -0.577, 0.295 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1704)");
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
   double alpha = 0;
   double A[] = { -0.577, 0.861, -0.439, -0.916 };
   int lda = 2;
   double B[] = { -0.933, -0.582, 0.528, 0.268, -0.804, 0.62 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1705)");
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
   double alpha = 0;
   double A[] = { -0.824, -0.119, -0.399, -0.653 };
   int lda = 2;
   double B[] = { 0.452, -0.168, 0.256, 0.554, 0.342, 0.318 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1706)");
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
   double alpha = 0;
   double A[] = { -0.299, 0.837, -0.03, 0.552 };
   int lda = 2;
   double B[] = { -0.83, -0.82, -0.362, -0.252, -0.062, -0.942 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1707)");
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
   double alpha = 0;
   double A[] = { -0.545, -0.107, 0.096, 0.183 };
   int lda = 2;
   double B[] = { -0.43, 0.841, 0.035, 0.7, 0.637, 0.095 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1708)");
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
   double alpha = 0;
   double A[] = { 0.626, 0.123, -0.959, 0.971 };
   int lda = 2;
   double B[] = { 0.185, -0.218, -0.074, 0.49, 0.802, -0.454 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1709)");
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
   double alpha = 0;
   double A[] = { -0.131, 0.048, 0.148, 0.834, -0.98, -0.009, -0.727, 0.241, 0.276 };
   int lda = 3;
   double B[] = { 0.75, -0.664, -0.136, -0.793, -0.742, 0.126 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1710)");
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
   double alpha = 0;
   double A[] = { 0.431, -0.387, 0.427, 0.495, 0.282, 0.158, -0.335, 0.535, -0.978 };
   int lda = 3;
   double B[] = { 0.518, -0.489, 0.899, -0.375, 0.376, -0.831 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1711)");
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
   double alpha = 0;
   double A[] = { -0.669, -0.976, -0.2, 0.661, -0.975, -0.965, -0.861, -0.779, -0.73 };
   int lda = 3;
   double B[] = { 0.31, 0.023, -0.853, 0.632, -0.174, 0.608 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1712)");
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
   double alpha = 0;
   double A[] = { 0.153, -0.408, -0.127, -0.634, -0.384, -0.815, 0.051, -0.096, 0.476 };
   int lda = 3;
   double B[] = { 0.343, -0.665, -0.348, 0.748, 0.893, 0.91 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1713)");
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
   double alpha = 0;
   double A[] = { -0.918, -0.19, 0.829, 0.942, 0.885, 0.087, 0.321, 0.67, -0.475 };
   int lda = 3;
   double B[] = { 0.377, 0.931, 0.291, -0.603, -0.617, 0.402 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1714)");
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
   double alpha = 0;
   double A[] = { -0.598, -0.232, -0.64, 0.595, 0.642, -0.921, -0.679, -0.846, -0.921 };
   int lda = 3;
   double B[] = { 0.032, -0.036, -0.278, -0.83, 0.922, -0.701 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1715)");
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
   double alpha = 0;
   double A[] = { 0.341, -0.858, -0.559, 0.499, -0.114, 0.57, 0.847, -0.612, 0.593 };
   int lda = 3;
   double B[] = { 0.672, 0.292, 0.752, 0.842, 0.625, 0.967 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1716)");
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
   double alpha = 0;
   double A[] = { 0.958, 0.823, -0.181, 0.141, 0.932, 0.097, -0.636, 0.844, 0.205 };
   int lda = 3;
   double B[] = { 0.113, -0.658, 0.703, -0.023, -0.384, 0.439 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1717)");
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
   double alpha = 0;
   double A[] = { 0.675, -0.468, -0.564, 0.71 };
   int lda = 2;
   double B[] = { -0.401, -0.823, 0.342, -0.384, 0.344, 0.18 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1718)");
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
   double alpha = 0;
   double A[] = { 0.932, -0.388, 0.432, -0.167 };
   int lda = 2;
   double B[] = { -0.624, 0.023, 0.065, 0.678, 0.044, -0.472 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1719)");
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
   double alpha = 0;
   double A[] = { -0.738, 0.649, -0.171, -0.462 };
   int lda = 2;
   double B[] = { -0.277, -0.519, -0.501, -0.024, -0.767, -0.591 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1720)");
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
   double alpha = 0;
   double A[] = { -0.17, -0.184, -0.243, 0.907 };
   int lda = 2;
   double B[] = { 0.593, 0.131, -0.317, -0.254, -0.948, 0.002 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1721)");
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
   double alpha = 0;
   double A[] = { 0.06, -0.838, -0.455, -0.715 };
   int lda = 2;
   double B[] = { -0.423, 0.665, -0.023, -0.872, -0.313, -0.698 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1722)");
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
   double alpha = 0;
   double A[] = { -0.506, 0.792, 0.338, -0.155 };
   int lda = 2;
   double B[] = { -0.257, -0.19, 0.201, 0.685, 0.663, 0.302 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1723)");
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
   double alpha = 0;
   double A[] = { 0.739, -0.996, 0.182, 0.626 };
   int lda = 2;
   double B[] = { 0.009, 0.485, -0.633, -0.08, -0.579, 0.223 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1724)");
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
   double alpha = 0;
   double A[] = { 0.777, 0.723, 0.378, 0.98 };
   int lda = 2;
   double B[] = { 0.291, -0.267, -0.076, 0.103, -0.021, -0.866 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1725)");
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
   double alpha = 0;
   double A[] = { -0.771, 0.469, 0.822, -0.619, 0.953, -0.706, 0.318, 0.559, -0.68 };
   int lda = 3;
   double B[] = { -0.32, 0.362, 0.719, -0.661, -0.504, 0.595 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1726)");
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
   double alpha = 0;
   double A[] = { 0.073, -0.501, -0.561, -0.229, -0.533, -0.138, 0.924, -0.164, -0.023 };
   int lda = 3;
   double B[] = { -0.208, 0.49, 0.827, 0.641, -0.884, -0.624 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1727)");
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
   double alpha = 0;
   double A[] = { 0.33, -0.649, -0.43, -0.266, 0.787, 0.449, 0.435, -0.774, -0.447 };
   int lda = 3;
   double B[] = { -0.687, -0.459, 0.189, 0.762, -0.039, 0.047 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1728)");
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
   double alpha = 0;
   double A[] = { 0.981, 0.242, 0.581, 0.064, 0.792, -0.529, 0.461, 0.224, -0.419 };
   int lda = 3;
   double B[] = { 0.285, 0.274, -0.912, 0.601, 0.24, 0.06 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1729)");
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
   double alpha = 0;
   double A[] = { -0.582, 0.269, -0.587, 0.68, -0.59, -0.936, 0.236, -0.728, -0.434 };
   int lda = 3;
   double B[] = { 0.113, 0.468, 0.943, 0.48, 0.215, -0.525 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1730)");
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
   double alpha = 0;
   double A[] = { -0.344, -0.938, 0.556, -0.678, -0.612, -0.519, -0.578, -0.848, 0.699 };
   int lda = 3;
   double B[] = { 0.915, -0.118, 0.538, -0.186, -0.413, -0.216 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1731)");
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
   double alpha = 0;
   double A[] = { -0.843, 0.54, -0.892, -0.296, 0.786, 0.136, 0.731, -0.418, -0.118 };
   int lda = 3;
   double B[] = { -0.775, 0.5, -0.399, -0.709, 0.779, 0.774 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1732)");
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
   double alpha = 0;
   double A[] = { -0.765, 0.233, 0.318, 0.547, -0.469, 0.023, -0.867, 0.687, -0.912 };
   int lda = 3;
   double B[] = { 0.019, -0.145, 0.472, 0.333, 0.527, -0.224 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_dtrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[i], B_expected[i], dbleps, "dtrmm(case 1733)");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { -0.852f, -0.409f, 0.871f, -0.854f, -0.493f, 0.444f, 0.973f, 0.027f };
   int lda = 2;
   float B[] = { -0.561f, 0.132f, 0.689f, 0.653f, -0.758f, -0.109f, -0.596f, 0.395f, -0.561f, 0.378f, 0.21f, 0.51f };
   int ldb = 3;
   float B_expected[] = { -0.0970014f, 0.0350174f, 0.0029825f, -0.048577f, -0.066776f, 0.121969f, -0.0368243f, -0.0590573f, -0.0352647f, -0.0556059f, -0.05019f, 0.019056f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1734) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1734) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { 0.349f, 0.0f, -0.462f, 0.91f, -0.693f, 0.587f, -0.617f, 0.112f };
   int lda = 2;
   float B[] = { 0.842f, -0.473f, 0.825f, 0.866f, 0.986f, 0.686f, 0.346f, 0.299f, -0.659f, 0.009f, 0.007f, -0.478f };
   int ldb = 3;
   float B_expected[] = { 0.0296278f, 0.0410058f, -0.0262152f, 0.112127f, -0.0913206f, 0.141775f, -0.0299f, 0.0346f, -0.0009f, -0.0659f, 0.0478f, 0.0007f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1735) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1735) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { 0.661f, -0.823f, 0.28f, 0.171f, 0.267f, 0.66f, 0.844f, 0.472f };
   int lda = 2;
   float B[] = { -0.256f, -0.518f, -0.933f, 0.066f, -0.513f, -0.286f, 0.109f, 0.372f, -0.183f, 0.482f, 0.362f, -0.436f };
   int ldb = 3;
   float B_expected[] = { 0.013171f, -0.059553f, -0.0811485f, -0.0562395f, -0.0233153f, -0.0574471f, -0.005815f, 0.018994f, 0.0277726f, -0.0674627f, 0.0612062f, 0.0563109f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1736) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1736) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { 0.623f, 0.314f, -0.594f, 0.717f, 0.566f, 0.001f, -0.411f, -0.387f };
   int lda = 2;
   float B[] = { -0.083f, 0.937f, -0.814f, 0.9f, -0.042f, 0.678f, -0.928f, 0.228f, 0.965f, -0.16f, 0.006f, -0.281f };
   int ldb = 3;
   float B_expected[] = { -0.0937f, -0.0083f, -0.09f, -0.0814f, -0.0678f, -0.0042f, -0.0758259f, -0.0975915f, -0.0348586f, 0.0503376f, -0.0102706f, -0.001845f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1737) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1737) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { -0.247f, -0.582f, 0.651f, -0.534f, -0.491f, 0.346f, 0.936f, -0.227f };
   int lda = 2;
   float B[] = { -0.002f, -0.02f, 0.162f, -0.62f, 0.632f, -0.07f, 0.352f, 0.042f, 0.574f, 0.272f, -0.139f, 0.012f };
   int ldb = 2;
   float B_expected[] = { -0.0366576f, 0.0123832f, 0.0617094f, 0.0010892f, 0.0249364f, -0.0384208f, 0.0040592f, 0.0339006f, 0.0455238f, 0.0080623f, -0.0042785f, -0.012738f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1738) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1738) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { -0.152f, 0.395f, -0.077f, -0.191f, -0.757f, 0.858f, -0.494f, -0.734f };
   int lda = 2;
   float B[] = { -0.166f, -0.413f, -0.373f, 0.915f, -0.824f, -0.066f, -0.114f, -0.921f, 0.862f, 0.312f, 0.221f, 0.699f };
   int ldb = 2;
   float B_expected[] = { 0.142569f, -0.0668709f, -0.0915f, -0.0373f, -0.0533385f, 0.0052516f, 0.0921f, -0.0114f, 0.0027525f, 0.0094961f, -0.0699f, 0.0221f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1739) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1739) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { -0.426f, 0.817f, -0.993f, -0.882f, 0.615f, 0.627f, -0.238f, -0.903f };
   int lda = 2;
   float B[] = { 0.895f, 0.849f, 0.811f, 0.402f, 0.074f, -0.493f, -0.548f, -0.82f, 0.323f, 0.301f, 0.612f, -0.092f };
   int ldb = 2;
   float B_expected[] = { -0.0369541f, -0.10749f, 0.246046f, 0.0030071f, -0.0270476f, 0.0371257f, -0.111428f, -0.111834f, -0.0135665f, -0.0383515f, 0.111452f, -0.0283989f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1740) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1740) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { -0.451f, -0.754f, -0.673f, 0.433f, -0.712f, -0.033f, -0.588f, 0.116f };
   int lda = 2;
   float B[] = { 0.787f, -0.377f, -0.854f, -0.464f, 0.118f, 0.231f, 0.362f, -0.457f, -0.076f, 0.373f, -0.286f, -0.468f };
   int ldb = 2;
   float B_expected[] = { 0.0377f, 0.0787f, -0.0130492f, -0.122041f, -0.0231f, 0.0118f, 0.0561369f, 0.0182563f, -0.0373f, -0.0076f, 0.0751937f, -0.0396361f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1741) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1741) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { -0.454f, 0.494f, 0.424f, -0.907f, 0.339f, -0.141f, 0.169f, 0.364f, -0.607f, 0.955f, -0.156f, 0.962f, -0.254f, 0.079f, 0.209f, 0.946f, 0.93f, 0.677f };
   int lda = 3;
   float B[] = { -0.99f, -0.484f, 0.915f, -0.383f, 0.228f, 0.797f, 0.597f, 0.765f, -0.629f, 0.002f, -0.89f, 0.077f };
   int ldb = 3;
   float B_expected[] = { 0.0269324f, 0.0688556f, -0.179902f, -0.104839f, -0.181106f, -0.0505677f, 0.0052392f, -0.0648948f, 0.0819028f, 0.132688f, 0.0961172f, -0.0473381f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1742) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1742) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { -0.008f, -0.654f, 0.174f, 0.448f, 0.388f, -0.108f, -0.479f, -0.708f, -0.035f, 0.816f, 0.487f, 0.22f, -0.482f, 0.57f, -0.317f, 0.203f, -0.547f, -0.415f };
   int lda = 3;
   float B[] = { 0.651f, 0.187f, 0.591f, -0.007f, 0.171f, -0.923f, -0.029f, -0.685f, -0.049f, 0.135f, 0.578f, 0.979f };
   int ldb = 3;
   float B_expected[] = { -0.0187f, 0.0651f, -0.0317186f, 0.0620498f, 0.0794141f, 0.0733141f, 0.0685f, -0.0029f, -0.0002818f, 0.0252834f, -0.0771317f, 0.0439205f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1743) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1743) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { 0.952f, 0.29f, 0.944f, 0.294f, -0.762f, -0.7f, -0.949f, 0.167f, 0.307f, 0.904f, -0.428f, -0.411f, 0.496f, 0.004f, -0.611f, -0.09f, -0.846f, 0.081f };
   int lda = 3;
   float B[] = { 0.782f, -0.035f, -0.441f, -0.791f, -0.09f, -0.56f, -0.438f, -0.691f, 0.88f, 0.545f, -0.55f, 0.595f };
   int ldb = 3;
   float B_expected[] = { -0.0592352f, 0.126282f, 0.0291241f, 0.0584267f, -0.046647f, 0.01215f, 0.0862177f, -0.14179f, -0.064879f, 0.016708f, 0.054792f, 0.0417105f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1744) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1744) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { -0.519f, 0.708f, -0.934f, -0.219f, 0.376f, -0.967f, 0.322f, -0.355f, 0.972f, -0.156f, -0.735f, 0.928f, 0.084f, -0.267f, -0.152f, 0.434f, 0.267f, 0.983f };
   int lda = 3;
   float B[] = { -0.54f, 0.149f, 0.574f, 0.742f, 0.704f, 0.459f, -0.9f, 0.04f, 0.538f, -0.858f, 0.467f, 0.686f };
   int ldb = 3;
   float B_expected[] = { -0.0034742f, 0.0089927f, -0.0977768f, 0.0267786f, -0.0459f, 0.0704f, 0.0494331f, -0.0808964f, 0.0759594f, 0.0169292f, -0.0686f, 0.0467f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1745) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1745) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { -0.541f, 0.67f, 0.014f, 0.446f, 0.086f, -0.525f, 0.033f, -0.932f, 0.977f, 0.321f, -0.651f, 0.027f, 0.409f, 0.328f, 0.359f, -0.615f, 0.419f, -0.25f };
   int lda = 3;
   float B[] = { -0.156f, 0.666f, -0.231f, 0.691f, 0.935f, -0.481f, -0.142f, -0.117f, 0.529f, 0.526f, 0.266f, 0.417f };
   int ldb = 2;
   float B_expected[] = { 0.0464826f, -0.0361824f, 0.0528601f, -0.0337999f, 0.0002432f, 0.168346f, -0.0078204f, 0.0535212f, 0.0438334f, 0.0110749f, -0.0360401f, -0.0228356f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1746) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1746) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { -0.459f, -0.349f, -0.335f, 0.008f, 0.866f, 0.978f, -0.869f, -0.361f, -0.711f, 0.712f, 0.207f, 0.305f, 0.766f, -0.262f, 0.012f, -0.333f, 0.617f, 0.91f };
   int lda = 3;
   float B[] = { -0.138f, -0.256f, -0.319f, -0.771f, 0.674f, -0.565f, -0.779f, -0.516f, -0.017f, -0.097f, -0.555f, 0.308f };
   int ldb = 2;
   float B_expected[] = { 0.0256f, -0.0138f, 0.0771f, -0.0319f, 0.0292718f, 0.0701506f, -0.0269158f, -0.078012f, 0.0488162f, -0.0369837f, -0.0054207f, -0.118253f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1747) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1747) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { -0.825f, -0.785f, -0.605f, -0.508f, 0.763f, -0.578f, -0.167f, -0.233f, 0.011f, -0.853f, 0.24f, 0.192f, 0.293f, -0.72f, -0.348f, 0.023f, -0.145f, -0.493f };
   int lda = 3;
   float B[] = { 0.305f, -0.255f, 0.882f, 0.883f, 0.088f, -0.473f, 0.135f, -0.063f, -0.671f, 0.473f, 0.874f, 0.548f };
   int ldb = 2;
   float B_expected[] = { -0.0961148f, -0.0983903f, 0.153836f, 0.0835432f, 0.0095579f, -0.0654357f, -0.018348f, 0.005229f, -0.0262218f, 0.0330484f, 0.0510342f, 0.0143434f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1748) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1748) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { 0.63f, 0.353f, 0.445f, 0.845f, 0.273f, -0.135f, 0.03f, 0.936f, 0.141f, 0.638f, -0.399f, 0.343f, -0.037f, -0.335f, -0.089f, 0.081f, 0.987f, -0.256f };
   int lda = 3;
   float B[] = { -0.567f, 0.803f, 0.168f, 0.744f, -0.328f, 0.835f, -0.852f, 0.702f, 0.21f, -0.618f, 0.666f, -0.303f };
   int ldb = 2;
   float B_expected[] = { -0.0700351f, -0.144464f, -0.0163821f, -0.0663417f, -0.115361f, -0.0199816f, -0.105134f, -0.10138f, 0.0618f, 0.021f, 0.0303f, 0.0666f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1749) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1749) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { -0.741f, 0.904f, -0.599f, 0.753f, -0.297f, 0.38f, -0.056f, -0.715f };
   int lda = 2;
   float B[] = { 0.646f, -0.447f, -0.147f, 0.314f, -0.713f, 0.187f, -0.589f, 0.287f, -0.809f, -0.293f, 0.418f, 0.778f };
   int ldb = 3;
   float B_expected[] = { -0.0915211f, -0.0074598f, 0.0365562f, -0.0174929f, 0.0783119f, 0.0359285f, -0.115925f, 0.0187826f, -0.0296066f, -0.031258f, 0.099134f, 0.0819138f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1750) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1750) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { -0.645f, 0.756f, 0.709f, -0.657f, -0.023f, -0.714f, 0.03f, 0.239f };
   int lda = 2;
   float B[] = { -0.16f, 0.254f, -0.68f, 0.183f, -0.402f, -0.259f, 0.104f, -0.09f, 0.944f, 0.729f, -0.378f, -0.792f };
   int ldb = 3;
   float B_expected[] = { -0.0254f, -0.016f, -0.0183f, -0.068f, 0.0259f, -0.0402f, -0.0195206f, 0.0157438f, -0.130551f, 0.0582111f, 0.0711517f, -0.0833181f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1751) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1751) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { 0.25f, -0.038f, 0.377f, -0.209f, 0.166f, -0.073f, -0.24f, 0.938f };
   int lda = 2;
   float B[] = { 0.26f, 0.696f, -0.183f, 0.668f, -0.08f, -0.938f, -0.837f, -0.509f, 0.781f, -0.063f, -0.953f, 0.227f };
   int ldb = 3;
   float B_expected[] = { -0.0140727f, -0.0084651f, -0.0106483f, 0.0104681f, 0.0124209f, -0.0197271f, 0.0662946f, 0.0678322f, -0.0747698f, -0.0128346f, 0.0948394f, 0.0015794f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1752) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1752) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { -0.668f, 0.804f, 0.608f, -0.682f, -0.513f, 0.521f, 0.878f, -0.664f };
   int lda = 2;
   float B[] = { -0.871f, 0.699f, 0.561f, 0.823f, -0.787f, 0.055f, -0.686f, 0.361f, -0.662f, -0.192f, -0.301f, -0.167f };
   int ldb = 3;
   float B_expected[] = { -0.0156401f, -0.0707163f, -0.0576594f, 0.100064f, 0.001615f, -0.054558f, -0.0361f, -0.0686f, 0.0192f, -0.0662f, 0.0167f, -0.0301f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1753) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1753) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { 0.091f, 0.189f, -0.221f, 0.749f, 0.354f, -0.397f, 0.105f, -0.944f };
   int lda = 2;
   float B[] = { 0.731f, -0.446f, 0.983f, 0.793f, 0.533f, 0.386f, -0.781f, -0.063f, 0.875f, -0.128f, -0.179f, -0.079f };
   int ldb = 2;
   float B_expected[] = { -0.0097573f, 0.0150815f, 0.129278f, 0.0933519f, -0.0135863f, -0.0024451f, -0.0655692f, 0.0200447f, -0.0153727f, 0.0103817f, 0.0232006f, 0.0165563f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1754) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1754) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { 0.676f, 0.644f, 0.03f, 0.456f, 0.002f, -0.909f, 0.984f, 0.771f };
   int lda = 2;
   float B[] = { 0.65f, 0.005f, -0.883f, -0.154f, -0.137f, -0.137f, 0.531f, -0.49f, 0.052f, 0.273f, -0.602f, 0.655f };
   int ldb = 2;
   float B_expected[] = { -0.0005f, 0.065f, 0.074484f, -0.0877155f, 0.0137f, -0.0137f, 0.0365741f, 0.0406193f, -0.0273f, 0.0052f, -0.0608278f, -0.0353739f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1755) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1755) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { 0.832f, -0.559f, 0.188f, -0.488f, -0.051f, -0.057f, 0.909f, 0.006f };
   int lda = 2;
   float B[] = { -0.408f, 0.303f, 0.03f, 0.529f, -0.584f, -0.976f, 0.443f, -0.762f, 0.43f, 0.812f, -0.075f, 0.06f };
   int ldb = 2;
   float B_expected[] = { -0.056498f, 0.0093713f, -0.0481041f, 0.0024096f, 0.0845016f, -0.132004f, 0.069f, 0.0407259f, -0.0483094f, 0.0826848f, -0.005409f, -0.0068535f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1756) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1756) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { -0.15f, -0.297f, 0.821f, -0.576f, -0.572f, 0.924f, 0.106f, -0.131f };
   int lda = 2;
   float B[] = { -0.271f, 0.793f, -0.232f, -0.967f, -0.466f, 0.37f, -0.745f, -0.156f, -0.091f, -0.877f, 0.595f, 0.448f };
   int ldb = 2;
   float B_expected[] = { -0.0132725f, -0.101846f, 0.0967f, -0.0232f, -0.0671044f, -0.11675f, 0.0156f, -0.0745f, 0.0851912f, 0.0655543f, -0.0448f, 0.0595f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1757) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1757) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { 0.002f, 0.0f, 0.626f, -0.148f, 0.874f, 0.229f, -0.227f, -0.55f, -0.895f, 0.586f, 0.934f, 0.618f, 0.958f, -0.543f, 0.49f, 0.671f, -0.871f, 0.227f };
   int lda = 3;
   float B[] = { -0.415f, 0.156f, -0.539f, -0.247f, -0.725f, 0.932f, 0.565f, 0.454f, -0.118f, 0.693f, -0.968f, -0.601f };
   int ldb = 3;
   float B_expected[] = { -0.0574005f, -0.122188f, -0.0327649f, -0.0625979f, 0.0976347f, 0.0419911f, 0.0294756f, -0.0678577f, 0.184894f, -0.0833182f, -0.0303735f, 0.0979555f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1758) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1758) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { 0.89f, 0.309f, -0.786f, 0.999f, 0.511f, 0.599f, 0.385f, -0.615f, 0.527f, -0.328f, -0.078f, -0.666f, 0.004f, -0.69f, -0.281f, -0.438f, 0.456f, 0.524f };
   int lda = 3;
   float B[] = { -0.648f, -0.189f, -0.295f, 0.477f, 0.509f, 0.685f, 0.875f, 0.277f, -0.34f, -0.632f, -0.453f, -0.798f };
   int ldb = 3;
   float B_expected[] = { 0.0203701f, -0.104287f, -0.0084576f, 0.0121508f, -0.0685f, 0.0509f, 0.0245033f, 0.202013f, 0.0268058f, -0.0836134f, 0.0798f, -0.0453f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1759) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1759) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { 0.772f, 0.686f, 0.693f, 0.803f, -0.328f, -0.627f, -0.869f, -0.656f, -0.055f, -0.366f, -0.981f, -0.151f, 0.147f, -0.368f, -0.824f, -0.454f, -0.445f, -0.794f };
   int lda = 3;
   float B[] = { -0.268f, -0.521f, -0.685f, -0.618f, 0.508f, 0.525f, -0.492f, -0.502f, -0.997f, 0.28f, 0.63f, 0.664f };
   int ldb = 3;
   float B_expected[] = { 0.058606f, 0.015051f, -0.0913257f, -0.0297397f, -0.0205282f, 0.0243534f, 0.0725056f, -0.0035452f, -0.110849f, 0.0255551f, 0.046652f, 0.0938454f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1760) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1760) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { -0.317f, -0.822f, 0.732f, 0.383f, 0.457f, 0.443f, 0.529f, -0.949f, -0.927f, -0.65f, -0.471f, -0.624f, -0.731f, 0.107f, -0.142f, 0.623f, 0.159f, -0.419f };
   int lda = 3;
   float B[] = { 0.292f, -0.665f, -0.93f, 0.517f, 0.123f, -0.181f, 0.325f, 0.954f, -0.988f, -0.128f, 0.637f, -0.997f };
   int ldb = 3;
   float B_expected[] = { 0.0665f, 0.0292f, 0.0111893f, -0.140662f, 0.0316445f, -0.0209328f, -0.0954f, 0.0325f, -0.0068241f, 0.0089271f, 0.225695f, 0.0517387f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1761) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1761) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { -0.809f, 0.393f, -0.015f, -0.273f, -0.956f, 0.49f, 0.365f, -0.386f, 0.941f, 0.992f, 0.297f, 0.761f, 0.425f, -0.605f, 0.672f, 0.725f, -0.077f, -0.628f };
   int lda = 3;
   float B[] = { 0.21f, 0.153f, 0.218f, -0.129f, 0.736f, -0.006f, 0.502f, -0.165f, 0.242f, 0.915f, 0.67f, 0.07f };
   int ldb = 2;
   float B_expected[] = { 0.0085068f, 0.069273f, 0.0439562f, 0.0320975f, -0.15148f, 0.0197777f, -0.0875509f, 0.103555f, 0.0222431f, 0.0555986f, 0.042615f, -0.000763f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1762) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1762) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { -0.187f, -0.508f, -0.987f, -0.861f, 0.519f, 0.752f, -0.117f, 0.972f, 0.068f, -0.752f, 0.344f, 0.074f, -0.343f, 0.0f, -0.876f, 0.857f, -0.148f, -0.933f };
   int lda = 3;
   float B[] = { 0.827f, 0.958f, 0.395f, 0.878f, 0.88f, -0.896f, -0.771f, -0.355f, -0.979f, 0.329f, -0.166f, -0.644f };
   int ldb = 2;
   float B_expected[] = { -0.180535f, 0.193075f, -0.0391015f, 0.0887205f, 0.202321f, 0.145565f, -0.0066882f, -0.0073676f, -0.0329f, -0.0979f, 0.0644f, -0.0166f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1763) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1763) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { -0.622f, 0.022f, -0.966f, 0.704f, 0.43f, -0.451f, -0.221f, 0.969f, 0.977f, 0.021f, -0.725f, -0.382f, 0.779f, 0.957f, 0.25f, 0.832f, 0.029f, -0.903f };
   int lda = 3;
   float B[] = { 0.315f, -0.297f, -0.864f, 0.519f, -0.601f, -0.119f, 0.028f, 0.072f, -0.171f, 0.648f, 0.159f, -0.623f };
   int ldb = 2;
   float B_expected[] = { -0.0191664f, -0.0189396f, 0.0341826f, 0.052599f, -0.0379778f, -0.067988f, 0.103868f, 0.0495092f, -0.0219287f, 0.0971955f, -0.0388294f, -0.0688205f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1764) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1764) imag");
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
   float alpha[2] = {0.0f, 0.1f};
   float A[] = { 0.106f, 0.87f, 0.21f, 0.463f, -0.496f, -0.981f, -0.354f, -0.604f, -0.149f, -0.384f, -0.958f, -0.502f, -0.579f, 0.736f, -0.322f, 0.028f, 0.193f, 0.14f };
   int lda = 3;
   float B[] = { -0.812f, 0.518f, 0.085f, -0.447f, -0.443f, 0.928f, -0.972f, 0.889f, 0.605f, -0.258f, -0.025f, 0.98f };
   int ldb = 2;
   float B_expected[] = { -0.0518f, -0.0812f, 0.0447f, 0.0085f, -0.0660824f, -0.0853354f, -0.0834485f, -0.0747189f, 0.0384994f, 0.240616f, -0.0754609f, 0.0871787f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1765) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1765) imag");
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
   float alpha[2] = {0.0f, 0.0f};
   float A[] = { -0.553f, 0.204f, -0.793f, -0.558f, 0.741f, 0.26f, 0.945f, -0.757f };
   int lda = 2;
   float B[] = { -0.515f, 0.532f, -0.321f, 0.326f, -0.81f, -0.924f, 0.474f, 0.985f, -0.03f, 0.406f, 0.923f, -0.956f };
   int ldb = 3;
   float B_expected[] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1766) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1766) imag");
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
   float alpha[2] = {0.0f, 0.0f};
   float A[] = { -0.41f, -0.804f, 0.988f, -0.715f, -0.281f, -0.89f, 0.389f, -0.408f };
   int lda = 2;
   float B[] = { 0.917f, 0.541f, -0.108f, -0.965f, 0.524f, 0.04f, -0.736f, -0.643f, -0.202f, 0.86f, 0.346f, -0.017f };
   int ldb = 3;
   float B_expected[] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1767) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1767) imag");
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
   float alpha[2] = {0.0f, 0.0f};
   float A[] = { 0.153f, -0.812f, -0.742f, -0.18f, 0.473f, 0.023f, -0.433f, 0.559f };
   int lda = 2;
   float B[] = { 0.078f, -0.691f, -0.717f, -0.637f, -0.016f, 0.375f, -0.902f, -0.343f, 0.155f, 0.563f, 0.419f, 0.451f };
   int ldb = 3;
   float B_expected[] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1768) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1768) imag");
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
   float alpha[2] = {0.0f, 0.0f};
   float A[] = { -0.288f, 0.241f, 0.593f, -0.597f, -0.469f, 0.735f, 0.193f, -0.104f };
   int lda = 2;
   float B[] = { -0.835f, 0.037f, -0.762f, 0.782f, -0.874f, -0.867f, -0.81f, -0.577f, 0.352f, 0.827f, 0.237f, -0.861f };
   int ldb = 3;
   float B_expected[] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1769) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1769) imag");
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
   float alpha[2] = {0.0f, 0.0f};
   float A[] = { 0.441f, -0.217f, 0.679f, 0.106f, -0.76f, -0.258f, -0.956f, -0.858f };
   int lda = 2;
   float B[] = { -0.802f, 0.163f, 0.293f, 0.54f, 0.228f, 0.071f, 0.942f, 0.345f, 0.591f, 0.654f, 0.382f, -0.892f };
   int ldb = 2;
   float B_expected[] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1770) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1770) imag");
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
   float alpha[2] = {0.0f, 0.0f};
   float A[] = { 0.916f, 0.909f, 0.834f, 0.38f, 0.391f, -0.412f, -0.714f, -0.456f };
   int lda = 2;
   float B[] = { -0.151f, 0.818f, 0.717f, -0.812f, -0.649f, -0.107f, -0.454f, 0.785f, 0.86f, 0.992f, -0.244f, -0.242f };
   int ldb = 2;
   float B_expected[] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1771) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1771) imag");
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
   float alpha[2] = {0.0f, 0.0f};
   float A[] = { -0.992f, 0.284f, -0.01f, 0.182f, 0.527f, -0.348f, -0.509f, 0.839f };
   int lda = 2;
   float B[] = { 0.504f, -0.782f, -0.88f, 0.079f, 0.216f, 0.525f, 0.198f, 0.851f, -0.102f, -0.046f, 0.079f, -0.045f };
   int ldb = 2;
   float B_expected[] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1772) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1772) imag");
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
   float alpha[2] = {0.0f, 0.0f};
   float A[] = { 0.985f, 0.068f, -0.095f, -0.575f, -0.607f, 0.893f, 0.085f, 0.145f };
   int lda = 2;
   float B[] = { -0.149f, 0.592f, 0.588f, -0.62f, -0.409f, -0.344f, 0.263f, 0.759f, -0.026f, -0.609f, 0.507f, -0.084f };
   int ldb = 2;
   float B_expected[] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1773) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1773) imag");
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
   float alpha[2] = {0.0f, 0.0f};
   float A[] = { 0.36f, 0.508f, -0.771f, -0.442f, -0.671f, -0.691f, -0.771f, 0.113f, 0.282f, 0.312f, 0.564f, -0.568f, -0.743f, 0.912f, -0.395f, 0.503f, -0.167f, -0.581f };
   int lda = 3;
   float B[] = { -0.018f, 0.574f, -0.144f, -0.758f, 0.53f, 0.623f, -0.771f, -0.733f, 0.932f, -0.192f, 0.997f, 0.773f };
   int ldb = 3;
   float B_expected[] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1774) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1774) imag");
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
   float alpha[2] = {0.0f, 0.0f};
   float A[] = { -0.627f, 0.511f, -0.246f, -0.091f, 0.66f, -0.983f, 0.99f, 0.057f, -0.259f, 0.18f, 0.606f, 0.058f, -0.238f, 0.717f, 0.358f, -0.851f, -0.71f, -0.683f };
   int lda = 3;
   float B[] = { -0.907f, 0.956f, 0.56f, -0.057f, 0.054f, -0.77f, 0.868f, -0.843f, 0.645f, -0.554f, -0.958f, 0.988f };
   int ldb = 3;
   float B_expected[] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1775) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1775) imag");
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
   float alpha[2] = {0.0f, 0.0f};
   float A[] = { -0.882f, 0.431f, -0.868f, -0.098f, -0.006f, -0.639f, 0.757f, -0.009f, -0.821f, 0.45f, 0.347f, 0.801f, 0.314f, 0.936f, -0.725f, 0.956f, 0.536f, 0.771f };
   int lda = 3;
   float B[] = { 0.38f, -0.435f, 0.977f, 0.296f, -0.624f, -0.53f, 0.73f, -0.837f, 0.105f, 0.189f, 0.362f, -0.664f };
   int ldb = 3;
   float B_expected[] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1776) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1776) imag");
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
   float alpha[2] = {0.0f, 0.0f};
   float A[] = { 0.595f, -0.775f, 0.75f, 0.16f, -0.572f, 0.658f, 0.216f, 0.557f, -0.279f, 0.095f, -0.495f, 0.503f, 0.071f, -0.03f, -0.116f, 0.78f, -0.104f, 0.073f };
   int lda = 3;
   float B[] = { 0.948f, 0.749f, -0.854f, 0.972f, 0.704f, 0.187f, 0.347f, 0.303f, -0.865f, 0.123f, -0.041f, 0.152f };
   int ldb = 3;
   float B_expected[] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1777) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1777) imag");
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
   float alpha[2] = {0.0f, 0.0f};
   float A[] = { -0.617f, -0.331f, -0.074f, 0.719f, -0.469f, -0.852f, 0.25f, -0.175f, -0.719f, -0.613f, -0.321f, 0.973f, -0.337f, -0.35f, 0.607f, -0.553f, 0.688f, 0.463f };
   int lda = 3;
   float B[] = { 0.568f, -0.471f, -0.947f, -0.205f, 0.835f, -0.859f, 0.27f, -0.599f, 0.171f, -0.514f, 0.939f, 0.176f };
   int ldb = 2;
   float B_expected[] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1778) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1778) imag");
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
   float alpha[2] = {0.0f, 0.0f};
   float A[] = { 0.99f, -0.857f, 0.728f, -0.31f, -0.506f, -0.393f, 0.97f, 0.282f, 0.375f, -0.286f, -0.496f, -0.057f, 0.186f, -0.34f, 0.608f, -0.52f, 0.921f, -0.875f };
   int lda = 3;
   float B[] = { -0.929f, 0.885f, 0.864f, -0.548f, 0.393f, 0.391f, 0.033f, 0.186f, 0.949f, -0.435f, 0.986f, -0.995f };
   int ldb = 2;
   float B_expected[] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1779) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1779) imag");
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
   float alpha[2] = {0.0f, 0.0f};
   float A[] = { -0.101f, -0.92f, 0.969f, -0.017f, -0.016f, -0.024f, -0.11f, 0.219f, -0.287f, -0.937f, 0.619f, 0.166f, -0.068f, 0.753f, 0.374f, 0.076f, 0.79f, -0.64f };
   int lda = 3;
   float B[] = { 0.255f, 0.564f, -0.478f, -0.818f, -0.043f, 0.224f, -0.268f, 0.253f, 0.021f, 0.654f, 0.98f, -0.774f };
   int ldb = 2;
   float B_expected[] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1780) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1780) imag");
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
   float alpha[2] = {0.0f, 0.0f};
   float A[] = { -0.068f, -0.603f, -0.055f, 0.14f, 0.664f, 0.987f, 0.861f, -0.691f, -0.897f, -0.778f, 0.516f, -0.073f, -0.156f, -0.42f, 0.57f, 0.628f, 0.116f, 0.344f };
   int lda = 3;
   float B[] = { 0.922f, 0.39f, -0.724f, 0.421f, 0.418f, 0.92f, -0.222f, 0.835f, 0.417f, -0.392f, 0.012f, -0.346f };
   int ldb = 2;
   float B_expected[] = { 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f };
   cblas_ctrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], flteps, "ctrmm(case 1781) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], flteps, "ctrmm(case 1781) imag");
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
   double A[] = { 0.904, 0.243, 0.206, 0.68, -0.946, 0.946, -0.675, 0.729 };
   int lda = 2;
   double B[] = { 0.427, 0.116, 0.916, -0.384, -0.372, -0.754, 0.148, 0.089, -0.924, 0.974, -0.307, -0.55 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1782) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1782) imag");
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
   double A[] = { -0.898, 0.709, 0.719, -0.207, -0.841, -0.017, 0.202, -0.385 };
   int lda = 2;
   double B[] = { 0.308, 0.507, -0.838, 0.594, -0.811, 0.152, 0.118, -0.024, -0.632, 0.992, -0.942, 0.901 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1783) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1783) imag");
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
   double A[] = { -0.849, 0.455, -0.273, -0.668, 0.196, -0.985, -0.39, 0.564 };
   int lda = 2;
   double B[] = { -0.874, 0.188, -0.039, 0.692, 0.33, 0.119, 0.012, 0.425, 0.787, -0.918, 0.739, -0.871 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1784) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1784) imag");
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
   double A[] = { -0.325, 0.28, 0.902, -0.603, 0.091, -0.92, 0.209, -0.009 };
   int lda = 2;
   double B[] = { -0.202, -0.53, -0.88, -0.688, -0.215, 0.837, 0.917, 0.755, 0.477, 0.892, -0.524, -0.741 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1785) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1785) imag");
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
   double A[] = { -0.756, 0.874, 0.56, 0.157, -0.831, -0.991, -0.531, 0.813 };
   int lda = 2;
   double B[] = { 0.271, 0.783, -0.861, 0.635, -0.088, 0.434, 0.256, -0.34, -0.724, -0.277, -0.604, 0.986 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1786) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1786) imag");
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
   double A[] = { -0.371, -0.609, -0.812, -0.818, 0.45, -0.41, -0.704, -0.917 };
   int lda = 2;
   double B[] = { -0.268, 0.929, 0.82, 0.253, -0.883, 0.497, -0.265, 0.623, 0.131, -0.946, -0.365, 0.333 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1787) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1787) imag");
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
   double A[] = { -0.265, 0.8, -0.676, -0.592, 0.78, -0.838, -0.651, 0.115 };
   int lda = 2;
   double B[] = { 0.942, 0.692, -0.516, 0.378, 0.028, 0.265, 0.289, -0.721, -0.25, -0.952, 0.463, -0.34 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1788) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1788) imag");
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
   double A[] = { -0.852, -0.478, 0.16, 0.824, 0.073, 0.962, 0.509, -0.58 };
   int lda = 2;
   double B[] = { -0.789, 0.015, -0.779, -0.565, 0.048, -0.095, -0.272, 0.405, 0.272, 0.082, -0.693, -0.365 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1789) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1789) imag");
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
   double A[] = { 0.251, 0.28, -0.092, 0.724, 0.928, -0.309, -0.222, -0.791, 0.113, -0.528, 0.148, 0.421, -0.833, 0.371, 0.354, 0.616, 0.313, 0.323 };
   int lda = 3;
   double B[] = { -0.769, -0.059, -0.068, 0.945, 0.938, -0.358, -0.17, 0.751, -0.248, -0.321, -0.818, 0.183 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1790) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1790) imag");
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
   double A[] = { -0.707, -0.802, 0.13, -0.19, -0.564, -0.74, 0.118, -0.194, -0.124, -0.421, 0.665, 0.308, 0.505, -0.278, 0.588, 0.957, -0.727, 0.976 };
   int lda = 3;
   double B[] = { 0.153, -0.09, -0.4, 0.669, 0.689, -0.238, -0.259, 0.891, 0.993, 0.996, -0.829, -0.736 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1791) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1791) imag");
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
   double A[] = { 0.83, 0.316, -0.099, 0.824, 0.767, 0.662, 0.244, 0.872, 0.35, 0.969, -0.084, 0.907, -0.752, -0.675, 0.129, -0.649, -0.539, 0.969 };
   int lda = 3;
   double B[] = { -0.145, 0.254, -0.497, -0.713, -0.742, 0.183, 0.272, -0.858, -0.606, -0.605, -0.807, 0.686 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1792) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1792) imag");
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
   double A[] = { -0.091, 0.658, -0.834, -0.171, -0.126, -0.268, 0.879, -0.431, 0.678, -0.749, 0.136, -0.757, -0.578, 0.456, 0.978, -0.315, 0.333, 0.327 };
   int lda = 3;
   double B[] = { 0.963, -0.859, 0.599, 0.856, -0.924, 0.382, -0.531, 0.567, -0.454, 0.018, 0.97, 0.578 };
   int ldb = 3;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1793) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1793) imag");
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
   double A[] = { -0.849, -0.819, 0.673, 0.574, -0.869, -0.969, -0.338, -0.097, -0.601, 0.903, 0.634, 0.313, 0.228, -0.028, 0.419, -0.762, 0.21, -0.532 };
   int lda = 3;
   double B[] = { -0.283, 0.999, -0.356, -0.459, 0.508, -0.132, -0.804, 0.173, 0.779, -0.427, 0.019, 0.347 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1794) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1794) imag");
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
   double A[] = { -0.117, -0.663, -0.95, -0.273, -0.497, -0.037, 0.084, -0.831, 0.023, -0.241, 0.063, -0.023, -0.498, -0.137, -0.77, 0.457, -0.021, -0.69 };
   int lda = 3;
   double B[] = { 0.308, -0.004, 0.013, 0.354, 0.077, -0.944, -0.877, 0.741, -0.807, -0.3, 0.891, -0.056 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1795) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1795) imag");
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
   double A[] = { -0.964, -0.653, 0.379, 0.994, -0.378, -0.409, 0.24, -0.333, 0.558, -0.099, -0.402, -0.812, 0.421, 0.823, -0.771, 0.998, 0.697, 0.253 };
   int lda = 3;
   double B[] = { 0.34, 0.479, 0.539, -0.133, 0.876, -0.347, 0.706, -0.623, 0.399, 0.903, -0.7, -0.088 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1796) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1796) imag");
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
   double A[] = { -0.104, 0.643, -0.253, -0.988, -0.051, -0.805, 0.451, -0.421, -0.177, -0.534, -0.714, -0.581, -0.177, -0.582, -0.57, 0.259, -0.66, -0.864 };
   int lda = 3;
   double B[] = { 0.636, -0.365, -0.107, -0.279, 0.425, 0.976, 0.657, 0.294, 0.827, 0.187, 0.353, 0.31 };
   int ldb = 2;
   double B_expected[] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1797) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1797) imag");
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
   double A[] = { 0.273, 0.812, 0.295, -0.415, -0.227, 0.901, 0.623, 0.786 };
   int lda = 2;
   double B[] = { -0.539, -0.551, -0.969, 0.09, -0.581, -0.594, -0.833, 0.457, -0.284, 0.434, -0.459, -0.662 };
   int ldb = 3;
   double B_expected[] = { -0.0312704, 0.2064538, 0.1775109, 0.1949157, -0.0337211, 0.2225517, 0.410638, -0.033917, 0.182384, -0.219409, 0.1257905, 0.1938415 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1798) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1798) imag");
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
   double A[] = { 0.323, 0.02, 0.718, 0.152, 0.665, 0.289, 0.317, 0.705 };
   int lda = 2;
   double B[] = { 0.448, -0.75, 0.851, 0.172, -0.244, 0.398, 0.602, 0.31, -0.017, 0.181, -0.119, 0.402 };
   int ldb = 3;
   double B_expected[] = { -0.0594, 0.2698, -0.2725, 0.0335, 0.0334, -0.1438, -0.2952588, 0.1518876, -0.213747, -0.073367, 0.0413388, -0.2306716 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1799) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1799) imag");
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
   double A[] = { -0.578, 0.018, -0.093, 0.964, 0.414, -0.729, 0.696, 0.874 };
   int lda = 2;
   double B[] = { -0.735, 0.788, -0.942, -0.71, -0.254, 0.265, 0.304, 0.218, 0.247, -0.172, 0.419, 0.448 };
   int ldb = 3;
   double B_expected[] = { -0.1486214, 0.2495598, -0.1744531, 0.0107667, -0.1648579, 0.1475263, -0.048058, -0.123122, -0.1062886, 0.0033742, -0.037823, -0.213397 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1800) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1800) imag");
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
   double A[] = { 0.358, -0.773, -0.065, 0.532, -0.319, 0.455, 0.578, 0.493 };
   int lda = 2;
   double B[] = { 0.744, -0.958, 0.162, 0.555, -0.131, 0.971, -0.467, 0.175, -0.794, 0.191, 0.361, 0.882 };
   int ldb = 3;
   double B_expected[] = { -0.1213734, 0.4492278, -0.1117944, -0.0070022, 0.108851, -0.320916, 0.1226, -0.0992, 0.2191, -0.1367, -0.1965, -0.2285 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1801) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1801) imag");
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
   double A[] = { -0.354, -0.504, -0.177, 0.186, -0.762, -0.506, 0.758, -0.994 };
   int lda = 2;
   double B[] = { -0.944, 0.562, 0.142, 0.742, 0.632, -0.627, -0.101, 0.476, 0.476, 0.675, 0.912, -0.33 };
   int ldb = 2;
   double B_expected[] = { -0.21291, -0.021306, -0.601736, 0.043676, 0.1715778, -0.0250026, 0.0587596, -0.2259812, -0.0036234, 0.1608258, 0.0885532, 0.6077736 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1802) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1802) imag");
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
   double A[] = { -0.001, 0.015, 0.942, 0.497, -0.104, 0.803, 0.679, 0.026 };
   int lda = 2;
   double B[] = { 0.889, -0.216, -0.912, -0.263, -0.329, 0.681, 0.332, -0.5, -0.484, 0.741, -0.728, -0.912 };
   int ldb = 2;
   double B_expected[] = { -0.2451, 0.1537, 0.2019693, -0.2251001, 0.0306, -0.2372, 0.1376892, 0.2324406, 0.0711, -0.2707, 0.5195777, 0.2860461 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1803) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1803) imag");
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
   double A[] = { -0.563, 0.394, -0.902, -0.27, 0.461, 0.939, -0.597, 0.803 };
   int lda = 2;
   double B[] = { 0.535, -0.111, 0.379, -0.036, 0.803, -0.341, 0.667, 0.001, 0.775, 0.714, 0.908, -0.508 };
   int ldb = 2;
   double B_expected[] = { 0.1623722, -0.1219324, 0.0266236, -0.1174842, 0.2429924, -0.1901218, 0.0662002, -0.2004014, 0.4905027, -0.2023089, -0.0629944, -0.3231352 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1804) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1804) imag");
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
   double A[] = { 0.159, 0.032, 0.785, 0.049, -0.128, 0.132, -0.735, -0.235 };
   int lda = 2;
   double B[] = { -0.331, -0.257, -0.725, 0.689, -0.793, 0.398, 0.127, -0.098, -0.498, -0.307, -0.019, 0.517 };
   int ldb = 2;
   double B_expected[] = { 0.2553318, -0.1678906, 0.1486, -0.2792, 0.1738216, -0.1670382, -0.0283, 0.0421, 0.151683, -0.083199, -0.046, -0.157 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1805) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1805) imag");
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
   double A[] = { -0.416, -0.424, -0.088, 0.614, -0.371, 0.983, -0.737, -0.647, 0.321, -0.518, 0.058, -0.533, 0.153, 0.283, 0.342, 0.993, -0.071, 0.225 };
   int lda = 3;
   double B[] = { -0.09, -0.844, -0.707, 0.903, 0.632, -0.294, -0.558, 0.74, -0.99, -0.855, -0.189, 0.543 };
   int ldb = 3;
   double B_expected[] = { 0.1668304, -0.2576208, -0.0664464, -0.0785782, -0.0226908, -0.0467944, -0.1091876, 0.3667652, 0.1076073, -0.1594011, 0.0407346, 0.0134478 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1806) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1806) imag");
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
   double A[] = { -0.67, -0.423, -0.165, 0.157, -0.43, 0.674, -0.35, 0.434, 0.972, -0.116, -0.029, 0.316, 0.914, 0.321, 0.132, 0.034, -0.907, -0.401 };
   int lda = 3;
   double B[] = { -0.396, 0.71, -0.588, 0.709, -0.024, -0.704, -0.988, 0.656, 0.665, -0.085, -0.778, 0.264 };
   int ldb = 3;
   double B_expected[] = { -0.1010812, -0.2287206, 0.0372688, -0.2530336, 0.0776, 0.2088, 0.264679, -0.133739, -0.147391, 0.161965, 0.207, -0.157 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1807) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1807) imag");
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
   double A[] = { 0.756, -0.149, -0.706, -0.162, -0.145, 0.67, 0.416, -0.27, -0.916, 0.995, -0.863, -0.25, -0.079, 0.248, -0.191, -0.195, 0.981, 0.834 };
   int lda = 3;
   double B[] = { 0.329, 0.921, -0.018, -0.02, 0.095, -0.892, -0.105, -0.799, -0.583, 0.564, -0.436, 0.965 };
   int ldb = 3;
   double B_expected[] = { -0.1805114, -0.1555812, -0.1560482, -0.0462226, -0.0967127, 0.2921239, 0.1183692, 0.1566766, 0.2260429, 0.3915667, 0.1788155, -0.2682995 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1808) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1808) imag");
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
   double A[] = { 0.552, -0.668, -0.013, 0.088, -0.766, 0.977, 0.088, -0.06, -0.311, 0.872, -0.328, -0.01, 0.659, -0.327, -0.276, 0.553, -0.734, -0.079 };
   int lda = 3;
   double B[] = { -0.87, 0.728, 0.997, -0.36, -0.046, -0.505, 0.082, -0.787, 0.414, 0.965, -0.048, -0.591 };
   int ldb = 3;
   double B_expected[] = { 0.1882, -0.3054, -0.2648624, 0.1695328, 0.0462155, -0.3187195, 0.0541, 0.2443, -0.2012812, -0.2298476, 0.3871505, 0.2622315 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1809) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1809) imag");
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
   double A[] = { 0.349, -0.072, 0.545, 0.212, -0.306, -0.009, 0.757, -0.925, 0.159, 0.308, 0.476, 0.1, 0.725, -0.757, -0.245, 0.571, 0.515, 0.993 };
   int lda = 3;
   double B[] = { 0.865, 0.501, 0.165, -0.63, -0.513, 0.351, -0.521, -0.062, 0.54, -0.634, -0.719, 0.216 };
   int ldb = 2;
   double B_expected[] = { -0.054193, 0.023274, 0.1487731, -0.3509657, -0.0481592, -0.1044386, 0.0666567, 0.1890461, -0.2932696, 0.0278532, 0.2357046, 0.1223408 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1810) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1810) imag");
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
   double A[] = { 0.941, -0.496, 0.492, 0.356, 0.353, 0.346, -0.519, -0.86, -0.677, -0.154, 0.313, 0.228, -0.56, -0.451, -0.78, 0.174, -0.663, 0.22 };
   int lda = 3;
   double B[] = { 0.162, -0.345, 0.188, 0.578, -0.675, 0.775, -0.018, 0.198, -0.222, -0.52, 0.672, -0.438 };
   int ldb = 2;
   double B_expected[] = { -0.3430472, 0.0394834, 0.0185782, -0.1505014, 0.0092108, -0.3837276, 0.0741276, -0.2435652, 0.1186, 0.1338, -0.1578, 0.1986 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1811) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1811) imag");
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
   double A[] = { 0.592, 0.708, 0.442, 0.212, 0.815, -0.638, 0.55, -0.512, -0.487, 0.181, 0.708, -0.126, 0.408, -0.51, 0.175, 0.114, -0.919, -0.268 };
   int lda = 3;
   double B[] = { 0.858, -0.004, 0.59, -0.395, -0.943, 0.824, 0.01, 0.455, -0.775, 0.062, -0.644, 0.03 };
   int ldb = 2;
   double B_expected[] = { -0.21374, -0.130452, -0.20707, 0.00773, -0.16787, 0.186571, -0.05026, 0.106515, -0.2887485, -0.0045065, -0.2446935, 0.1590455 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1812) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1812) imag");
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
   double A[] = { -0.988, -0.915, 0.963, 0.103, 0.921, 0.555, 0.846, 0.148, -0.43, 0.336, -0.371, 0.381, -0.487, 0.717, 0.881, -0.777, 0.774, -0.962 };
   int lda = 3;
   double B[] = { -0.805, 0.605, 0.481, 0.163, -0.057, -0.017, -0.886, 0.809, 0.875, 0.905, 0.095, 0.894 };
   int ldb = 2;
   double B_expected[] = { 0.181, -0.262, -0.1606, -0.0008, 0.220089, -0.234263, 0.0303246, -0.3486122, -0.0476352, -0.3174616, -0.2077412, -0.1552106 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1813) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1813) imag");
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
   double A[] = { -0.513, -0.385, -0.524, 0.726, 0.823, 0.839, -0.355, -0.881 };
   int lda = 2;
   double B[] = { -0.707, 0.016, 0.481, 0.935, 0.052, 0.719, 0.277, 0.169, 0.894, 0.352, -0.216, -0.741 };
   int ldb = 3;
   double B_expected[] = { -0.078919, 0.119774, 0.2114654, 0.0276682, 0.12593, 0.074299, -0.109352, -0.193196, 0.077864, 0.032876, -0.3330992, 0.2249494 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1814) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1814) imag");
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
   double A[] = { -0.136, -0.37, 0.669, -0.731, -0.4, 0.638, 0.833, -0.29 };
   int lda = 2;
   double B[] = { -0.861, -0.278, 0.941, 0.822, 0.88, 0.501, 0.911, -0.502, 0.573, -0.498, -0.517, -0.518 };
   int ldb = 3;
   double B_expected[] = { 0.2861, -0.0027, -0.3645, -0.1525, -0.3141, -0.0623, -0.0297254, 0.4490328, -0.254473, -0.161772, 0.0423084, -0.1675858 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1815) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1815) imag");
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
   double A[] = { 0.641, -0.058, 0.246, 0.884, -0.686, 0.123, -0.869, 0.891 };
   int lda = 2;
   double B[] = { 0.107, -0.333, 0.556, 0.124, 0.206, 0.049, -0.573, -0.9, -0.417, -0.734, -0.719, 0.76 };
   int ldb = 3;
   double B_expected[] = { -0.1591469, -0.1071617, -0.2301499, -0.1454657, -0.1758188, 0.1884616, -0.0380754, -0.4181892, -0.013453, -0.33198, -0.3886102, 0.1361404 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1816) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1816) imag");
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
   double A[] = { 0.083, 0.441, 0.995, 0.338, -0.988, -0.828, -0.254, -0.036 };
   int lda = 2;
   double B[] = { -0.792, 0.552, 0.033, -0.178, -0.225, 0.553, 0.348, 0.229, -0.151, -0.594, 0.711, -0.335 };
   int ldb = 3;
   double B_expected[] = { 0.3362416, -0.3167112, -0.2305904, -0.0177512, 0.0477576, -0.5068152, -0.1273, -0.0339, 0.1047, 0.1631, -0.1798, 0.1716 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1817) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1817) imag");
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
   double A[] = { 0.105, 0.584, -0.33, -0.182, -0.096, -0.257, 0.327, -0.123 };
   int lda = 2;
   double B[] = { -0.249, -0.274, -0.197, -0.899, 0.85, -0.318, 0.596, -0.237, 0.179, 0.046, -0.859, -0.459 };
   int ldb = 2;
   double B_expected[] = { 0.0441837, -0.0536099, -0.0065547, 0.1208159, 0.0819176, 0.1492908, -0.0917294, -0.0510192, -0.0037271, 0.0344777, 0.0974489, 0.0389047 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1818) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1818) imag");
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
   double A[] = { -0.972, 0.794, -0.968, -0.406, -0.2, -0.512, 0.436, 0.161 };
   int lda = 2;
   double B[] = { 0.817, -0.17, -0.613, -0.565, -0.494, 0.129, -0.593, -0.516, -0.695, -0.42, 0.848, 0.122 };
   int ldb = 2;
   double B_expected[] = { -0.2281, 0.1327, 0.2180776, -0.0351272, 0.1353, -0.0881, 0.2475472, 0.1823936, 0.2505, 0.0565, -0.345628, 0.165156 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1819) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1819) imag");
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
   double A[] = { 0.373, -0.316, -0.052, 0.025, -0.878, 0.612, 0.486, 0.953 };
   int lda = 2;
   double B[] = { -0.626, 0.408, 0.536, 0.66, -0.666, -0.127, 0.622, 0.036, -0.761, 0.773, -0.137, 0.074 };
   int ldb = 2;
   double B_expected[] = { 0.1214746, -0.0093742, -0.247838, 0.145962, 0.0994439, 0.0586017, -0.043453, 0.206241, 0.1510011, -0.0661437, -0.0178345, -0.0495635 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1820) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1820) imag");
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
   double A[] = { 0.621, -0.252, -0.942, 0.073, 0.416, -0.724, -0.972, 0.028 };
   int lda = 2;
   double B[] = { -0.006, 0.427, 0.292, -0.212, -0.319, -0.08, -0.401, 0.465, -0.493, -0.529, 0.003, -0.19 };
   int ldb = 2;
   double B_expected[] = { 0.0284232, -0.2112704, -0.0664, 0.0928, 0.0210696, 0.1558958, 0.0738, -0.1796, 0.1879327, 0.0541021, 0.0181, 0.0573 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1821) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1821) imag");
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
   double A[] = { -0.415, 0.215, 0.507, 0.094, 0.697, 0.633, 0.206, -0.383, -0.974, 0.734, -0.533, -0.15, -0.982, -0.232, -0.297, 0.501, -0.092, 0.663 };
   int lda = 3;
   double B[] = { 0.812, 0.323, 0.294, -0.423, -0.85, 0.043, -0.338, -0.568, 0.976, -0.375, 0.913, -0.119 };
   int ldb = 3;
   double B_expected[] = { 0.2153111, -0.0775367, 0.0404927, -0.0287599, -0.0879721, -0.1572073, -0.2481947, 0.2941819, 0.5234716, -0.1242382, 0.108305, 0.162022 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1822) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1822) imag");
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
   double A[] = { 0.827, -0.754, 0.719, 0.88, -0.942, -0.152, 0.051, 0.033, -0.603, -0.557, 0.668, 0.024, 0.082, 0.458, 0.733, 0.669, 0.722, -0.661 };
   int lda = 3;
   double B[] = { -0.523, 0.365, -0.811, -0.632, -0.06, 0.151, -0.962, -0.71, -0.543, 0.8, -0.264, 0.994 };
   int ldb = 3;
   double B_expected[] = { 0.4413193, -0.3047431, 0.307206, 0.074162, 0.0029, -0.0513, 0.2285887, 0.1349491, 0.061616, -0.510648, -0.0202, -0.3246 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1823) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1823) imag");
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
   double A[] = { 0.958, 0.948, -0.161, -0.34, -0.184, 0.43, -0.045, -0.465, -0.278, 0.461, 0.584, 0.003, -0.794, -0.778, -0.65, -0.91, 0.24, -0.944 };
   int lda = 3;
   double B[] = { 0.279, 0.041, -0.033, 0.332, 0.788, 0.611, -0.644, -0.133, 0.247, 0.06, 0.125, -0.407 };
   int ldb = 3;
   double B_expected[] = { -0.0693236, 0.0981792, -0.0442625, -0.0021815, 0.1936084, -0.3409328, 0.174601, -0.219233, 0.0274565, 0.1321885, -0.2252264, 0.1381888 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1824) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1824) imag");
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
   double A[] = { -0.983, -0.795, -0.115, -0.542, 0.837, 0.518, -0.164, 0.776, -0.453, -0.28, 0.135, -0.377, -0.199, -0.965, 0.784, -0.39, -0.499, 0.257 };
   int lda = 3;
   double B[] = { -0.712, 0.364, -0.28, 0.05, 0.314, 0.748, -0.719, 0.619, 0.474, -0.906, -0.859, 0.943 };
   int ldb = 3;
   double B_expected[] = { 0.1772, -0.1804, -0.0900512, -0.1509216, 0.0485292, 0.0109956, 0.1538, -0.2576, -0.2767208, 0.2420976, 0.2164354, 0.0610082 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1825) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1825) imag");
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
   double A[] = { 0.105, 0.503, -0.17, 0.2, -0.861, -0.279, -0.231, 0.058, 0.699, 0.437, 0.578, 0.462, 0.473, -0.793, -0.34, -0.162, -0.128, -0.844 };
   int lda = 3;
   double B[] = { -0.802, 0.292, -0.155, -0.916, -0.099, -0.082, 0.057, 0.215, 0.94, 0.911, -0.714, 0.41 };
   int ldb = 2;
   double B_expected[] = { -0.1044001, -0.5102243, 0.3865174, 0.0189802, 0.1888166, -0.0057672, -0.0800722, 0.0699214, 0.199086, -0.291946, 0.141904, 0.171064 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1826) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1826) imag");
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
   double A[] = { 0.468, 0.378, -0.498, 0.251, 0.777, -0.543, -0.913, 0.095, 0.779, -0.933, 0.068, -0.669, 0.715, 0.03, 0.012, 0.392, -0.785, -0.056 };
   int lda = 3;
   double B[] = { 0.143, -0.242, -0.379, -0.831, -0.46, -0.663, -0.735, -0.098, -0.861, -0.894, 0.772, -0.059 };
   int ldb = 2;
   double B_expected[] = { 0.0633681, 0.0476643, -0.1761819, 0.3044093, 0.2798556, 0.0187868, 0.2647924, 0.0455132, 0.3477, 0.1821, -0.2257, 0.0949 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1827) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1827) imag");
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
   double A[] = { -0.764, 0.908, 0.899, 0.119, -0.447, 0.279, 0.338, 0.73, -0.74, -0.366, -0.572, 0.583, 0.75, 0.519, 0.603, 0.831, 0.697, 0.822 };
   int lda = 3;
   double B[] = { 0.399, 0.572, -0.489, 0.964, -0.167, -0.104, 0.75, -0.199, 0.777, 0.503, -0.025, -0.386 };
   int ldb = 2;
   double B_expected[] = { 0.015568, 0.261244, -0.345424, 0.212636, -0.2247824, -0.0859342, 0.1074596, -0.4846822, -0.2415227, 0.2465939, 0.2042976, 0.2206978 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1828) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1828) imag");
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
   double A[] = { 0.432, 0.063, 0.065, -0.546, 0.099, 0.892, 0.48, -0.085, 0.746, -0.541, -0.739, -0.207, 0.695, 0.765, 0.197, -0.86, 0.621, -0.653 };
   int lda = 3;
   double B[] = { 0.182, 0.731, 0.571, 0.01, -0.357, -0.612, 0.581, 0.756, -0.911, -0.225, 0.438, 0.546 };
   int ldb = 2;
   double B_expected[] = { -0.1277, -0.2011, -0.1723, 0.0541, 0.2698001, 0.0651043, -0.2906381, -0.2592593, -0.0512125, -0.0040605, 0.0647965, 0.1119875 };
   cblas_ztrmm(order, side, uplo, trans, diag, M, N, alpha, A, lda, B, ldb);
   {
     int i;
     for (i = 0; i < 6; i++) {
       gsl_test_rel(B[2*i], B_expected[2*i], dbleps, "ztrmm(case 1829) real");
       gsl_test_rel(B[2*i+1], B_expected[2*i+1], dbleps, "ztrmm(case 1829) imag");
     };
   };
  };


}
