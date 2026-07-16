#include <stdint.h>

void CEIGEN_LITE_smatmul(const long m, const long n, const long k,
                    float* a, char a_layout,
                    float* b, char b_layout,
                    float* c, char c_layout);
void CEIGEN_LITE_dmatmul(const long m, const long n, const long k,
                    double* a, char a_layout,
                    double* b, char b_layout,
                    double* c, char c_layout);

void CEIGEN_LITE_spartialPivLu(const long m, const long n,
                               float* A, char A_layout,
                               float* b, float* x);
void CEIGEN_LITE_dpartialPivLu(const long m, const long n,
                               double* A, char A_layout,
                               double* b, double* x);
void CEIGEN_LITE_scolPivHouseholderQr(const long m, const long n,
                                      float* A, char A_layout,
                                      float* b, float* x);
void CEIGEN_LITE_dcolPivHouseholderQr(const long m, const long n,
                                      double* A, char A_layout,
                                      double* b, double* x);
void CEIGEN_LITE_shouseholderQr(const long m, const long n,
                                float* A, char A_layout,
                                float* b, float* x);
void CEIGEN_LITE_dhouseholderQr(const long m, const long n,
                                double* A, char A_layout,
                                double* b, double* x);

int CEIGEN_LITE_srank(const long m, const long n,
                      float* A, char A_layout,
                      float threshold);
int CEIGEN_LITE_drank(const long m, const long n,
                      double* A, char A_layout,
                      double threshold);
float CEIGEN_LITE_snorm2v(const long m, float* A);
double CEIGEN_LITE_dnorm2v(const long m, double* A);
float CEIGEN_LITE_snorm2m(const long m, const long n, float* A, char A_layout);
double CEIGEN_LITE_dnorm2m(const long m, const long n, double* A, char A_layout);

float CEIGEN_LITE_sdeterminant(const long m, float* A, char A_layout);
double CEIGEN_LITE_ddeterminant(const long m, double* A, char A_layout);

void CEIGEN_LITE_sinverse(const long m, float* A, char A_layout,
                          float* Ainv, char Ainv_layout);
void CEIGEN_LITE_dinverse(const long m, double* A, char A_layout,
                          double* Ainv, char Ainv_layout);
void CEIGEN_LITE_spinverse(const long m, const long n,
                           float* A, char A_layout,
                           float* Ainv, char Ainv_layout);
void CEIGEN_LITE_dpinverse(const long m, const long n,
                           double* A, char A_layout,
                           double* Ainv, char Ainv_layout);
void CEIGEN_LITE_scholesky(const long n,
                           float* A, char A_layout,
                           float* L, char L_layout);
void CEIGEN_LITE_dcholesky(const long m, const long n,
                           double* A, char A_layout,
                           double* L, char L_layout);

void CEIGEN_LITE_sqr(const long m, const long n,
                     float* a, char a_layout,
                     float* q, char q_layout,
                     float* r, char r_layout);
void CEIGEN_LITE_dqr(const long m, const long n,
                     double* a, char a_layout,
                     double* q, char q_layout,
                     double* r, char r_layout);

void CEIGEN_LITE_slu(const long m, const long n,
                     char layout,
                     float* a,  float* p,
                     float* lu, float* q);
void CEIGEN_LITE_dlu(const long m, const long n,
                     char layout,
                     double* a,  double* p,
                     double* lu, double* q);

void CEIGEN_LITE_ssvd(const long m, const long n,
                      float* a, char a_layout,
                      float* u, char u_layout,
                      float* v, char v_layout,
                      float* s);
void CEIGEN_LITE_dsvd(const long m, const long n,
                      double* a, char a_layout,
                      double* u, char u_layout,
                      double* v, char v_layout,
                      double* s);

void CEIGEN_LITE_seigvals(const long n, float* a, char a_layout, float _Complex* eigvals);
void CEIGEN_LITE_deigvals(const long n, double* a, char a_layout, double _Complex* eigvals);

void CEIGEN_LITE_seigvecs(const long n, float* a, char a_layout,
                          float _Complex* eigvals,
                          float _Complex* eigvecs, char ev_layout);
void CEIGEN_LITE_deigvecs(const long n, double* a, char a_layout,
                          double _Complex* eigvals,
                          double _Complex* eigvecs, char ev_layout);

void* CEIGEN_LITE_seed_r(unsigned long seed); // Returns the seed object
void CEIGEN_LITE_rng_destroy(void *rng);  // Deletes the seed object
void CEIGEN_LITE_seed(const unsigned long n);


#define rand_fn(name, type, ...)                                \
  void CEIGEN_LITE_##name(const long len, type* x, __VA_ARGS__);	\
  void CEIGEN_LITE_##name##_r(void* rng, const long len, type* x, __VA_ARGS__);

/* ---------- Real-valued distributions (float/double) ---------- */

rand_fn(sbeta, float,  float a, float b);
rand_fn(dbeta, double, double a, double b);

rand_fn(scauchy, float,  float loc, float scale);
rand_fn(dcauchy, double, double loc, double scale);

rand_fn(schiSquared, float,  float ndof);
rand_fn(dchiSquared, double, double ndof);

rand_fn(sexponential, float,  float lambda);
rand_fn(dexponential, double, double lambda);

rand_fn(sextremeValue, float,  float loc, float scale);
rand_fn(dextremeValue, double, double loc, double scale);

rand_fn(sfisherF, float,  float m, float n);
rand_fn(dfisherF, double, double m, double n);

rand_fn(sgamma, float,  float alpha, float beta);
rand_fn(dgamma, double, double alpha, double beta);

rand_fn(slognormal, float,  float mean, float stdev);
rand_fn(dlognormal, double, double mean, double stdev);

rand_fn(snormal, float,  float mean, float stdev);
rand_fn(dnormal, double, double mean, double stdev);

rand_fn(sstudentT, float,  float ndof);
rand_fn(dstudentT, double, double ndof);

rand_fn(suniformReal, float,  float min, float max);
rand_fn(duniformReal, double, double min, double max);

rand_fn(sweibull, float,  float shape, float scale);
rand_fn(dweibull, double, double shape, double scale);

/* ---------- Integer distributions ---------- */
//   Only 32 bit integers are supported, see: https://github.com/bab2min/EigenRand/issues/58

rand_fn(i32bernoulli, int32_t, double p);
/* rand_fn(i64bernoulli, int64_t, double p); */

rand_fn(i32binomial, int32_t, int32_t trials, double p);
/* rand_fn(i64binomial, int64_t, int64_t trials, double p); */

rand_fn(i32geometric, int32_t, double p);
/* rand_fn(i64geometric, int64_t, double p); */

rand_fn(i32negativeBinomial, int32_t, int32_t trials, double p);
/* rand_fn(i64negativeBinomial, int64_t, int64_t trials, double p); */

rand_fn(i32poisson, int32_t, double mean);
/* rand_fn(i64poisson, int64_t, double mean); */

rand_fn(i32uniformInt, int32_t, int32_t min, int32_t max);
/* rand_fn(i64uniformInt, int64_t, int64_t min, int64_t max); */
