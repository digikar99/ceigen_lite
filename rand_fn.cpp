
// Strips the outer parens off a tuple: CEIGEN_LITE_ARGS (a, b) -> a, b
#define CEIGEN_LITE_ARGS(...) __VA_ARGS__

#define CEIGEN_LITE_DEFINE_DIST(CTYPE, SUF, DIST, PARAMS, ARGVALS)	\
  void CEIGEN_LITE_##SUF##DIST##_r(CEIGEN_LITE_rng *urng64, const long len, \
				   CTYPE *x, CEIGEN_LITE_ARGS PARAMS) {	\
    Map<MatrixX<CTYPE>> eigmat(x, len, 1);				\
    eigmat = Rand::DIST<MatrixX<CTYPE>>(len, 1, urng64->engine,		\
					CEIGEN_LITE_ARGS ARGVALS);	\
  }									\
  void CEIGEN_LITE_##SUF##DIST(const long len, CTYPE *x,		\
			       CEIGEN_LITE_ARGS PARAMS) {		\
    CEIGEN_LITE_##SUF##DIST##_r(&global_urng64, len, x,			\
				CEIGEN_LITE_ARGS ARGVALS);		\
  }

extern "C"{
  struct CEIGEN_LITE_rng { Rand::Vmt19937_64 engine; };
  CEIGEN_LITE_rng global_urng64{Rand::Vmt19937_64{42}};

  CEIGEN_LITE_rng* CEIGEN_LITE_seed_r(unsigned long seed) {
    return new CEIGEN_LITE_rng{Rand::Vmt19937_64(seed)};
  }
  void CEIGEN_LITE_rng_destroy(CEIGEN_LITE_rng *rng) { delete rng; }
  void CEIGEN_LITE_seed(const unsigned long seed){
    global_urng64.engine = Rand::Vmt19937_64(seed);
  }

  CEIGEN_LITE_DEFINE_DIST(float,  s, beta,   (float a, float b), (a, b));
  CEIGEN_LITE_DEFINE_DIST(double, d, beta,   (double a, double b), (a, b));

  CEIGEN_LITE_DEFINE_DIST(float,  s, cauchy, (float loc = 0, float scale = 1), (loc, scale));
  CEIGEN_LITE_DEFINE_DIST(double, d, cauchy, (double loc = 0, double scale = 1), (loc, scale));

  CEIGEN_LITE_DEFINE_DIST(float,  s, chiSquared, (float ndof), (ndof));
  CEIGEN_LITE_DEFINE_DIST(double, d, chiSquared, (double ndof), (ndof));

  CEIGEN_LITE_DEFINE_DIST(float,  s, exponential, (float lambda), (lambda));
  CEIGEN_LITE_DEFINE_DIST(double, d, exponential, (double lambda), (lambda));

  CEIGEN_LITE_DEFINE_DIST(float,  s, extremeValue, (float loc = 0, float scale = 1), (loc, scale));
  CEIGEN_LITE_DEFINE_DIST(double, d, extremeValue, (double loc = 0, double scale = 1), (loc, scale));

  CEIGEN_LITE_DEFINE_DIST(float,  s, fisherF, (float m, float n), (m, n));
  CEIGEN_LITE_DEFINE_DIST(double, d, fisherF, (double m, double n), (m, n));

  CEIGEN_LITE_DEFINE_DIST(float,  s, gamma, (float alpha = 1, float beta = 1), (alpha, beta));
  CEIGEN_LITE_DEFINE_DIST(double, d, gamma, (double alpha = 1, double beta = 1), (alpha, beta));

  CEIGEN_LITE_DEFINE_DIST(float,  s, lognormal, (float mean = 0, float stdev = 1), (mean, stdev));
  CEIGEN_LITE_DEFINE_DIST(double, d, lognormal, (double mean = 0, double stdev = 1), (mean, stdev));

  CEIGEN_LITE_DEFINE_DIST(float,  s, normal, (float mean, float stdev), (mean, stdev));
  CEIGEN_LITE_DEFINE_DIST(double, d, normal, (double mean, double stdev), (mean, stdev));

  CEIGEN_LITE_DEFINE_DIST(float,  s, studentT, (float ndof = 1), (ndof));
  CEIGEN_LITE_DEFINE_DIST(double, d, studentT, (double ndof = 1), (ndof));

  CEIGEN_LITE_DEFINE_DIST(float,  s, uniformReal, (float min = 0, float max = 1), (min, max));
  CEIGEN_LITE_DEFINE_DIST(double, d, uniformReal, (double min = 0, double max = 1), (min, max));

  CEIGEN_LITE_DEFINE_DIST(float,  s, weibull, (float shape = 1, float scale = 1), (shape, scale));
  CEIGEN_LITE_DEFINE_DIST(double, d, weibull, (double shape = 1, double scale = 1), (shape, scale));

  // Integer distributions
  //   Only 32 bit integers are supported, see: https://github.com/bab2min/EigenRand/issues/58

  CEIGEN_LITE_DEFINE_DIST(int32_t, i32, bernoulli, (double p = 0.5), (p));
  // CEIGEN_LITE_DEFINE_DIST(int64_t, i64, bernoulli, (double p = 0.5), (p));

  CEIGEN_LITE_DEFINE_DIST(int32_t, i32, binomial, (int trials = 1, double p = 0.5), (trials, p));
  // CEIGEN_LITE_DEFINE_DIST(int64_t, i64, binomial, (int trials = 1, double p = 0.5), (trials, p));

  CEIGEN_LITE_DEFINE_DIST(int32_t, i32, geometric, (double p = 0.5), (p));
  // CEIGEN_LITE_DEFINE_DIST(int64_t, i64, geometric, (double p = 0.5), (p));

  CEIGEN_LITE_DEFINE_DIST(int32_t, i32, negativeBinomial, (int32_t trials = 1, double p = 0.5), (trials, p));
  // CEIGEN_LITE_DEFINE_DIST(int64_t, i64, negativeBinomial, (int64_t trials = 1, double p = 0.5), (trials, p));

  CEIGEN_LITE_DEFINE_DIST(int32_t, i32, poisson, (double mean = 1), (mean));
  // CEIGEN_LITE_DEFINE_DIST(int64_t, i64, poisson, (double mean = 1), (mean));

  CEIGEN_LITE_DEFINE_DIST(int32_t, i32, uniformInt, (int32_t min, int32_t max), (min, max));
  // CEIGEN_LITE_DEFINE_DIST(int64_t, i64, uniformInt, (int64_t min = 0, int64_t max = 0), (min, max));
}
