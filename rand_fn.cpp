
// Strips the outer parens off a tuple: CEIGEN_LITE_ARGS (a, b) -> a, b
#define CEIGEN_LITE_ARGS(...) __VA_ARGS__

#define CEIGEN_LITE_DEFINE_DIST(CTYPE, SUF, DIST, PARAMS, ARGVALS)	\
  void CEIGEN_LITE_##SUF##DIST##_r(CEIGEN_LITE_rng *urng64, const long n, \
				   CTYPE *x, CEIGEN_LITE_ARGS PARAMS) {	\
    Map<MatrixX<CTYPE>> eigmat(x, n, 1);				\
    eigmat = Rand::DIST<MatrixX<CTYPE>>(n, 1, urng64->engine,		\
					CEIGEN_LITE_ARGS ARGVALS);	\
  }									\
  void CEIGEN_LITE_##SUF##DIST(const long n, CTYPE *x,			\
			       CEIGEN_LITE_ARGS PARAMS) {		\
    CEIGEN_LITE_##SUF##DIST##_r(&global_urng64, n, x,			\
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

  CEIGEN_LITE_DEFINE_DIST(float,  s, normal, (float mean, float stdev), (mean, stdev));
  CEIGEN_LITE_DEFINE_DIST(double, d, normal, (double mean, double stdev), (mean, stdev));
  CEIGEN_LITE_DEFINE_DIST(float,  s, beta,   (float a, float b), (a, b));
  CEIGEN_LITE_DEFINE_DIST(double, d, beta,   (double a, double b), (a, b));
  CEIGEN_LITE_DEFINE_DIST(float,  s, chiSquared, (float ndof), (ndof));
  CEIGEN_LITE_DEFINE_DIST(double, d, chiSquared, (double ndof), (ndof));
}
