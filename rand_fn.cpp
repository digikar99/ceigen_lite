
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

  void CEIGEN_LITE_snormal_r(CEIGEN_LITE_rng* urng64,
			     const long n, float *x, float mean, float stdev){
    Map<MatrixX<float>> eigmat(x, n, 1);
    eigmat = Rand::normal<MatrixX<float>>(n, 1, urng64->engine, mean, stdev);
  }
  void CEIGEN_LITE_snormal(const long n, float *x, float mean, float stdev){
    CEIGEN_LITE_snormal_r(&global_urng64, n, x, mean, stdev);
  }

  void CEIGEN_LITE_dnormal_r(CEIGEN_LITE_rng* urng64,
			     const long n, double *x, double mean, double stdev){
    Map<MatrixX<double>> eigmat(x, n, 1);
    eigmat = Rand::normal<MatrixX<double>>(n, 1, urng64->engine, mean, stdev);
  }
  void CEIGEN_LITE_dnormal(const long n, double *x, double mean, double stdev){
    CEIGEN_LITE_dnormal_r(&global_urng64, n, x, mean, stdev);
  }


  void CEIGEN_LITE_sbeta_r(CEIGEN_LITE_rng* urng64,
			   const long n, float *x, float a, float b){
    Map<MatrixX<float>> eigmat(x, n, 1);
    eigmat = Rand::beta<MatrixX<float>>(n, 1, urng64->engine, a, b);
    return;
  }
  void CEIGEN_LITE_dbeta_r(CEIGEN_LITE_rng* urng64,
			   const long n, double *x, double a, double b){
    Map<MatrixX<double>> eigmat(x, n, 1);
    eigmat = Rand::beta<MatrixX<double>>(n, 1, urng64->engine, a, b);
    return;
  }

  void CEIGEN_LITE_schiSquared_r(CEIGEN_LITE_rng* urng64,
				 const long n, float *x, float ndof){
    Map<MatrixX<float>> eigmat(x, n, 1);
    eigmat = Rand::chiSquared<MatrixX<float>>(n, 1, urng64->engine, ndof);
    return;
  }
  void CEIGEN_LITE_dchiSquared_r(CEIGEN_LITE_rng* urng64,
				 const long n, double *x, double ndof){
    Map<MatrixX<double>> eigmat(x, n, 1);
    eigmat = Rand::chiSquared<MatrixX<double>>(n, 1, urng64->engine, ndof);
    return;
  }
};
