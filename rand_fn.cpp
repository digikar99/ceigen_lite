

Rand::P8_mt19937_64 urng64{ 42 };
Rand::P8_mt19937_64_32 urng32{ 42 };

extern "C"{
  void CEIGEN_LITE_snormal(const long n, float *x, float mean, float stdev){
    Map<MatrixX<float>> eigmat(x, n, 1);
    eigmat = Rand::normal<MatrixX<float>>(n, 1, urng32, mean, stdev);
    return;
  }
  void CEIGEN_LITE_dnormal(const long n, double *x, double mean, double stdev){
    Map<MatrixX<double>> eigmat(x, n, 1);
    eigmat = Rand::normal<MatrixX<double>>(n, 1, urng64, mean, stdev);
    return;
  }

  void CEIGEN_LITE_sbeta(const long n, float *x, float a, float b){
    Map<MatrixX<float>> eigmat(x, n, 1);
    eigmat = Rand::beta<MatrixX<float>>(n, 1, urng32, a, b);
    return;
  }
  void CEIGEN_LITE_dbeta(const long n, double *x, double a, double b){
    Map<MatrixX<double>> eigmat(x, n, 1);
    eigmat = Rand::beta<MatrixX<double>>(n, 1, urng64, a, b);
    return;
  }

  void CEIGEN_LITE_schiSquared(const long n, float *x, float ndof){
    Map<MatrixX<float>> eigmat(x, n, 1);
    eigmat = Rand::chiSquared<MatrixX<float>>(n, 1, urng32, ndof);
    return;
  }
  void CEIGEN_LITE_dchiSquared(const long n, double *x, double ndof){
    Map<MatrixX<double>> eigmat(x, n, 1);
    eigmat = Rand::chiSquared<MatrixX<double>>(n, 1, urng64, ndof);
    return;
  }
};
