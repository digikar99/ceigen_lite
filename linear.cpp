
template <typename T, int options>
using MatrixRC = Matrix<T, Dynamic, Dynamic, options>;

extern "C"{

  // matrix multiplication
  void CEIGEN_LITE_smatmul(const long m, const long n, const long k,
                           float* a, char a_layout,
                           float* b, char b_layout,
                           float* c, char c_layout){
    if (a_layout == 'C' || a_layout == 'c'){
      Map<MatrixRC<float, ColMajor>> amat(a, m, n);
      if (b_layout == 'C' || b_layout == 'c'){
        Map<MatrixRC<float, ColMajor>> bmat(b, n, k);
        if (c_layout == 'C' || c_layout == 'c'){
          Map<MatrixRC<float, ColMajor>> cmat(c, m, k);
          cmat = amat*bmat;
        }else if (c_layout == 'R' || c_layout == 'r'){
          Map<MatrixRC<float, RowMajor>> cmat(c, m, k);
          cmat = amat*bmat;
        }
      }else if (b_layout == 'R' || b_layout == 'r'){
        Map<MatrixRC<float, RowMajor>> bmat(b, n, k);
        if (c_layout == 'C' || c_layout == 'c'){
          Map<MatrixRC<float, ColMajor>> cmat(c, m, k);
          cmat = amat*bmat;
        }else if (c_layout == 'R' || c_layout == 'r'){
          Map<MatrixRC<float, RowMajor>> cmat(c, m, k);
          cmat = amat*bmat;
        }
      }
    }else if (a_layout == 'R' || a_layout == 'r'){
      Map<MatrixRC<float, RowMajor>> amat(a, m, n);
      if (b_layout == 'C' || b_layout == 'c'){
        Map<MatrixRC<float, ColMajor>> bmat(b, n, k);
        if (c_layout == 'C' || c_layout == 'c'){
          Map<MatrixRC<float, ColMajor>> cmat(c, m, k);
          cmat = amat*bmat;
        }else if (c_layout == 'R' || c_layout == 'r'){
          Map<MatrixRC<float, RowMajor>> cmat(c, m, k);
          cmat = amat*bmat;
        }
      }else if (b_layout == 'R' || b_layout == 'r'){
        Map<MatrixRC<float, RowMajor>> bmat(b, n, k);
        if (c_layout == 'C' || c_layout == 'c'){
          Map<MatrixRC<float, ColMajor>> cmat(c, m, k);
          cmat = amat*bmat;
        }else if (c_layout == 'R' || c_layout == 'r'){
          Map<MatrixRC<float, RowMajor>> cmat(c, m, k);
          cmat = amat*bmat;
        }
      }
    }
  }
  void CEIGEN_LITE_dmatmul(const long m, const long n, const long k,
                           double* a, char a_layout,
                           double* b, char b_layout,
                           double* c, char c_layout){
    if (a_layout == 'C' || a_layout == 'c'){
      Map<MatrixRC<double, ColMajor>> amat(a, m, n);
      if (b_layout == 'C' || b_layout == 'c'){
        Map<MatrixRC<double, ColMajor>> bmat(b, n, k);
        if (c_layout == 'C' || c_layout == 'c'){
          Map<MatrixRC<double, ColMajor>> cmat(c, m, k);
          cmat = amat*bmat;
        }else if (c_layout == 'R' || c_layout == 'r'){
          Map<MatrixRC<double, RowMajor>> cmat(c, m, k);
          cmat = amat*bmat;
        }
      }else if (b_layout == 'R' || b_layout == 'r'){
        Map<MatrixRC<double, RowMajor>> bmat(b, n, k);
        if (c_layout == 'C' || c_layout == 'c'){
          Map<MatrixRC<double, ColMajor>> cmat(c, m, k);
          cmat = amat*bmat;
        }else if (c_layout == 'R' || c_layout == 'r'){
          Map<MatrixRC<double, RowMajor>> cmat(c, m, k);
          cmat = amat*bmat;
        }
      }
    }else if (a_layout == 'R' || a_layout == 'r'){
      Map<MatrixRC<double, RowMajor>> amat(a, m, n);
      if (b_layout == 'C' || b_layout == 'c'){
        Map<MatrixRC<double, ColMajor>> bmat(b, n, k);
        if (c_layout == 'C' || c_layout == 'c'){
          Map<MatrixRC<double, ColMajor>> cmat(c, m, k);
          cmat = amat*bmat;
        }else if (c_layout == 'R' || c_layout == 'r'){
          Map<MatrixRC<double, RowMajor>> cmat(c, m, k);
          cmat = amat*bmat;
        }
      }else if (b_layout == 'R' || b_layout == 'r'){
        Map<MatrixRC<double, RowMajor>> bmat(b, n, k);
        if (c_layout == 'C' || c_layout == 'c'){
          Map<MatrixRC<double, ColMajor>> cmat(c, m, k);
          cmat = amat*bmat;
        }else if (c_layout == 'R' || c_layout == 'r'){
          Map<MatrixRC<double, RowMajor>> cmat(c, m, k);
          cmat = amat*bmat;
        }
      }
    }
  }

  // Equation Solving: A is a matrix with m rows and n cols
  void CEIGEN_LITE_spartialPivLu(const long m, const long n,
                                 float* A, char A_layout,
                                 float* b, float* x){
    Map<VectorX<float>> bvec(b, n);
    Map<VectorX<float>> xvec(x, n);
    if (A_layout == 'C' || A_layout == 'c'){
      Map<MatrixRC<float, ColMajor>> amat(A, m, n);
      xvec = amat.partialPivLu().solve(bvec);
    }else if (A_layout == 'R' || A_layout == 'r'){
      Map<MatrixRC<float, RowMajor>> amat(A, m, n);
      xvec = amat.partialPivLu().solve(bvec);
    }
  }
  void CEIGEN_LITE_dpartialPivLu(const long m, const long n,
                                 double* A, char A_layout,
                                 double* b, double* x){
    Map<VectorX<double>> bvec(b, n);
    Map<VectorX<double>> xvec(x, n);
    if (A_layout == 'C' || A_layout == 'c'){
      Map<MatrixRC<double, ColMajor>> amat(A, m, n);
      xvec = amat.partialPivLu().solve(bvec);
    }else if (A_layout == 'R' || A_layout == 'r'){
      Map<MatrixRC<double, RowMajor>> amat(A, m, n);
      xvec = amat.partialPivLu().solve(bvec);
    }
  }

  void CEIGEN_LITE_scolPivHouseholderQr(const long m, const long n,
                                        float* A, char A_layout,
                                        float* b, float* x){
    Map<VectorX<float>> bvec(b, n);
    Map<VectorX<float>> xvec(x, n);
    if (A_layout == 'C' || A_layout == 'c'){
      Map<MatrixRC<float, ColMajor>> amat(A, m, n);
      xvec = amat.colPivHouseholderQr().solve(bvec);
    }else if (A_layout == 'R' || A_layout == 'r'){
      Map<MatrixRC<float, RowMajor>> amat(A, m, n);
      xvec = amat.colPivHouseholderQr().solve(bvec);
    }
  }
  void CEIGEN_LITE_dcolPivHouseholderQr(const long m, const long n,
                                        double* A, char A_layout,
                                        double* b, double* x){
    Map<VectorX<double>> bvec(b, n);
    Map<VectorX<double>> xvec(x, n);
    if (A_layout == 'C' || A_layout == 'c'){
      Map<MatrixRC<double, ColMajor>> amat(A, m, n);
      xvec = amat.colPivHouseholderQr().solve(bvec);
    }else if (A_layout == 'R' || A_layout == 'r'){
      Map<MatrixRC<double, RowMajor>> amat(A, m, n);
      xvec = amat.colPivHouseholderQr().solve(bvec);
    }
  }

  void CEIGEN_LITE_shouseholderQr(const long m, const long n,
                                  float* A, char A_layout,
                                  float* b, float* x){
    Map<VectorX<float>> bvec(b, n);
    Map<VectorX<float>> xvec(x, n);
    if (A_layout == 'C' || A_layout == 'c'){
      Map<MatrixRC<float, ColMajor>> amat(A, m, n);
      xvec = amat.householderQr().solve(bvec);
    }else if (A_layout == 'R' || A_layout == 'r'){
      Map<MatrixRC<float, RowMajor>> amat(A, m, n);
      xvec = amat.householderQr().solve(bvec);
    }
  }
  void CEIGEN_LITE_dhouseholderQr(const long m, const long n,
                                  double* A, char A_layout,
                                  double* b, double* x){
    Map<VectorX<double>> bvec(b, n);
    Map<VectorX<double>> xvec(x, n);
    if (A_layout == 'C' || A_layout == 'c'){
      Map<MatrixRC<double, ColMajor>> amat(A, m, n);
      xvec = amat.householderQr().solve(bvec);
    }else if (A_layout == 'R' || A_layout == 'r'){
      Map<MatrixRC<double, RowMajor>> amat(A, m, n);
      xvec = amat.householderQr().solve(bvec);
    }
  }

  // Computing Determinants
  float CEIGEN_LITE_sdeterminant(const long m, float* A, char A_layout){
    if (A_layout == 'C' || A_layout == 'c'){
      Map<MatrixRC<float, ColMajor>> amat(A, m, m);
      return amat.determinant();
    }else if (A_layout == 'R' || A_layout == 'r'){
      Map<MatrixRC<float, RowMajor>> amat(A, m, m);
      return amat.determinant();
    }else return -1;
  }
  double CEIGEN_LITE_ddeterminant(const long m, double* A, char A_layout){
    if (A_layout == 'C' || A_layout == 'c'){
      Map<MatrixRC<double, ColMajor>> amat(A, m, m);
      return amat.determinant();
    }else if (A_layout == 'R' || A_layout == 'r'){
      Map<MatrixRC<double, RowMajor>> amat(A, m, m);
      return amat.determinant();
    }else return -1;
  }

  // Computing ranks
  int CEIGEN_LITE_srank(const long m, const long n,
                        float* A, char A_layout,
                        float threshold=0){
    if (A_layout == 'C' || A_layout == 'c'){
      Map<MatrixRC<float, ColMajor>> amat(A, m, n);
      auto qr = amat.colPivHouseholderQr();
      if (threshold != 0) qr.setThreshold(threshold);
      return qr.rank();
    }else if (A_layout == 'R' || A_layout == 'r'){
      Map<MatrixRC<float, RowMajor>> amat(A, m, n);
      auto qr = amat.colPivHouseholderQr();
      if (threshold != 0) qr.setThreshold(threshold);
      return qr.rank();
    }else return -1;
  }
  int CEIGEN_LITE_drank(const long m, const long n,
                        double* A, char A_layout,
                        double threshold=0){
    if (A_layout == 'C' || A_layout == 'c'){
      Map<MatrixRC<double, ColMajor>> amat(A, m, n);
      auto qr = amat.colPivHouseholderQr();
      if (threshold != 0) qr.setThreshold(threshold);
      return qr.rank();
    }else if (A_layout == 'R' || A_layout == 'r'){
      Map<MatrixRC<double, RowMajor>> amat(A, m, n);
      auto qr = amat.colPivHouseholderQr();
      if (threshold != 0) qr.setThreshold(threshold);
      return qr.rank();
    }else return -1;
  }

  // Computing L2 norms
  float CEIGEN_LITE_snorm2v(const long m, float* A){
    Map<VectorX<float>> avec(A, m);
    return avec.norm();
  }
  double CEIGEN_LITE_dnorm2v(const long m, double* A){
    Map<VectorX<double>> avec(A, m);
    return avec.norm();
  }
  float CEIGEN_LITE_snorm2m(const long m, const long n, float* A, char A_layout){
    if (A_layout == 'C' || A_layout == 'c'){
      Map<MatrixRC<float, ColMajor>> amat(A, m, n);
      return amat.norm();
    }else if (A_layout == 'R' || A_layout == 'r'){
      Map<MatrixRC<float, RowMajor>> amat(A, m, n);
      return amat.norm();
    }else return -1;
  }
  double CEIGEN_LITE_dnorm2m(const long m, const long n, double* A, char A_layout){
    if (A_layout == 'C' || A_layout == 'c'){
      Map<MatrixRC<double, ColMajor>> amat(A, m, n);
      return amat.norm();
    }else if (A_layout == 'R' || A_layout == 'r'){
      Map<MatrixRC<double, RowMajor>> amat(A, m, n);
      return amat.norm();
    }else return -1;
  }

  // Computing inverses
  void CEIGEN_LITE_sinverse(const long m, float* A, char A_layout,
                            float* Ainv, char Ainv_layout){
    if (A_layout == 'C' || A_layout == 'c'){
      Map<MatrixRC<float, ColMajor>> amat(A, m, m);
      if (Ainv_layout == 'C' || Ainv_layout == 'c'){
        Map<MatrixRC<float, ColMajor>> invmat(Ainv, m, m);
        invmat = amat.inverse();
      }else if (Ainv_layout == 'R' || Ainv_layout == 'r'){
        Map<MatrixRC<float, RowMajor>> invmat(Ainv, m, m);
        invmat = amat.inverse();
      }
    }else if (A_layout == 'R' || A_layout == 'r'){
      Map<MatrixRC<float, RowMajor>> amat(A, m, m);
      if (Ainv_layout == 'C' || Ainv_layout == 'c'){
        Map<MatrixRC<float, ColMajor>> invmat(Ainv, m, m);
        invmat = amat.inverse();
      }else if (Ainv_layout == 'R' || Ainv_layout == 'r'){
        Map<MatrixRC<float, RowMajor>> invmat(Ainv, m, m);
        invmat = amat.inverse();
      }
    }
  }
  void CEIGEN_LITE_dinverse(const long m, double* A, char A_layout,
                              double* Ainv, char Ainv_layout){
    if (A_layout == 'C' || A_layout == 'c'){
      Map<MatrixRC<double, ColMajor>> amat(A, m, m);
      if (Ainv_layout == 'C' || Ainv_layout == 'c'){
        Map<MatrixRC<double, ColMajor>> invmat(Ainv, m, m);
        invmat = amat.inverse();
      }else if (Ainv_layout == 'R' || Ainv_layout == 'r'){
        Map<MatrixRC<double, RowMajor>> invmat(Ainv, m, m);
        invmat = amat.inverse();
      }
    }else if (A_layout == 'R' || A_layout == 'r'){
      Map<MatrixRC<double, RowMajor>> amat(A, m, m);
      if (Ainv_layout == 'C' || Ainv_layout == 'c'){
        Map<MatrixRC<double, ColMajor>> invmat(Ainv, m, m);
        invmat = amat.inverse();
      }else if (Ainv_layout == 'R' || Ainv_layout == 'r'){
        Map<MatrixRC<double, RowMajor>> invmat(Ainv, m, m);
        invmat = amat.inverse();
      }
    }
  }

  // A is a matrix with m rows and n cols
  void CEIGEN_LITE_spinverse(const long m, const long n,
                             float* A, char A_layout,
                             float* Ainv, char Ainv_layout){
    if (A_layout == 'C' || A_layout == 'c'){
      Map<MatrixRC<float, ColMajor>> amat(A, m, n);
      if (Ainv_layout == 'C' || Ainv_layout == 'c'){
        Map<MatrixRC<float, ColMajor>> invmat(Ainv, n, m);
        invmat = amat.completeOrthogonalDecomposition().pseudoInverse();
      }else if (Ainv_layout == 'R' || Ainv_layout == 'r'){
        Map<MatrixRC<float, RowMajor>> invmat(Ainv, n, m);
        invmat = amat.completeOrthogonalDecomposition().pseudoInverse();
      }
    }else if (A_layout == 'R' || A_layout == 'r'){
      Map<MatrixRC<float, RowMajor>> amat(A, m, n);
      if (Ainv_layout == 'C' || Ainv_layout == 'c'){
        Map<MatrixRC<float, ColMajor>> invmat(Ainv, n, m);
        invmat = amat.completeOrthogonalDecomposition().pseudoInverse();
      }else if (Ainv_layout == 'R' || Ainv_layout == 'r'){
        Map<MatrixRC<float, RowMajor>> invmat(Ainv, n, m);
        invmat = amat.completeOrthogonalDecomposition().pseudoInverse();
      }
    }
  }
  void CEIGEN_LITE_dpinverse(const long m, const long n,
                             double* A, char A_layout,
                             double* Ainv, char Ainv_layout){
    if (A_layout == 'C' || A_layout == 'c'){
      Map<MatrixRC<double, ColMajor>> amat(A, m, n);
      if (Ainv_layout == 'C' || Ainv_layout == 'c'){
        Map<MatrixRC<double, ColMajor>> invmat(Ainv, n, m);
        invmat = amat.completeOrthogonalDecomposition().pseudoInverse();
      }else if (Ainv_layout == 'R' || Ainv_layout == 'r'){
        Map<MatrixRC<double, RowMajor>> invmat(Ainv, n, m);
        invmat = amat.completeOrthogonalDecomposition().pseudoInverse();
      }
    }else if (A_layout == 'R' || A_layout == 'r'){
      Map<MatrixRC<double, RowMajor>> amat(A, m, n);
      if (Ainv_layout == 'C' || Ainv_layout == 'c'){
        Map<MatrixRC<double, ColMajor>> invmat(Ainv, n, m);
        invmat = amat.completeOrthogonalDecomposition().pseudoInverse();
      }else if (Ainv_layout == 'R' || Ainv_layout == 'r'){
        Map<MatrixRC<double, RowMajor>> invmat(Ainv, n, m);
        invmat = amat.completeOrthogonalDecomposition().pseudoInverse();
      }
    }
  }

  // Computing Choleskey
  void CEIGEN_LITE_scholesky(const long n,
                             float* A, char A_layout,
                             float* L, char L_layout){
    if (A_layout == 'C' || A_layout == 'c'){
      Map<MatrixRC<float, ColMajor>> amat(A, n, n);
      if (L_layout == 'C' || L_layout == 'c'){
        Map<MatrixRC<float, ColMajor>> lmat(L, n, n);
        lmat = amat.llt().matrixL();
      }else if (L_layout == 'R' || L_layout == 'r'){
        Map<MatrixRC<float, RowMajor>> lmat(L, n, n);
        lmat = amat.llt().matrixL();
      }
    }else if (A_layout == 'R' || A_layout == 'r'){
      Map<MatrixRC<float, RowMajor>> amat(A, n, n);
      if (L_layout == 'C' || L_layout == 'c'){
        Map<MatrixRC<float, ColMajor>> lmat(L, n, n);
        lmat = amat.llt().matrixL();
      }else if (L_layout == 'R' || L_layout == 'r'){
        Map<MatrixRC<float, RowMajor>> lmat(L, n, n);
        lmat = amat.llt().matrixL();
      }
    }
  }
  void CEIGEN_LITE_dcholesky(const long m, const long n,
                             double* A, char A_layout,
                             double* L, char L_layout){
    if (A_layout == 'C' || A_layout == 'c'){
      Map<MatrixRC<double, ColMajor>> amat(A, n, n);
      if (L_layout == 'C' || L_layout == 'c'){
        Map<MatrixRC<double, ColMajor>> lmat(L, n, n);
        lmat = amat.llt().matrixL();
      }else if (L_layout == 'R' || L_layout == 'r'){
        Map<MatrixRC<double, RowMajor>> lmat(L, n, n);
        lmat = amat.llt().matrixL();
      }
    }else if (A_layout == 'R' || A_layout == 'r'){
      Map<MatrixRC<double, RowMajor>> amat(A, n, n);
      if (L_layout == 'C' || L_layout == 'c'){
        Map<MatrixRC<double, ColMajor>> lmat(L, n, n);
        lmat = amat.llt().matrixL();
      }else if (L_layout == 'R' || L_layout == 'r'){
        Map<MatrixRC<double, RowMajor>> lmat(L, n, n);
        lmat = amat.llt().matrixL();
      }
    }
  }

  // Computing QR
  void CEIGEN_LITE_sqr(const long m, const long n,
                       float* a, char a_layout,
                       float* q, char q_layout,
                       float* r, char r_layout){
    if (a_layout == 'C' || a_layout == 'c'){
      Map<MatrixRC<float, ColMajor>> amat(a, m, n);
      auto qr = amat.colPivHouseholderQr();
      if (q_layout == 'C' || q_layout == 'c'){
        Map<MatrixRC<float, ColMajor>> qmat(q, m, m);
        qmat = qr.matrixQ();
        if (r_layout == 'C' || r_layout == 'c'){
          Map<MatrixRC<float, ColMajor>> rmat(r, m, n);
          rmat = qr.matrixQR().triangularView<Eigen::Upper>();
        }else if (r_layout == 'R' || r_layout == 'r'){
          Map<MatrixRC<float, RowMajor>> rmat(r, m, n);
          rmat = qr.matrixQR().triangularView<Eigen::Upper>();
        }
      }else if (q_layout == 'R' || q_layout == 'r'){
        Map<MatrixRC<float, RowMajor>> qmat(q, m, m);
        qmat = qr.matrixQ();
        if (r_layout == 'C' || r_layout == 'c'){
          Map<MatrixRC<float, ColMajor>> rmat(r, m, n);
          rmat = qr.matrixQR().triangularView<Eigen::Upper>();
        }else if (r_layout == 'R' || r_layout == 'r'){
          Map<MatrixRC<float, RowMajor>> rmat(r, m, n);
          rmat = qr.matrixQR().triangularView<Eigen::Upper>();
        }
      }
    }else if (a_layout == 'R' || a_layout == 'r'){
      Map<MatrixRC<float, RowMajor>> amat(a, m, n);
      auto qr = amat.colPivHouseholderQr();
      if (q_layout == 'C' || q_layout == 'c'){
        Map<MatrixRC<float, ColMajor>> qmat(q, m, m);
        qmat = qr.matrixQ();
        if (r_layout == 'C' || r_layout == 'c'){
          Map<MatrixRC<float, ColMajor>> rmat(r, m, n);
          rmat = qr.matrixQR().triangularView<Eigen::Upper>();
        }else if (r_layout == 'R' || r_layout == 'r'){
          Map<MatrixRC<float, RowMajor>> rmat(r, m, n);
          rmat = qr.matrixQR().triangularView<Eigen::Upper>();
        }
      }else if (q_layout == 'R' || q_layout == 'r'){
        Map<MatrixRC<float, RowMajor>> qmat(q, m, m);
        qmat = qr.matrixQ();
        if (r_layout == 'C' || r_layout == 'c'){
          Map<MatrixRC<float, ColMajor>> rmat(r, m, n);
          rmat = qr.matrixQR().triangularView<Eigen::Upper>();
        }else if (r_layout == 'R' || r_layout == 'r'){
          Map<MatrixRC<float, RowMajor>> rmat(r, m, n);
          rmat = qr.matrixQR().triangularView<Eigen::Upper>();
        }
      }
    }
  }
  void CEIGEN_LITE_dqr(const long m, const long n,
                       double* a, char a_layout,
                       double* q, char q_layout,
                       double* r, char r_layout){
    if (a_layout == 'C' || a_layout == 'c'){
      Map<MatrixRC<double, ColMajor>> amat(a, m, n);
      auto qr = amat.colPivHouseholderQr();
      if (q_layout == 'C' || q_layout == 'c'){
        Map<MatrixRC<double, ColMajor>> qmat(q, m, m);
        qmat = qr.matrixQ();
        if (r_layout == 'C' || r_layout == 'c'){
          Map<MatrixRC<double, ColMajor>> rmat(r, m, n);
          rmat = qr.matrixQR().triangularView<Eigen::Upper>();
        }else if (r_layout == 'R' || r_layout == 'r'){
          Map<MatrixRC<double, RowMajor>> rmat(r, m, n);
          rmat = qr.matrixQR().triangularView<Eigen::Upper>();
        }
      }else if (q_layout == 'R' || q_layout == 'r'){
        Map<MatrixRC<double, RowMajor>> qmat(q, m, m);
        qmat = qr.matrixQ();
        if (r_layout == 'C' || r_layout == 'c'){
          Map<MatrixRC<double, ColMajor>> rmat(r, m, n);
          rmat = qr.matrixQR().triangularView<Eigen::Upper>();
        }else if (r_layout == 'R' || r_layout == 'r'){
          Map<MatrixRC<double, RowMajor>> rmat(r, m, n);
          rmat = qr.matrixQR().triangularView<Eigen::Upper>();
        }
      }
    }else if (a_layout == 'R' || a_layout == 'r'){
      Map<MatrixRC<double, RowMajor>> amat(a, m, n);
      auto qr = amat.colPivHouseholderQr();
      if (q_layout == 'C' || q_layout == 'c'){
        Map<MatrixRC<double, ColMajor>> qmat(q, m, m);
        qmat = qr.matrixQ();
        if (r_layout == 'C' || r_layout == 'c'){
          Map<MatrixRC<double, ColMajor>> rmat(r, m, n);
          rmat = qr.matrixQR().triangularView<Eigen::Upper>();
        }else if (r_layout == 'R' || r_layout == 'r'){
          Map<MatrixRC<double, RowMajor>> rmat(r, m, n);
          rmat = qr.matrixQR().triangularView<Eigen::Upper>();
        }
      }else if (q_layout == 'R' || q_layout == 'r'){
        Map<MatrixRC<double, RowMajor>> qmat(q, m, m);
        qmat = qr.matrixQ();
        if (r_layout == 'C' || r_layout == 'c'){
          Map<MatrixRC<double, ColMajor>> rmat(r, m, n);
          rmat = qr.matrixQR().triangularView<Eigen::Upper>();
        }else if (r_layout == 'R' || r_layout == 'r'){
          Map<MatrixRC<double, RowMajor>> rmat(r, m, n);
          rmat = qr.matrixQR().triangularView<Eigen::Upper>();
        }
      }
    }
  }

  // Compute LU decomposition
  void CEIGEN_LITE_slu(const long m, const long n,
                       char layout,
                       float* a,  float* p,
                       float* lu, float* q){

    if (layout == 'C' || layout == 'c'){
      Map<MatrixRC<float, ColMajor>> amat(a, m, n);
      Map<MatrixRC<float, ColMajor>> pmat(p, m, m);
      Map<MatrixRC<float, ColMajor>> lumat(lu, m, n);
      Map<MatrixRC<float, ColMajor>> qmat(q, n, n);
      auto lu = amat.fullPivLu();
      lumat = lu.matrixLU();
      pmat = lu.permutationP();
      qmat = lu.permutationQ();
    }else if (layout == 'R' || layout == 'r'){
      Map<MatrixRC<float, RowMajor>> amat(a, m, n);
      Map<MatrixRC<float, RowMajor>> pmat(p, m, m);
      Map<MatrixRC<float, RowMajor>> lumat(lu, m, n);
      Map<MatrixRC<float, RowMajor>> qmat(q, n, n);
      auto lu = amat.fullPivLu();
      lumat = lu.matrixLU();
      pmat = lu.permutationP();
      qmat = lu.permutationQ();
    }
  }
  void CEIGEN_LITE_dlu(const long m, const long n,
                       char layout,
                       double* a,  double* p,
                       double* lu, double* q){

    if (layout == 'C' || layout == 'c'){
      Map<MatrixRC<double, ColMajor>> amat(a, m, n);
      Map<MatrixRC<double, ColMajor>> pmat(p, m, m);
      Map<MatrixRC<double, ColMajor>> lumat(lu, m, n);
      Map<MatrixRC<double, ColMajor>> qmat(q, n, n);
      auto lu = amat.fullPivLu();
      lumat = lu.matrixLU();
      pmat = lu.permutationP();
      qmat = lu.permutationQ();
    }else if (layout == 'R' || layout == 'r'){
      Map<MatrixRC<double, RowMajor>> amat(a, m, n);
      Map<MatrixRC<double, RowMajor>> pmat(p, m, m);
      Map<MatrixRC<double, RowMajor>> lumat(lu, m, n);
      Map<MatrixRC<double, RowMajor>> qmat(q, n, n);
      auto lu = amat.fullPivLu();
      lumat = lu.matrixLU();
      pmat = lu.permutationP();
      qmat = lu.permutationQ();
    }
  }

  // Computing SVD
  void CEIGEN_LITE_ssvd(const long m, const long n,
                        float* a, char a_layout,
                        float* u, char u_layout,
                        float* v, char v_layout,
                        float* s){

    const long len = (m<n ? m : n);
    Map<VectorX<float>> svec(s, len);
    if (a_layout == 'C' || a_layout == 'c'){
      Map<MatrixRC<float, ColMajor>> amat(a, m, n);
      auto svd = amat.bdcSvd(ComputeFullU | ComputeFullV);
      svec = svd.singularValues();
      if (u_layout == 'C' || u_layout == 'c'){
        Map<MatrixRC<float, ColMajor>> umat(u, m, m);
        umat = svd.matrixU();
        if (v_layout == 'C' || v_layout == 'c'){
          Map<MatrixRC<float, ColMajor>> vmat(v, n, n);
          vmat = svd.matrixV();
        }else if (v_layout == 'R' || v_layout == 'r'){
          Map<MatrixRC<float, RowMajor>> vmat(v, n, n);
          vmat = svd.matrixV();
        }
      }else if (u_layout == 'R' || u_layout == 'r'){
        Map<MatrixRC<float, RowMajor>> umat(u, m, m);
        umat = svd.matrixU();
        if (v_layout == 'C' || v_layout == 'c'){
          Map<MatrixRC<float, ColMajor>> vmat(v, n, n);
          vmat = svd.matrixV();
        }else if (v_layout == 'R' || v_layout == 'r'){
          Map<MatrixRC<float, RowMajor>> vmat(v, n, n);
          vmat = svd.matrixV();
        }
      }
    }else if (a_layout == 'R' || a_layout == 'r'){
      Map<MatrixRC<float, RowMajor>> amat(a, m, n);
      auto svd = amat.bdcSvd(ComputeFullU | ComputeFullV);
      svec = svd.singularValues();
      if (u_layout == 'C' || u_layout == 'c'){
        Map<MatrixRC<float, ColMajor>> umat(u, m, m);
        umat = svd.matrixU();
        if (v_layout == 'C' || v_layout == 'c'){
          Map<MatrixRC<float, ColMajor>> vmat(v, n, n);
          vmat = svd.matrixV();
        }else if (v_layout == 'R' || v_layout == 'r'){
          Map<MatrixRC<float, RowMajor>> vmat(v, n, n);
          vmat = svd.matrixV();
        }
      }else if (u_layout == 'R' || u_layout == 'r'){
        Map<MatrixRC<float, RowMajor>> umat(u, m, m);
        umat = svd.matrixU();
        if (v_layout == 'C' || v_layout == 'c'){
          Map<MatrixRC<float, ColMajor>> vmat(v, n, n);
          vmat = svd.matrixV();
        }else if (v_layout == 'R' || v_layout == 'r'){
          Map<MatrixRC<float, RowMajor>> vmat(v, n, n);
          vmat = svd.matrixV();
        }
      }
    }
  }
  void CEIGEN_LITE_dsvd(const long m, const long n,
                        double* a, char a_layout,
                        double* u, char u_layout,
                        double* v, char v_layout,
                        double* s){

    const long len = (m<n ? m : n);
    Map<VectorX<double>> svec(s, len);
    if (a_layout == 'C' || a_layout == 'c'){
      Map<MatrixRC<double, ColMajor>> amat(a, m, n);
      auto svd = amat.bdcSvd(ComputeFullU | ComputeFullV);
      svec = svd.singularValues();
      if (u_layout == 'C' || u_layout == 'c'){
        Map<MatrixRC<double, ColMajor>> umat(u, m, m);
        umat = svd.matrixU();
        if (v_layout == 'C' || v_layout == 'c'){
          Map<MatrixRC<double, ColMajor>> vmat(v, n, n);
          vmat = svd.matrixV();
        }else if (v_layout == 'R' || v_layout == 'r'){
          Map<MatrixRC<double, RowMajor>> vmat(v, n, n);
          vmat = svd.matrixV();
        }
      }else if (u_layout == 'R' || u_layout == 'r'){
        Map<MatrixRC<double, RowMajor>> umat(u, m, m);
        umat = svd.matrixU();
        if (v_layout == 'C' || v_layout == 'c'){
          Map<MatrixRC<double, ColMajor>> vmat(v, n, n);
          vmat = svd.matrixV();
        }else if (v_layout == 'R' || v_layout == 'r'){
          Map<MatrixRC<double, RowMajor>> vmat(v, n, n);
          vmat = svd.matrixV();
        }
      }
    }else if (a_layout == 'R' || a_layout == 'r'){
      Map<MatrixRC<double, RowMajor>> amat(a, m, n);
      auto svd = amat.bdcSvd(ComputeFullU | ComputeFullV);
      svec = svd.singularValues();
      if (u_layout == 'C' || u_layout == 'c'){
        Map<MatrixRC<double, ColMajor>> umat(u, m, m);
        umat = svd.matrixU();
        if (v_layout == 'C' || v_layout == 'c'){
          Map<MatrixRC<double, ColMajor>> vmat(v, n, n);
          vmat = svd.matrixV();
        }else if (v_layout == 'R' || v_layout == 'r'){
          Map<MatrixRC<double, RowMajor>> vmat(v, n, n);
          vmat = svd.matrixV();
        }
      }else if (u_layout == 'R' || u_layout == 'r'){
        Map<MatrixRC<double, RowMajor>> umat(u, m, m);
        umat = svd.matrixU();
        if (v_layout == 'C' || v_layout == 'c'){
          Map<MatrixRC<double, ColMajor>> vmat(v, n, n);
          vmat = svd.matrixV();
        }else if (v_layout == 'R' || v_layout == 'r'){
          Map<MatrixRC<double, RowMajor>> vmat(v, n, n);
          vmat = svd.matrixV();
        }
      }
    }
  }

  // Eigendecomposition: Compute eigenvalues and eigenvectors
  void CEIGEN_LITE_seigvals(const long n, float* a, char a_layout, std::complex<float>* eigvals){
    Map<VectorX<std::complex<float>>> evals(eigvals, n);
    if (a_layout == 'C' || a_layout == 'c'){
      Map<MatrixRC<float, ColMajor>> amat(a, n, n);
      EigenSolver<MatrixRC<float, ColMajor>> eig(amat, false);
      evals = eig.eigenvalues();

    }else if (a_layout == 'R' || a_layout == 'r'){
      Map<MatrixRC<float, RowMajor>> amat(a, n, n);
      EigenSolver<MatrixRC<float, RowMajor>> eig(amat, false);
      evals = eig.eigenvalues();
    }
  }
  void CEIGEN_LITE_deigvals(const long n, double* a, char a_layout, std::complex<double>* eigvals){
    Map<VectorX<std::complex<double>>> evals(eigvals, n);
    if (a_layout == 'C' || a_layout == 'c'){
      Map<MatrixRC<double, ColMajor>> amat(a, n, n);
      EigenSolver<MatrixRC<double, ColMajor>> eig(amat, false);
      evals = eig.eigenvalues();

    }else if (a_layout == 'R' || a_layout == 'r'){
      Map<MatrixRC<double, RowMajor>> amat(a, n, n);
      EigenSolver<MatrixRC<double, RowMajor>> eig(amat, false);
      evals = eig.eigenvalues();
    }
  }

  void CEIGEN_LITE_seigvecs(const long n, float* a, char a_layout,
                            std::complex<float>* eigvals,
                            std::complex<float>* eigvecs, char ev_layout){
    Map<VectorX<std::complex<float>>> evals(eigvals, n);
    if (a_layout == 'C' || a_layout == 'c'){
      Map<MatrixRC<float, ColMajor>> amat(a, n, n);
      EigenSolver<MatrixRC<float, ColMajor>> eig(amat, true);
      evals = eig.eigenvalues();
      if (ev_layout == 'C' || ev_layout == 'c'){
        Map<MatrixRC<std::complex<float>, ColMajor>> evecs(eigvecs, n, n);
        evecs = eig.eigenvectors();
      }else if (ev_layout == 'R' || ev_layout == 'r'){
        Map<MatrixRC<std::complex<float>, RowMajor>> evecs(eigvecs, n, n);
        evecs = eig.eigenvectors();
      }
    }else if (a_layout == 'R' || a_layout == 'r'){
      Map<MatrixRC<float, RowMajor>> amat(a, n, n);
      EigenSolver<MatrixRC<float, RowMajor>> eig(amat, true);
      evals = eig.eigenvalues();
      if (ev_layout == 'C' || ev_layout == 'c'){
        Map<MatrixRC<std::complex<float>, ColMajor>> evecs(eigvecs, n, n);
        evecs = eig.eigenvectors();
      }else if (ev_layout == 'R' || ev_layout == 'r'){
        Map<MatrixRC<std::complex<float>, RowMajor>> evecs(eigvecs, n, n);
        evecs = eig.eigenvectors();
      }
    }
  }
  void CEIGEN_LITE_deigvecs(const long n, double* a, char a_layout,
                            std::complex<double>* eigvals,
                            std::complex<double>* eigvecs, char ev_layout){
    Map<VectorX<std::complex<double>>> evals(eigvals, n);
    if (a_layout == 'C' || a_layout == 'c'){
      Map<MatrixRC<double, ColMajor>> amat(a, n, n);
      EigenSolver<MatrixRC<double, ColMajor>> eig(amat, true);
      evals = eig.eigenvalues();
      if (ev_layout == 'C' || ev_layout == 'c'){
        Map<MatrixRC<std::complex<double>, ColMajor>> evecs(eigvecs, n, n);
        evecs = eig.eigenvectors();
      }else if (ev_layout == 'R' || ev_layout == 'r'){
        Map<MatrixRC<std::complex<double>, RowMajor>> evecs(eigvecs, n, n);
        evecs = eig.eigenvectors();
      }
    }else if (a_layout == 'R' || a_layout == 'r'){
      Map<MatrixRC<double, RowMajor>> amat(a, n, n);
      EigenSolver<MatrixRC<double, RowMajor>> eig(amat, true);
      evals = eig.eigenvalues();
      if (ev_layout == 'C' || ev_layout == 'c'){
        Map<MatrixRC<std::complex<double>, ColMajor>> evecs(eigvecs, n, n);
        evecs = eig.eigenvectors();
      }else if (ev_layout == 'R' || ev_layout == 'r'){
        Map<MatrixRC<std::complex<double>, RowMajor>> evecs(eigvecs, n, n);
        evecs = eig.eigenvectors();
      }
    }
  }
};
