#ifndef _smirnov272_22_198_approx_HPP_
#define _smirnov272_22_198_approx_HPP_

// This is an automatically generated file from gen.py.
#include "common.hpp"

namespace smirnov272_22_198_approx {

template <typename Scalar>
void S_Add1(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& S5, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideS4 = S4.stride();
    const int strideS5 = S5.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    const Scalar *dataS4 = S4.data();
    const Scalar *dataS5 = S5.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataS1[i + j * strideS1] + Scalar(-(0.02 * (x))) * dataS2[i + j * strideS2] + Scalar(1.0 / (x * x)) * dataS3[i + j * strideS3] + Scalar(-(0.28 * (x))) * dataS4[i + j * strideS4] + Scalar(1.0 / (x * x)) * dataS5[i + j * strideS5];
        }
    }
}

template <typename Scalar>
void S_Add2(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideS4 = S4.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    const Scalar *dataS4 = S4.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(0.125) * dataS1[i + j * strideS1] + Scalar(-(0.02 * (x))) * dataS2[i + j * strideS2] + Scalar(0.125 * (1.0 / (x * x))) * dataS3[i + j * strideS3] + Scalar(1.0 / (x * x)) * dataS4[i + j * strideS4];
        }
    }
}

template <typename Scalar>
void S_Add3(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(1.0 / (x * x)) * dataS1[i + j * strideS1] + Scalar(-(0.343333333333 * (x))) * dataS2[i + j * strideS2] + Scalar(1.0 / (x * x)) * dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add4(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideS4 = S4.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    const Scalar *dataS4 = S4.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataS1[i + j * strideS1] + Scalar(0.02 * (x)) * dataS2[i + j * strideS2] + Scalar(1.0 / (x * x)) * dataS3[i + j * strideS3] + Scalar(-(1.0 / (x * x))) * dataS4[i + j * strideS4];
        }
    }
}

template <typename Scalar>
void S_Add5(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(0.02 * (x)) * dataS1[i + j * strideS1] + Scalar(1.0 / (x * x)) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add6(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideS4 = S4.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    const Scalar *dataS4 = S4.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(-(0.125 * (1.0 / (x * x)))) * dataS1[i + j * strideS1] + Scalar(1.0 / (x * x)) * dataS2[i + j * strideS2] + Scalar(0.07625 * (x)) * dataS3[i + j * strideS3] + Scalar(-(0.125 * (1.0 / (x * x)))) * dataS4[i + j * strideS4];
        }
    }
}

template <typename Scalar>
void S_Add7(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& S5, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideS4 = S4.stride();
    const int strideS5 = S5.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    const Scalar *dataS4 = S4.data();
    const Scalar *dataS5 = S5.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataS1[i + j * strideS1] + Scalar(-(1.0 / (x * x))) * dataS2[i + j * strideS2] + dataS3[i + j * strideS3] + Scalar(0.03 * (x)) * dataS4[i + j * strideS4] + Scalar(1.0 / (x * x)) * dataS5[i + j * strideS5];
        }
    }
}

template <typename Scalar>
void S_Add8(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(-(1.0 / (x * x))) * dataS1[i + j * strideS1] + Scalar(0.0633333333333 * (x)) * dataS2[i + j * strideS2] + Scalar(1.0 / (x * x)) * dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add9(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& S5, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideS4 = S4.stride();
    const int strideS5 = S5.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    const Scalar *dataS4 = S4.data();
    const Scalar *dataS5 = S5.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataS1[i + j * strideS1] + Scalar(1.0 / (x * x)) * dataS2[i + j * strideS2] + dataS3[i + j * strideS3] + Scalar(-(0.29 * (x))) * dataS4[i + j * strideS4] + Scalar(1.0 / (x * x)) * dataS5[i + j * strideS5];
        }
    }
}

template <typename Scalar>
void S_Add10(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(0.125) * dataS1[i + j * strideS1] -dataS2[i + j * strideS2] + Scalar(-0.0625) * dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add11(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(0.0333333333333 * (x)) * dataS1[i + j * strideS1] + Scalar(-(1.0 / (x * x))) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add12(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& S5, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideS4 = S4.stride();
    const int strideS5 = S5.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    const Scalar *dataS4 = S4.data();
    const Scalar *dataS5 = S5.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(1.0 / (x * x)) * dataS1[i + j * strideS1] + Scalar(-(0.28 * (x))) * dataS2[i + j * strideS2] + Scalar(1.0 / (x * x)) * dataS3[i + j * strideS3] + Scalar(-(0.02 * (x))) * dataS4[i + j * strideS4] + dataS5[i + j * strideS5];
        }
    }
}

template <typename Scalar>
void S_Add13(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideS4 = S4.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    const Scalar *dataS4 = S4.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(1.0 / (x * x)) * dataS1[i + j * strideS1] + Scalar(0.125 * (1.0 / (x * x))) * dataS2[i + j * strideS2] + Scalar(-(0.02 * (x))) * dataS3[i + j * strideS3] + Scalar(0.125) * dataS4[i + j * strideS4];
        }
    }
}

template <typename Scalar>
void S_Add14(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(1.0 / (x * x)) * dataS1[i + j * strideS1] + Scalar(-(0.343333333333 * (x))) * dataS2[i + j * strideS2] + Scalar(1.0 / (x * x)) * dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add15(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideS4 = S4.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    const Scalar *dataS4 = S4.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(-(1.0 / (x * x))) * dataS1[i + j * strideS1] + Scalar(1.0 / (x * x)) * dataS2[i + j * strideS2] + Scalar(0.02 * (x)) * dataS3[i + j * strideS3] + dataS4[i + j * strideS4];
        }
    }
}

template <typename Scalar>
void S_Add16(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(1.0 / (x * x)) * dataS1[i + j * strideS1] + Scalar(0.02 * (x)) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add17(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideS4 = S4.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    const Scalar *dataS4 = S4.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(-(0.125 * (1.0 / (x * x)))) * dataS1[i + j * strideS1] + Scalar(0.07625 * (x)) * dataS2[i + j * strideS2] + Scalar(1.0 / (x * x)) * dataS3[i + j * strideS3] + Scalar(-(0.125 * (1.0 / (x * x)))) * dataS4[i + j * strideS4];
        }
    }
}

template <typename Scalar>
void S_Add18(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& S5, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideS4 = S4.stride();
    const int strideS5 = S5.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    const Scalar *dataS4 = S4.data();
    const Scalar *dataS5 = S5.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(1.0 / (x * x)) * dataS1[i + j * strideS1] + Scalar(0.03 * (x)) * dataS2[i + j * strideS2] + dataS3[i + j * strideS3] + Scalar(-(1.0 / (x * x))) * dataS4[i + j * strideS4] -dataS5[i + j * strideS5];
        }
    }
}

template <typename Scalar>
void S_Add19(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(1.0 / (x * x)) * dataS1[i + j * strideS1] + Scalar(0.0633333333333 * (x)) * dataS2[i + j * strideS2] + Scalar(-(1.0 / (x * x))) * dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add20(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& S5, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideS4 = S4.stride();
    const int strideS5 = S5.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    const Scalar *dataS4 = S4.data();
    const Scalar *dataS5 = S5.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(1.0 / (x * x)) * dataS1[i + j * strideS1] + Scalar(-(0.29 * (x))) * dataS2[i + j * strideS2] + dataS3[i + j * strideS3] + Scalar(1.0 / (x * x)) * dataS4[i + j * strideS4] + dataS5[i + j * strideS5];
        }
    }
}

template <typename Scalar>
void S_Add21(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(-0.0625) * dataS1[i + j * strideS1] -dataS2[i + j * strideS2] + Scalar(0.125) * dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add22(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(-(1.0 / (x * x))) * dataS1[i + j * strideS1] + Scalar(0.0333333333333 * (x)) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void T_Add1(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideT3 = T3.stride();
    const int strideT4 = T4.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    const Scalar *dataT3 = T3.data();
    const Scalar *dataT4 = T4.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(x * x) * dataT1[i + j * strideT1] + Scalar(0.380952380952) * dataT2[i + j * strideT2] + Scalar(3.57142857143 * (x)) * dataT3[i + j * strideT3] + Scalar(-(0.419642857143 * (x * x * x * x))) * dataT4[i + j * strideT4];
        }
    }
}

template <typename Scalar>
void T_Add2(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideT3 = T3.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    const Scalar *dataT3 = T3.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(0.677248677249) * dataT1[i + j * strideT1] + Scalar(50.7936507937 * (x)) * dataT2[i + j * strideT2] + Scalar(1.14285714286 * (x * x * x * x)) * dataT3[i + j * strideT3];
        }
    }
}

template <typename Scalar>
void T_Add3(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& T5, Matrix<Scalar>& T6, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideT3 = T3.stride();
    const int strideT4 = T4.stride();
    const int strideT5 = T5.stride();
    const int strideT6 = T6.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    const Scalar *dataT3 = T3.data();
    const Scalar *dataT4 = T4.data();
    const Scalar *dataT5 = T5.data();
    const Scalar *dataT6 = T6.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(-(x * x)) * dataT1[i + j * strideT1] + Scalar(-0.333333333333) * dataT2[i + j * strideT2] + Scalar(1.25 * (1.0 / (x))) * dataT3[i + j * strideT3] + Scalar(1.5625 * (x * x * x * x)) * dataT4[i + j * strideT4] + Scalar(0.333333333333 * (x * x)) * dataT5[i + j * strideT5] + Scalar(0.125 * (x * x * x * x)) * dataT6[i + j * strideT6];
        }
    }
}

template <typename Scalar>
void T_Add4(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(0.296296296296) * dataT1[i + j * strideT1] + Scalar(-(2.77777777778 * (x))) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add5(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideT3 = T3.stride();
    const int strideT4 = T4.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    const Scalar *dataT3 = T3.data();
    const Scalar *dataT4 = T4.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataT1[i + j * strideT1] + Scalar(-(75.0 * (x))) * dataT2[i + j * strideT2] + Scalar(2.34375 * (1.0 / (x))) * dataT3[i + j * strideT3] + Scalar(x * x) * dataT4[i + j * strideT4];
        }
    }
}

template <typename Scalar>
void T_Add6(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(5.0 * (1.0 / (x))) * dataT1[i + j * strideT1] + Scalar(x * x * x * x) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add7(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(0.333333333333) * dataT1[i + j * strideT1] + Scalar(-(0.4375 * (x * x * x * x))) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add8(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideT3 = T3.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    const Scalar *dataT3 = T3.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(-0.333333333333) * dataT1[i + j * strideT1] + Scalar(1.25 * (1.0 / (x))) * dataT2[i + j * strideT2] + Scalar(0.333333333333 * (x * x)) * dataT3[i + j * strideT3];
        }
    }
}

template <typename Scalar>
void T_Add9(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(x * x) * dataT1[i + j * strideT1] + Scalar(0.333333333333) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add10(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(5.0 * (1.0 / (x * x * x))) * dataT1[i + j * strideT1];
        }
    }
}

template <typename Scalar>
void T_Add11(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(30.0 * (1.0 / (x))) * dataT1[i + j * strideT1] + Scalar(-(x * x)) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add12(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideT3 = T3.stride();
    const int strideT4 = T4.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    const Scalar *dataT3 = T3.data();
    const Scalar *dataT4 = T4.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(-(0.419642857143 * (x * x * x * x))) * dataT1[i + j * strideT1] + Scalar(3.57142857143 * (x)) * dataT2[i + j * strideT2] + Scalar(0.380952380952) * dataT3[i + j * strideT3] + Scalar(x * x) * dataT4[i + j * strideT4];
        }
    }
}

template <typename Scalar>
void T_Add13(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideT3 = T3.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    const Scalar *dataT3 = T3.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(1.14285714286 * (x * x * x * x)) * dataT1[i + j * strideT1] + Scalar(50.7936507937 * (x)) * dataT2[i + j * strideT2] + Scalar(0.677248677249) * dataT3[i + j * strideT3];
        }
    }
}

template <typename Scalar>
void T_Add14(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& T5, Matrix<Scalar>& T6, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideT3 = T3.stride();
    const int strideT4 = T4.stride();
    const int strideT5 = T5.stride();
    const int strideT6 = T6.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    const Scalar *dataT3 = T3.data();
    const Scalar *dataT4 = T4.data();
    const Scalar *dataT5 = T5.data();
    const Scalar *dataT6 = T6.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(0.125 * (x * x * x * x)) * dataT1[i + j * strideT1] + Scalar(0.333333333333 * (x * x)) * dataT2[i + j * strideT2] + Scalar(1.5625 * (x * x * x * x)) * dataT3[i + j * strideT3] + Scalar(1.25 * (1.0 / (x))) * dataT4[i + j * strideT4] + Scalar(-0.333333333333) * dataT5[i + j * strideT5] + Scalar(-(x * x)) * dataT6[i + j * strideT6];
        }
    }
}

template <typename Scalar>
void T_Add15(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(-(2.77777777778 * (x))) * dataT1[i + j * strideT1] + Scalar(0.296296296296) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add16(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideT3 = T3.stride();
    const int strideT4 = T4.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    const Scalar *dataT3 = T3.data();
    const Scalar *dataT4 = T4.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(x * x) * dataT1[i + j * strideT1] + Scalar(2.34375 * (1.0 / (x))) * dataT2[i + j * strideT2] + Scalar(-(75.0 * (x))) * dataT3[i + j * strideT3] -dataT4[i + j * strideT4];
        }
    }
}

template <typename Scalar>
void T_Add17(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(x * x * x * x) * dataT1[i + j * strideT1] + Scalar(5.0 * (1.0 / (x))) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add18(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(-(0.4375 * (x * x * x * x))) * dataT1[i + j * strideT1] + Scalar(0.333333333333) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add19(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideT3 = T3.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    const Scalar *dataT3 = T3.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(0.333333333333 * (x * x)) * dataT1[i + j * strideT1] + Scalar(1.25 * (1.0 / (x))) * dataT2[i + j * strideT2] + Scalar(-0.333333333333) * dataT3[i + j * strideT3];
        }
    }
}

template <typename Scalar>
void T_Add20(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(0.333333333333) * dataT1[i + j * strideT1] + Scalar(x * x) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add21(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(5.0 * (1.0 / (x * x * x))) * dataT1[i + j * strideT1];
        }
    }
}

template <typename Scalar>
void T_Add22(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(-(x * x)) * dataT1[i + j * strideT1] + Scalar(30.0 * (1.0 / (x))) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void M_Add1(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& M8, Matrix<Scalar>& M9, Matrix<Scalar>& M10, Matrix<Scalar>& M11, Matrix<Scalar>& M12, Matrix<Scalar>& M13, Matrix<Scalar>& M14, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
    const int strideM1 = M1.stride();
    const int strideM2 = M2.stride();
    const int strideM3 = M3.stride();
    const int strideM4 = M4.stride();
    const int strideM5 = M5.stride();
    const int strideM6 = M6.stride();
    const int strideM7 = M7.stride();
    const int strideM8 = M8.stride();
    const int strideM9 = M9.stride();
    const int strideM10 = M10.stride();
    const int strideM11 = M11.stride();
    const int strideM12 = M12.stride();
    const int strideM13 = M13.stride();
    const int strideM14 = M14.stride();
    const int strideC = C.stride();
    const Scalar *dataM1 = M1.data();
    const Scalar *dataM2 = M2.data();
    const Scalar *dataM3 = M3.data();
    const Scalar *dataM4 = M4.data();
    const Scalar *dataM5 = M5.data();
    const Scalar *dataM6 = M6.data();
    const Scalar *dataM7 = M7.data();
    const Scalar *dataM8 = M8.data();
    const Scalar *dataM9 = M9.data();
    const Scalar *dataM10 = M10.data();
    const Scalar *dataM11 = M11.data();
    const Scalar *dataM12 = M12.data();
    const Scalar *dataM13 = M13.data();
    const Scalar *dataM14 = M14.data();
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(1.0 / (x * x)) * dataM1[i + j * strideM1] + Scalar(-(1.0 / (x * x))) * dataM2[i + j * strideM2] + Scalar(1.0 / (x * x)) * dataM3[i + j * strideM3] + Scalar(-(1.0 / (x * x))) * dataM4[i + j * strideM4] + Scalar(-(0.666666666667 * (1.0 / (x * x)))) * dataM5[i + j * strideM5] + Scalar(1.0 / (x * x)) * dataM6[i + j * strideM6] + Scalar(1.0 / (x * x)) * dataM7[i + j * strideM7] + Scalar(1.0 / (x * x)) * dataM8[i + j * strideM8] + Scalar(1.3125) * dataM9[i + j * strideM9] + Scalar(3.0) * dataM10[i + j * strideM10] + dataM11[i + j * strideM11] + Scalar(2.5) * dataM12[i + j * strideM12] + Scalar(3.0) * dataM13[i + j * strideM13] + Scalar(2.5) * dataM14[i + j * strideM14] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(1.0 / (x * x)) * dataM1[i + j * strideM1] + Scalar(-(1.0 / (x * x))) * dataM2[i + j * strideM2] + Scalar(1.0 / (x * x)) * dataM3[i + j * strideM3] + Scalar(-(1.0 / (x * x))) * dataM4[i + j * strideM4] + Scalar(-(0.666666666667 * (1.0 / (x * x)))) * dataM5[i + j * strideM5] + Scalar(1.0 / (x * x)) * dataM6[i + j * strideM6] + Scalar(1.0 / (x * x)) * dataM7[i + j * strideM7] + Scalar(1.0 / (x * x)) * dataM8[i + j * strideM8] + Scalar(1.3125) * dataM9[i + j * strideM9] + Scalar(3.0) * dataM10[i + j * strideM10] + dataM11[i + j * strideM11] + Scalar(2.5) * dataM12[i + j * strideM12] + Scalar(3.0) * dataM13[i + j * strideM13] + Scalar(2.5) * dataM14[i + j * strideM14];
            }
        }
    }
}

template <typename Scalar>
void M_Add2(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& M8, Matrix<Scalar>& M9, Matrix<Scalar>& M10, Matrix<Scalar>& M11, Matrix<Scalar>& M12, Matrix<Scalar>& M13, Matrix<Scalar>& M14, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
    const int strideM1 = M1.stride();
    const int strideM2 = M2.stride();
    const int strideM3 = M3.stride();
    const int strideM4 = M4.stride();
    const int strideM5 = M5.stride();
    const int strideM6 = M6.stride();
    const int strideM7 = M7.stride();
    const int strideM8 = M8.stride();
    const int strideM9 = M9.stride();
    const int strideM10 = M10.stride();
    const int strideM11 = M11.stride();
    const int strideM12 = M12.stride();
    const int strideM13 = M13.stride();
    const int strideM14 = M14.stride();
    const int strideC = C.stride();
    const Scalar *dataM1 = M1.data();
    const Scalar *dataM2 = M2.data();
    const Scalar *dataM3 = M3.data();
    const Scalar *dataM4 = M4.data();
    const Scalar *dataM5 = M5.data();
    const Scalar *dataM6 = M6.data();
    const Scalar *dataM7 = M7.data();
    const Scalar *dataM8 = M8.data();
    const Scalar *dataM9 = M9.data();
    const Scalar *dataM10 = M10.data();
    const Scalar *dataM11 = M11.data();
    const Scalar *dataM12 = M12.data();
    const Scalar *dataM13 = M13.data();
    const Scalar *dataM14 = M14.data();
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(0.5) * dataM1[i + j * strideM1] + Scalar(0.8125) * dataM2[i + j * strideM2] + Scalar(0.5) * dataM3[i + j * strideM3] + Scalar(2.5) * dataM4[i + j * strideM4] + Scalar(0.666666666667) * dataM5[i + j * strideM5] + dataM6[i + j * strideM6] + Scalar(-2.5) * dataM7[i + j * strideM7] + Scalar(-5.0) * dataM8[i + j * strideM8] + dataM9[i + j * strideM9] + Scalar(-(1.0 / (x * x))) * dataM10[i + j * strideM10] + Scalar(0.125 * (1.0 / (x * x))) * dataM11[i + j * strideM11] + Scalar(-(1.0 / (x * x))) * dataM12[i + j * strideM12] + Scalar(-(1.0 / (x * x))) * dataM13[i + j * strideM13] + Scalar(1.0 / (x * x)) * dataM14[i + j * strideM14] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(0.5) * dataM1[i + j * strideM1] + Scalar(0.8125) * dataM2[i + j * strideM2] + Scalar(0.5) * dataM3[i + j * strideM3] + Scalar(2.5) * dataM4[i + j * strideM4] + Scalar(0.666666666667) * dataM5[i + j * strideM5] + dataM6[i + j * strideM6] + Scalar(-2.5) * dataM7[i + j * strideM7] + Scalar(-5.0) * dataM8[i + j * strideM8] + dataM9[i + j * strideM9] + Scalar(-(1.0 / (x * x))) * dataM10[i + j * strideM10] + Scalar(0.125 * (1.0 / (x * x))) * dataM11[i + j * strideM11] + Scalar(-(1.0 / (x * x))) * dataM12[i + j * strideM12] + Scalar(-(1.0 / (x * x))) * dataM13[i + j * strideM13] + Scalar(1.0 / (x * x)) * dataM14[i + j * strideM14];
            }
        }
    }
}

template <typename Scalar>
void M_Add3(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& M8, Matrix<Scalar>& M9, Matrix<Scalar>& M10, Matrix<Scalar>& M11, Matrix<Scalar>& M12, Matrix<Scalar>& M13, Matrix<Scalar>& M14, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
    const int strideM1 = M1.stride();
    const int strideM2 = M2.stride();
    const int strideM3 = M3.stride();
    const int strideM4 = M4.stride();
    const int strideM5 = M5.stride();
    const int strideM6 = M6.stride();
    const int strideM7 = M7.stride();
    const int strideM8 = M8.stride();
    const int strideM9 = M9.stride();
    const int strideM10 = M10.stride();
    const int strideM11 = M11.stride();
    const int strideM12 = M12.stride();
    const int strideM13 = M13.stride();
    const int strideM14 = M14.stride();
    const int strideC = C.stride();
    const Scalar *dataM1 = M1.data();
    const Scalar *dataM2 = M2.data();
    const Scalar *dataM3 = M3.data();
    const Scalar *dataM4 = M4.data();
    const Scalar *dataM5 = M5.data();
    const Scalar *dataM6 = M6.data();
    const Scalar *dataM7 = M7.data();
    const Scalar *dataM8 = M8.data();
    const Scalar *dataM9 = M9.data();
    const Scalar *dataM10 = M10.data();
    const Scalar *dataM11 = M11.data();
    const Scalar *dataM12 = M12.data();
    const Scalar *dataM13 = M13.data();
    const Scalar *dataM14 = M14.data();
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(-(1.0 / (x * x))) * dataM1[i + j * strideM1] + Scalar(0.125 * (1.0 / (x * x))) * dataM2[i + j * strideM2] + Scalar(-(1.0 / (x * x))) * dataM3[i + j * strideM3] + Scalar(-(1.0 / (x * x))) * dataM4[i + j * strideM4] + Scalar(1.0 / (x * x)) * dataM5[i + j * strideM5] + Scalar(0.5) * dataM6[i + j * strideM6] + Scalar(0.8125) * dataM7[i + j * strideM7] + Scalar(0.5) * dataM8[i + j * strideM8] + Scalar(2.5) * dataM9[i + j * strideM9] + Scalar(0.666666666667) * dataM10[i + j * strideM10] + dataM11[i + j * strideM11] + Scalar(-2.5) * dataM12[i + j * strideM12] + Scalar(-5.0) * dataM13[i + j * strideM13] + dataM14[i + j * strideM14] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(-(1.0 / (x * x))) * dataM1[i + j * strideM1] + Scalar(0.125 * (1.0 / (x * x))) * dataM2[i + j * strideM2] + Scalar(-(1.0 / (x * x))) * dataM3[i + j * strideM3] + Scalar(-(1.0 / (x * x))) * dataM4[i + j * strideM4] + Scalar(1.0 / (x * x)) * dataM5[i + j * strideM5] + Scalar(0.5) * dataM6[i + j * strideM6] + Scalar(0.8125) * dataM7[i + j * strideM7] + Scalar(0.5) * dataM8[i + j * strideM8] + Scalar(2.5) * dataM9[i + j * strideM9] + Scalar(0.666666666667) * dataM10[i + j * strideM10] + dataM11[i + j * strideM11] + Scalar(-2.5) * dataM12[i + j * strideM12] + Scalar(-5.0) * dataM13[i + j * strideM13] + dataM14[i + j * strideM14];
            }
        }
    }
}

template <typename Scalar>
void M_Add4(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& M8, Matrix<Scalar>& M9, Matrix<Scalar>& M10, Matrix<Scalar>& M11, Matrix<Scalar>& M12, Matrix<Scalar>& M13, Matrix<Scalar>& M14, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
    const int strideM1 = M1.stride();
    const int strideM2 = M2.stride();
    const int strideM3 = M3.stride();
    const int strideM4 = M4.stride();
    const int strideM5 = M5.stride();
    const int strideM6 = M6.stride();
    const int strideM7 = M7.stride();
    const int strideM8 = M8.stride();
    const int strideM9 = M9.stride();
    const int strideM10 = M10.stride();
    const int strideM11 = M11.stride();
    const int strideM12 = M12.stride();
    const int strideM13 = M13.stride();
    const int strideM14 = M14.stride();
    const int strideC = C.stride();
    const Scalar *dataM1 = M1.data();
    const Scalar *dataM2 = M2.data();
    const Scalar *dataM3 = M3.data();
    const Scalar *dataM4 = M4.data();
    const Scalar *dataM5 = M5.data();
    const Scalar *dataM6 = M6.data();
    const Scalar *dataM7 = M7.data();
    const Scalar *dataM8 = M8.data();
    const Scalar *dataM9 = M9.data();
    const Scalar *dataM10 = M10.data();
    const Scalar *dataM11 = M11.data();
    const Scalar *dataM12 = M12.data();
    const Scalar *dataM13 = M13.data();
    const Scalar *dataM14 = M14.data();
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(1.3125) * dataM1[i + j * strideM1] + Scalar(3.0) * dataM2[i + j * strideM2] + dataM3[i + j * strideM3] + Scalar(2.5) * dataM4[i + j * strideM4] + Scalar(3.0) * dataM5[i + j * strideM5] + Scalar(2.5) * dataM6[i + j * strideM6] + Scalar(1.0 / (x * x)) * dataM7[i + j * strideM7] + Scalar(-(1.0 / (x * x))) * dataM8[i + j * strideM8] + Scalar(1.0 / (x * x)) * dataM9[i + j * strideM9] + Scalar(-(1.0 / (x * x))) * dataM10[i + j * strideM10] + Scalar(-(0.666666666667 * (1.0 / (x * x)))) * dataM11[i + j * strideM11] + Scalar(1.0 / (x * x)) * dataM12[i + j * strideM12] + Scalar(1.0 / (x * x)) * dataM13[i + j * strideM13] + Scalar(1.0 / (x * x)) * dataM14[i + j * strideM14] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(1.3125) * dataM1[i + j * strideM1] + Scalar(3.0) * dataM2[i + j * strideM2] + dataM3[i + j * strideM3] + Scalar(2.5) * dataM4[i + j * strideM4] + Scalar(3.0) * dataM5[i + j * strideM5] + Scalar(2.5) * dataM6[i + j * strideM6] + Scalar(1.0 / (x * x)) * dataM7[i + j * strideM7] + Scalar(-(1.0 / (x * x))) * dataM8[i + j * strideM8] + Scalar(1.0 / (x * x)) * dataM9[i + j * strideM9] + Scalar(-(1.0 / (x * x))) * dataM10[i + j * strideM10] + Scalar(-(0.666666666667 * (1.0 / (x * x)))) * dataM11[i + j * strideM11] + Scalar(1.0 / (x * x)) * dataM12[i + j * strideM12] + Scalar(1.0 / (x * x)) * dataM13[i + j * strideM13] + Scalar(1.0 / (x * x)) * dataM14[i + j * strideM14];
            }
        }
    }
}

template <typename Scalar>
void FastMatmulRecursive(LockAndCounter& locker, MemoryManager<Scalar>& mem_mngr, Matrix<Scalar>& A, Matrix<Scalar>& B, Matrix<Scalar>& C, int total_steps, int steps_left, int start_index, double x, int num_threads, Scalar beta) {
    // Update multipliers
    C.UpdateMultiplier(A.multiplier());
    C.UpdateMultiplier(B.multiplier());
    A.set_multiplier(Scalar(1.0));
    B.set_multiplier(Scalar(1.0));
    // Base case for recursion
    if (steps_left == 0) {
        MatMul(A, B, C);
        return;
    }

    Matrix<Scalar> A11 = A.Subblock(2, 7, 1, 1);
    Matrix<Scalar> A12 = A.Subblock(2, 7, 1, 2);
    Matrix<Scalar> A13 = A.Subblock(2, 7, 1, 3);
    Matrix<Scalar> A14 = A.Subblock(2, 7, 1, 4);
    Matrix<Scalar> A15 = A.Subblock(2, 7, 1, 5);
    Matrix<Scalar> A16 = A.Subblock(2, 7, 1, 6);
    Matrix<Scalar> A17 = A.Subblock(2, 7, 1, 7);
    Matrix<Scalar> A21 = A.Subblock(2, 7, 2, 1);
    Matrix<Scalar> A22 = A.Subblock(2, 7, 2, 2);
    Matrix<Scalar> A23 = A.Subblock(2, 7, 2, 3);
    Matrix<Scalar> A24 = A.Subblock(2, 7, 2, 4);
    Matrix<Scalar> A25 = A.Subblock(2, 7, 2, 5);
    Matrix<Scalar> A26 = A.Subblock(2, 7, 2, 6);
    Matrix<Scalar> A27 = A.Subblock(2, 7, 2, 7);
    Matrix<Scalar> B11 = B.Subblock(7, 2, 1, 1);
    Matrix<Scalar> B12 = B.Subblock(7, 2, 1, 2);
    Matrix<Scalar> B21 = B.Subblock(7, 2, 2, 1);
    Matrix<Scalar> B22 = B.Subblock(7, 2, 2, 2);
    Matrix<Scalar> B31 = B.Subblock(7, 2, 3, 1);
    Matrix<Scalar> B32 = B.Subblock(7, 2, 3, 2);
    Matrix<Scalar> B41 = B.Subblock(7, 2, 4, 1);
    Matrix<Scalar> B42 = B.Subblock(7, 2, 4, 2);
    Matrix<Scalar> B51 = B.Subblock(7, 2, 5, 1);
    Matrix<Scalar> B52 = B.Subblock(7, 2, 5, 2);
    Matrix<Scalar> B61 = B.Subblock(7, 2, 6, 1);
    Matrix<Scalar> B62 = B.Subblock(7, 2, 6, 2);
    Matrix<Scalar> B71 = B.Subblock(7, 2, 7, 1);
    Matrix<Scalar> B72 = B.Subblock(7, 2, 7, 2);
    Matrix<Scalar> C11 = C.Subblock(2, 2, 1, 1);
    Matrix<Scalar> C12 = C.Subblock(2, 2, 1, 2);
    Matrix<Scalar> C21 = C.Subblock(2, 2, 2, 1);
    Matrix<Scalar> C22 = C.Subblock(2, 2, 2, 2);


    // Matrices to store the results of multiplications.
#ifdef _PARALLEL_
    Matrix<Scalar> M1(mem_mngr.GetMem(start_index, 1, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M2(mem_mngr.GetMem(start_index, 2, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M3(mem_mngr.GetMem(start_index, 3, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M4(mem_mngr.GetMem(start_index, 4, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M5(mem_mngr.GetMem(start_index, 5, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M6(mem_mngr.GetMem(start_index, 6, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M7(mem_mngr.GetMem(start_index, 7, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M8(mem_mngr.GetMem(start_index, 8, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M9(mem_mngr.GetMem(start_index, 9, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M10(mem_mngr.GetMem(start_index, 10, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M11(mem_mngr.GetMem(start_index, 11, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M12(mem_mngr.GetMem(start_index, 12, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M13(mem_mngr.GetMem(start_index, 13, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M14(mem_mngr.GetMem(start_index, 14, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M15(mem_mngr.GetMem(start_index, 15, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M16(mem_mngr.GetMem(start_index, 16, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M17(mem_mngr.GetMem(start_index, 17, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M18(mem_mngr.GetMem(start_index, 18, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M19(mem_mngr.GetMem(start_index, 19, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M20(mem_mngr.GetMem(start_index, 20, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M21(mem_mngr.GetMem(start_index, 21, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M22(mem_mngr.GetMem(start_index, 22, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
#else
    Matrix<Scalar> M1(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M2(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M3(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M4(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M5(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M6(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M7(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M8(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M9(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M10(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M11(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M12(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M13(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M14(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M15(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M16(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M17(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M18(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M19(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M20(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M21(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M22(C11.m(), C11.n(), C.multiplier());
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
    bool sequential1 = should_launch_task(22, total_steps, steps_left, start_index, 1, num_threads);
    bool sequential2 = should_launch_task(22, total_steps, steps_left, start_index, 2, num_threads);
    bool sequential3 = should_launch_task(22, total_steps, steps_left, start_index, 3, num_threads);
    bool sequential4 = should_launch_task(22, total_steps, steps_left, start_index, 4, num_threads);
    bool sequential5 = should_launch_task(22, total_steps, steps_left, start_index, 5, num_threads);
    bool sequential6 = should_launch_task(22, total_steps, steps_left, start_index, 6, num_threads);
    bool sequential7 = should_launch_task(22, total_steps, steps_left, start_index, 7, num_threads);
    bool sequential8 = should_launch_task(22, total_steps, steps_left, start_index, 8, num_threads);
    bool sequential9 = should_launch_task(22, total_steps, steps_left, start_index, 9, num_threads);
    bool sequential10 = should_launch_task(22, total_steps, steps_left, start_index, 10, num_threads);
    bool sequential11 = should_launch_task(22, total_steps, steps_left, start_index, 11, num_threads);
    bool sequential12 = should_launch_task(22, total_steps, steps_left, start_index, 12, num_threads);
    bool sequential13 = should_launch_task(22, total_steps, steps_left, start_index, 13, num_threads);
    bool sequential14 = should_launch_task(22, total_steps, steps_left, start_index, 14, num_threads);
    bool sequential15 = should_launch_task(22, total_steps, steps_left, start_index, 15, num_threads);
    bool sequential16 = should_launch_task(22, total_steps, steps_left, start_index, 16, num_threads);
    bool sequential17 = should_launch_task(22, total_steps, steps_left, start_index, 17, num_threads);
    bool sequential18 = should_launch_task(22, total_steps, steps_left, start_index, 18, num_threads);
    bool sequential19 = should_launch_task(22, total_steps, steps_left, start_index, 19, num_threads);
    bool sequential20 = should_launch_task(22, total_steps, steps_left, start_index, 20, num_threads);
    bool sequential21 = should_launch_task(22, total_steps, steps_left, start_index, 21, num_threads);
    bool sequential22 = should_launch_task(22, total_steps, steps_left, start_index, 22, num_threads);
#else
    bool sequential1 = false;
    bool sequential2 = false;
    bool sequential3 = false;
    bool sequential4 = false;
    bool sequential5 = false;
    bool sequential6 = false;
    bool sequential7 = false;
    bool sequential8 = false;
    bool sequential9 = false;
    bool sequential10 = false;
    bool sequential11 = false;
    bool sequential12 = false;
    bool sequential13 = false;
    bool sequential14 = false;
    bool sequential15 = false;
    bool sequential16 = false;
    bool sequential17 = false;
    bool sequential18 = false;
    bool sequential19 = false;
    bool sequential20 = false;
    bool sequential21 = false;
    bool sequential22 = false;
#endif



    // M1 = (1.0 * A11 + -(0.02 * (x)) * A12 + 1.0 / (x * x) * A13 + -(0.28 * (x)) * A22 + 1.0 / (x * x) * A23) * (x * x * B11 + 0.380952380952 * B12 + 3.57142857143 * (x) * B21 + -(0.419642857143 * (x * x * x * x)) * B31)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential1) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S1(mem_mngr.GetMem(start_index, 1, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S1(A11.m(), A11.n());
#endif
    S_Add1(A11, A12, A13, A22, A23, S1, x, sequential1);
#ifdef _PARALLEL_
    Matrix<Scalar> T1(mem_mngr.GetMem(start_index, 1, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T1(B11.m(), B11.n());
#endif
    T_Add1(B11, B12, B21, B31, T1, x, sequential1);
    FastMatmulRecursive(locker, mem_mngr, S1, T1, M1, total_steps, steps_left - 1, (start_index + 1 - 1) * 22, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S1.deallocate();
    T1.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(22, total_steps, steps_left, start_index, 1, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M2 = (0.125 * A11 + -(0.02 * (x)) * A12 + 0.125 * (1.0 / (x * x)) * A13 + 1.0 / (x * x) * A23) * (0.677248677249 * B12 + 50.7936507937 * (x) * B21 + 1.14285714286 * (x * x * x * x) * B31)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential2) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S2(mem_mngr.GetMem(start_index, 2, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S2(A11.m(), A11.n());
#endif
    S_Add2(A11, A12, A13, A23, S2, x, sequential2);
#ifdef _PARALLEL_
    Matrix<Scalar> T2(mem_mngr.GetMem(start_index, 2, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T2(B11.m(), B11.n());
#endif
    T_Add2(B12, B21, B31, T2, x, sequential2);
    FastMatmulRecursive(locker, mem_mngr, S2, T2, M2, total_steps, steps_left - 1, (start_index + 2 - 1) * 22, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S2.deallocate();
    T2.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(22, total_steps, steps_left, start_index, 2, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M3 = (1.0 / (x * x) * A13 + -(0.343333333333 * (x)) * A22 + 1.0 / (x * x) * A23) * (-(x * x) * B11 + -0.333333333333 * B12 + 1.25 * (1.0 / (x)) * B22 + 1.5625 * (x * x * x * x) * B31 + 0.333333333333 * (x * x) * B32 + 0.125 * (x * x * x * x) * B41)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential3) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S3(mem_mngr.GetMem(start_index, 3, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S3(A11.m(), A11.n());
#endif
    S_Add3(A13, A22, A23, S3, x, sequential3);
#ifdef _PARALLEL_
    Matrix<Scalar> T3(mem_mngr.GetMem(start_index, 3, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T3(B11.m(), B11.n());
#endif
    T_Add3(B11, B12, B22, B31, B32, B41, T3, x, sequential3);
    FastMatmulRecursive(locker, mem_mngr, S3, T3, M3, total_steps, steps_left - 1, (start_index + 3 - 1) * 22, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S3.deallocate();
    T3.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(22, total_steps, steps_left, start_index, 3, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M4 = (1.0 * A11 + 0.02 * (x) * A12 + 1.0 / (x * x) * A13 + -(1.0 / (x * x)) * A23) * (0.296296296296 * B12 + -(2.77777777778 * (x)) * B21)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential4) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S4(mem_mngr.GetMem(start_index, 4, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S4(A11.m(), A11.n());
#endif
    S_Add4(A11, A12, A13, A23, S4, x, sequential4);
#ifdef _PARALLEL_
    Matrix<Scalar> T4(mem_mngr.GetMem(start_index, 4, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T4(B11.m(), B11.n());
#endif
    T_Add4(B12, B21, T4, x, sequential4);
    FastMatmulRecursive(locker, mem_mngr, S4, T4, M4, total_steps, steps_left - 1, (start_index + 4 - 1) * 22, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S4.deallocate();
    T4.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(22, total_steps, steps_left, start_index, 4, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M5 = (0.02 * (x) * A22 + 1.0 / (x * x) * A23) * (-1.0 * B12 + -(75.0 * (x)) * B21 + 2.34375 * (1.0 / (x)) * B22 + x * x * B32)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential5) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S5(mem_mngr.GetMem(start_index, 5, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S5(A11.m(), A11.n());
#endif
    S_Add5(A22, A23, S5, x, sequential5);
#ifdef _PARALLEL_
    Matrix<Scalar> T5(mem_mngr.GetMem(start_index, 5, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T5(B11.m(), B11.n());
#endif
    T_Add5(B12, B21, B22, B32, T5, x, sequential5);
    FastMatmulRecursive(locker, mem_mngr, S5, T5, M5, total_steps, steps_left - 1, (start_index + 5 - 1) * 22, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S5.deallocate();
    T5.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(22, total_steps, steps_left, start_index, 5, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M6 = (-(0.125 * (1.0 / (x * x))) * A13 + 1.0 / (x * x) * A14 + 0.07625 * (x) * A22 + -(0.125 * (1.0 / (x * x))) * A23) * (5.0 * (1.0 / (x)) * B22 + x * x * x * x * B41)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential6) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S6(mem_mngr.GetMem(start_index, 6, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S6(A11.m(), A11.n());
#endif
    S_Add6(A13, A14, A22, A23, S6, x, sequential6);
#ifdef _PARALLEL_
    Matrix<Scalar> T6(mem_mngr.GetMem(start_index, 6, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T6(B11.m(), B11.n());
#endif
    T_Add6(B22, B41, T6, x, sequential6);
    FastMatmulRecursive(locker, mem_mngr, S6, T6, M6, total_steps, steps_left - 1, (start_index + 6 - 1) * 22, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S6.deallocate();
    T6.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(22, total_steps, steps_left, start_index, 6, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M7 = (-1.0 * A11 + -(1.0 / (x * x)) * A13 + 1.0 * A21 + 0.03 * (x) * A22 + 1.0 / (x * x) * A23) * (0.333333333333 * B12 + -(0.4375 * (x * x * x * x)) * B31)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential7) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S7(mem_mngr.GetMem(start_index, 7, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S7(A11.m(), A11.n());
#endif
    S_Add7(A11, A13, A21, A22, A23, S7, x, sequential7);
#ifdef _PARALLEL_
    Matrix<Scalar> T7(mem_mngr.GetMem(start_index, 7, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T7(B11.m(), B11.n());
#endif
    T_Add7(B12, B31, T7, x, sequential7);
    FastMatmulRecursive(locker, mem_mngr, S7, T7, M7, total_steps, steps_left - 1, (start_index + 7 - 1) * 22, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S7.deallocate();
    T7.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(22, total_steps, steps_left, start_index, 7, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M8 = (-(1.0 / (x * x)) * A13 + 0.0633333333333 * (x) * A22 + 1.0 / (x * x) * A23) * (-0.333333333333 * B12 + 1.25 * (1.0 / (x)) * B22 + 0.333333333333 * (x * x) * B32)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential8) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S8(mem_mngr.GetMem(start_index, 8, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S8(A11.m(), A11.n());
#endif
    S_Add8(A13, A22, A23, S8, x, sequential8);
#ifdef _PARALLEL_
    Matrix<Scalar> T8(mem_mngr.GetMem(start_index, 8, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T8(B11.m(), B11.n());
#endif
    T_Add8(B12, B22, B32, T8, x, sequential8);
    FastMatmulRecursive(locker, mem_mngr, S8, T8, M8, total_steps, steps_left - 1, (start_index + 8 - 1) * 22, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S8.deallocate();
    T8.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(22, total_steps, steps_left, start_index, 8, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M9 = (1.0 * A11 + 1.0 / (x * x) * A13 + 1.0 * A21 + -(0.29 * (x)) * A22 + 1.0 / (x * x) * A23) * (x * x * B11 + 0.333333333333 * B12)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential9) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S9(mem_mngr.GetMem(start_index, 9, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S9(A11.m(), A11.n());
#endif
    S_Add9(A11, A13, A21, A22, A23, S9, x, sequential9);
#ifdef _PARALLEL_
    Matrix<Scalar> T9(mem_mngr.GetMem(start_index, 9, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T9(B11.m(), B11.n());
#endif
    T_Add9(B11, B12, T9, x, sequential9);
    FastMatmulRecursive(locker, mem_mngr, S9, T9, M9, total_steps, steps_left - 1, (start_index + 9 - 1) * 22, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S9.deallocate();
    T9.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(22, total_steps, steps_left, start_index, 9, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M10 = (0.125 * A13 + -1.0 * A14 + -0.0625 * A23) * (5.0 * (1.0 / (x * x * x)) * B22)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential10) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S10(mem_mngr.GetMem(start_index, 10, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S10(A11.m(), A11.n());
#endif
    S_Add10(A13, A14, A23, S10, x, sequential10);
    M10.UpdateMultiplier(Scalar(5.0 * (1.0 / (x * x * x))));
    FastMatmulRecursive(locker, mem_mngr, S10, B22, M10, total_steps, steps_left - 1, (start_index + 10 - 1) * 22, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S10.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(22, total_steps, steps_left, start_index, 10, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M11 = (0.0333333333333 * (x) * A12 + -(1.0 / (x * x)) * A14) * (30.0 * (1.0 / (x)) * B22 + -(x * x) * B42)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential11) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S11(mem_mngr.GetMem(start_index, 11, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S11(A11.m(), A11.n());
#endif
    S_Add11(A12, A14, S11, x, sequential11);
#ifdef _PARALLEL_
    Matrix<Scalar> T11(mem_mngr.GetMem(start_index, 11, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T11(B11.m(), B11.n());
#endif
    T_Add11(B22, B42, T11, x, sequential11);
    FastMatmulRecursive(locker, mem_mngr, S11, T11, M11, total_steps, steps_left - 1, (start_index + 11 - 1) * 22, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S11.deallocate();
    T11.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(22, total_steps, steps_left, start_index, 11, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M12 = (1.0 / (x * x) * A15 + -(0.28 * (x)) * A16 + 1.0 / (x * x) * A25 + -(0.02 * (x)) * A26 + 1.0 * A27) * (-(0.419642857143 * (x * x * x * x)) * B52 + 3.57142857143 * (x) * B62 + 0.380952380952 * B71 + x * x * B72)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential12) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S12(mem_mngr.GetMem(start_index, 12, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S12(A11.m(), A11.n());
#endif
    S_Add12(A15, A16, A25, A26, A27, S12, x, sequential12);
#ifdef _PARALLEL_
    Matrix<Scalar> T12(mem_mngr.GetMem(start_index, 12, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T12(B11.m(), B11.n());
#endif
    T_Add12(B52, B62, B71, B72, T12, x, sequential12);
    FastMatmulRecursive(locker, mem_mngr, S12, T12, M12, total_steps, steps_left - 1, (start_index + 12 - 1) * 22, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S12.deallocate();
    T12.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(22, total_steps, steps_left, start_index, 12, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M13 = (1.0 / (x * x) * A15 + 0.125 * (1.0 / (x * x)) * A25 + -(0.02 * (x)) * A26 + 0.125 * A27) * (1.14285714286 * (x * x * x * x) * B52 + 50.7936507937 * (x) * B62 + 0.677248677249 * B71)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential13) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S13(mem_mngr.GetMem(start_index, 13, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S13(A11.m(), A11.n());
#endif
    S_Add13(A15, A25, A26, A27, S13, x, sequential13);
#ifdef _PARALLEL_
    Matrix<Scalar> T13(mem_mngr.GetMem(start_index, 13, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T13(B11.m(), B11.n());
#endif
    T_Add13(B52, B62, B71, T13, x, sequential13);
    FastMatmulRecursive(locker, mem_mngr, S13, T13, M13, total_steps, steps_left - 1, (start_index + 13 - 1) * 22, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S13.deallocate();
    T13.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(22, total_steps, steps_left, start_index, 13, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M14 = (1.0 / (x * x) * A15 + -(0.343333333333 * (x)) * A16 + 1.0 / (x * x) * A25) * (0.125 * (x * x * x * x) * B42 + 0.333333333333 * (x * x) * B51 + 1.5625 * (x * x * x * x) * B52 + 1.25 * (1.0 / (x)) * B61 + -0.333333333333 * B71 + -(x * x) * B72)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential14) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S14(mem_mngr.GetMem(start_index, 14, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S14(A11.m(), A11.n());
#endif
    S_Add14(A15, A16, A25, S14, x, sequential14);
#ifdef _PARALLEL_
    Matrix<Scalar> T14(mem_mngr.GetMem(start_index, 14, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T14(B11.m(), B11.n());
#endif
    T_Add14(B42, B51, B52, B61, B71, B72, T14, x, sequential14);
    FastMatmulRecursive(locker, mem_mngr, S14, T14, M14, total_steps, steps_left - 1, (start_index + 14 - 1) * 22, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S14.deallocate();
    T14.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(22, total_steps, steps_left, start_index, 14, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M15 = (-(1.0 / (x * x)) * A15 + 1.0 / (x * x) * A25 + 0.02 * (x) * A26 + 1.0 * A27) * (-(2.77777777778 * (x)) * B62 + 0.296296296296 * B71)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential15) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S15(mem_mngr.GetMem(start_index, 15, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S15(A11.m(), A11.n());
#endif
    S_Add15(A15, A25, A26, A27, S15, x, sequential15);
#ifdef _PARALLEL_
    Matrix<Scalar> T15(mem_mngr.GetMem(start_index, 15, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T15(B11.m(), B11.n());
#endif
    T_Add15(B62, B71, T15, x, sequential15);
    FastMatmulRecursive(locker, mem_mngr, S15, T15, M15, total_steps, steps_left - 1, (start_index + 15 - 1) * 22, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S15.deallocate();
    T15.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(22, total_steps, steps_left, start_index, 15, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M16 = (1.0 / (x * x) * A15 + 0.02 * (x) * A16) * (x * x * B51 + 2.34375 * (1.0 / (x)) * B61 + -(75.0 * (x)) * B62 + -1.0 * B71)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential16) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S16(mem_mngr.GetMem(start_index, 16, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S16(A11.m(), A11.n());
#endif
    S_Add16(A15, A16, S16, x, sequential16);
#ifdef _PARALLEL_
    Matrix<Scalar> T16(mem_mngr.GetMem(start_index, 16, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T16(B11.m(), B11.n());
#endif
    T_Add16(B51, B61, B62, B71, T16, x, sequential16);
    FastMatmulRecursive(locker, mem_mngr, S16, T16, M16, total_steps, steps_left - 1, (start_index + 16 - 1) * 22, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S16.deallocate();
    T16.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(22, total_steps, steps_left, start_index, 16, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M17 = (-(0.125 * (1.0 / (x * x))) * A15 + 0.07625 * (x) * A16 + 1.0 / (x * x) * A24 + -(0.125 * (1.0 / (x * x))) * A25) * (x * x * x * x * B42 + 5.0 * (1.0 / (x)) * B61)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential17) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S17(mem_mngr.GetMem(start_index, 17, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S17(A11.m(), A11.n());
#endif
    S_Add17(A15, A16, A24, A25, S17, x, sequential17);
#ifdef _PARALLEL_
    Matrix<Scalar> T17(mem_mngr.GetMem(start_index, 17, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T17(B11.m(), B11.n());
#endif
    T_Add17(B42, B61, T17, x, sequential17);
    FastMatmulRecursive(locker, mem_mngr, S17, T17, M17, total_steps, steps_left - 1, (start_index + 17 - 1) * 22, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S17.deallocate();
    T17.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(22, total_steps, steps_left, start_index, 17, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M18 = (1.0 / (x * x) * A15 + 0.03 * (x) * A16 + 1.0 * A17 + -(1.0 / (x * x)) * A25 + -1.0 * A27) * (-(0.4375 * (x * x * x * x)) * B52 + 0.333333333333 * B71)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential18) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S18(mem_mngr.GetMem(start_index, 18, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S18(A11.m(), A11.n());
#endif
    S_Add18(A15, A16, A17, A25, A27, S18, x, sequential18);
#ifdef _PARALLEL_
    Matrix<Scalar> T18(mem_mngr.GetMem(start_index, 18, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T18(B11.m(), B11.n());
#endif
    T_Add18(B52, B71, T18, x, sequential18);
    FastMatmulRecursive(locker, mem_mngr, S18, T18, M18, total_steps, steps_left - 1, (start_index + 18 - 1) * 22, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S18.deallocate();
    T18.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(22, total_steps, steps_left, start_index, 18, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M19 = (1.0 / (x * x) * A15 + 0.0633333333333 * (x) * A16 + -(1.0 / (x * x)) * A25) * (0.333333333333 * (x * x) * B51 + 1.25 * (1.0 / (x)) * B61 + -0.333333333333 * B71)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential19) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S19(mem_mngr.GetMem(start_index, 19, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S19(A11.m(), A11.n());
#endif
    S_Add19(A15, A16, A25, S19, x, sequential19);
#ifdef _PARALLEL_
    Matrix<Scalar> T19(mem_mngr.GetMem(start_index, 19, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T19(B11.m(), B11.n());
#endif
    T_Add19(B51, B61, B71, T19, x, sequential19);
    FastMatmulRecursive(locker, mem_mngr, S19, T19, M19, total_steps, steps_left - 1, (start_index + 19 - 1) * 22, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S19.deallocate();
    T19.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(22, total_steps, steps_left, start_index, 19, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M20 = (1.0 / (x * x) * A15 + -(0.29 * (x)) * A16 + 1.0 * A17 + 1.0 / (x * x) * A25 + 1.0 * A27) * (0.333333333333 * B71 + x * x * B72)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential20) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S20(mem_mngr.GetMem(start_index, 20, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S20(A11.m(), A11.n());
#endif
    S_Add20(A15, A16, A17, A25, A27, S20, x, sequential20);
#ifdef _PARALLEL_
    Matrix<Scalar> T20(mem_mngr.GetMem(start_index, 20, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T20(B11.m(), B11.n());
#endif
    T_Add20(B71, B72, T20, x, sequential20);
    FastMatmulRecursive(locker, mem_mngr, S20, T20, M20, total_steps, steps_left - 1, (start_index + 20 - 1) * 22, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S20.deallocate();
    T20.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(22, total_steps, steps_left, start_index, 20, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M21 = (-0.0625 * A15 + -1.0 * A24 + 0.125 * A25) * (5.0 * (1.0 / (x * x * x)) * B61)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential21) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S21(mem_mngr.GetMem(start_index, 21, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S21(A11.m(), A11.n());
#endif
    S_Add21(A15, A24, A25, S21, x, sequential21);
    M21.UpdateMultiplier(Scalar(5.0 * (1.0 / (x * x * x))));
    FastMatmulRecursive(locker, mem_mngr, S21, B61, M21, total_steps, steps_left - 1, (start_index + 21 - 1) * 22, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S21.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(22, total_steps, steps_left, start_index, 21, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M22 = (-(1.0 / (x * x)) * A24 + 0.0333333333333 * (x) * A26) * (-(x * x) * B41 + 30.0 * (1.0 / (x)) * B61)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential22) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S22(mem_mngr.GetMem(start_index, 22, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S22(A11.m(), A11.n());
#endif
    S_Add22(A24, A26, S22, x, sequential22);
#ifdef _PARALLEL_
    Matrix<Scalar> T22(mem_mngr.GetMem(start_index, 22, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T22(B11.m(), B11.n());
#endif
    T_Add22(B41, B61, T22, x, sequential22);
    FastMatmulRecursive(locker, mem_mngr, S22, T22, M22, total_steps, steps_left - 1, (start_index + 22 - 1) * 22, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S22.deallocate();
    T22.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(22, total_steps, steps_left, start_index, 22, num_threads)) {
# pragma omp taskwait
    }
#endif

    M_Add1(M1, M2, M3, M4, M5, M6, M8, M10, M13, M15, M16, M17, M18, M21, C11, x, false, beta);
    M_Add2(M1, M2, M3, M4, M5, M6, M8, M10, M11, M12, M13, M15, M18, M20, C12, x, false, beta);
    M_Add3(M1, M2, M4, M7, M9, M12, M13, M14, M15, M16, M17, M19, M21, M22, C21, x, false, beta);
    M_Add4(M2, M4, M5, M6, M7, M10, M12, M13, M14, M15, M16, M17, M19, M21, C22, x, false, beta);

    // Handle edge cases with dynamic peeling
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
    if (total_steps == steps_left) {
        mkl_set_num_threads_local(num_threads);
        mkl_set_dynamic(0);
    }
#endif
    DynamicPeeling(A, B, C, 2, 7, 2, beta);
}

// C := alpha * A * B + beta * C
template <typename Scalar>
double FastMatmul(Matrix<Scalar>& A, Matrix<Scalar>& B, Matrix<Scalar>& C,
    int num_steps, double x=1e-8, int num_threads=12, Scalar alpha=Scalar(1.0), Scalar beta=Scalar(0.0)) {
    MemoryManager<Scalar> mem_mngr;
#ifdef _PARALLEL_
    mem_mngr.Allocate(2, 7, 2, 22, num_steps, A.m(), A.n(), B.n());
#endif
    A.set_multiplier(alpha);
    int num_multiplies_per_step = 22;
    int total_multiplies = pow(num_multiplies_per_step, num_steps);

    // Set parameters needed for all types of parallelism.
    // int num_threads = 0;
    omp_set_num_threads(num_threads);
#ifdef _PARALLEL_
# pragma omp parallel
    {
        if (omp_get_thread_num() == 0) { num_threads = omp_get_num_threads(); }
    }
    omp_set_max_active_levels(2);
#endif

#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_)
# pragma omp parallel
    {
        mkl_set_num_threads_local(1);
        mkl_set_dynamic(0);
    }
#endif

#if defined(_PARALLEL_) && (_PARALLEL_ == _DFS_PAR_)
    mkl_set_dynamic(0);
#endif

#if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    if (num_threads > total_multiplies) {
        mkl_set_dynamic(0);
    } else {
# pragma omp parallel
        {
            mkl_set_num_threads_local(1);
            mkl_set_dynamic(0);
        }
    }
#endif

    LockAndCounter locker(total_multiplies - (total_multiplies % num_threads));
    using FpMilliseconds = std::chrono::duration<float, std::chrono::milliseconds::period>;
    auto t1 = std::chrono::high_resolution_clock::now();

#ifdef _PARALLEL_
# pragma omp parallel
    {
# pragma omp single
#endif
        FastMatmulRecursive(locker, mem_mngr, A, B, C, num_steps, num_steps, 0, x, num_threads, beta);
#ifdef _PARALLEL_
    }
#endif
    auto t2 = std::chrono::high_resolution_clock::now();
    return FpMilliseconds(t2 - t1).count();
}

}  // namespace smirnov272_22_198_approx

#endif  // _smirnov272_22_198_approx_HPP_
