#ifndef _smirnov444_46_352_approx_HPP_
#define _smirnov444_46_352_approx_HPP_

// This is an automatically generated file from gen.py.
#include "common.hpp"

namespace smirnov444_46_352_approx {

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
            dataC[i + j * strideC] = Scalar(x) * dataS1[i + j * strideS1] + Scalar(x) * dataS2[i + j * strideS2] + Scalar(-(x)) * dataS3[i + j * strideS3] + Scalar(1.0 / (x)) * dataS4[i + j * strideS4] + Scalar(1.0 / (x)) * dataS5[i + j * strideS5];
        }
    }
}

template <typename Scalar>
void S_Add2(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = dataS1[i + j * strideS1] + Scalar(1.0 / (x)) * dataS2[i + j * strideS2] + Scalar(-(x)) * dataS3[i + j * strideS3];
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
            dataC[i + j * strideC] = Scalar(-(x * x)) * dataS1[i + j * strideS1] + Scalar(-(x)) * dataS2[i + j * strideS2] + Scalar(1.0 / (x)) * dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add4(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& S5, Matrix<Scalar>& S6, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideS4 = S4.stride();
    const int strideS5 = S5.stride();
    const int strideS6 = S6.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    const Scalar *dataS4 = S4.data();
    const Scalar *dataS5 = S5.data();
    const Scalar *dataS6 = S6.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataS1[i + j * strideS1] + Scalar(1.0 / (x)) * dataS2[i + j * strideS2] + Scalar(-(x)) * dataS3[i + j * strideS3] + Scalar(1.0 / (x)) * dataS4[i + j * strideS4] -dataS5[i + j * strideS5] + Scalar(1.0 / (x)) * dataS6[i + j * strideS6];
        }
    }
}

template <typename Scalar>
void S_Add5(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x * x) * dataS1[i + j * strideS1] + Scalar(x) * dataS2[i + j * strideS2] -dataS3[i + j * strideS3] + Scalar(1.0 / (x)) * dataS4[i + j * strideS4];
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
            dataC[i + j * strideC] = Scalar(x * x) * dataS1[i + j * strideS1] + Scalar(1.0 / (x)) * dataS2[i + j * strideS2] + dataS3[i + j * strideS3] + Scalar(x * x) * dataS4[i + j * strideS4];
        }
    }
}

template <typename Scalar>
void S_Add7(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x * x) * dataS1[i + j * strideS1] + Scalar(x * x) * dataS2[i + j * strideS2] + Scalar(x) * dataS3[i + j * strideS3] + Scalar(1.0 / (x)) * dataS4[i + j * strideS4];
        }
    }
}

template <typename Scalar>
void S_Add8(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1] + Scalar(x * x) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add9(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1];
        }
    }
}

template <typename Scalar>
void S_Add10(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& S5, Matrix<Scalar>& S6, Matrix<Scalar>& S7, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideS4 = S4.stride();
    const int strideS5 = S5.stride();
    const int strideS6 = S6.stride();
    const int strideS7 = S7.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    const Scalar *dataS4 = S4.data();
    const Scalar *dataS5 = S5.data();
    const Scalar *dataS6 = S6.data();
    const Scalar *dataS7 = S7.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(x) * dataS1[i + j * strideS1] + Scalar(x) * dataS2[i + j * strideS2] + Scalar(-(x)) * dataS3[i + j * strideS3] + Scalar(x) * dataS4[i + j * strideS4] + Scalar(1.0 / (x)) * dataS5[i + j * strideS5] + Scalar(1.0 / (x)) * dataS6[i + j * strideS6] + Scalar(1.0 / (x)) * dataS7[i + j * strideS7];
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
            dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataS1[i + j * strideS1] + Scalar(x) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add12(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x) * dataS1[i + j * strideS1] + Scalar(1.0 / (x)) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add13(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1] + Scalar(x) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add14(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = dataS1[i + j * strideS1] + Scalar(1.0 / (x)) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add15(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x * x) * dataS1[i + j * strideS1] + dataS2[i + j * strideS2] + Scalar(1.0 / (x)) * dataS3[i + j * strideS3];
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1] + Scalar(1.0 / (x)) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add17(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(x)) * dataS1[i + j * strideS1] + Scalar(1.0 / (x)) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add18(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x) * dataS1[i + j * strideS1] + Scalar(-(1.0 / (x))) * dataS2[i + j * strideS2];
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
            dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataS1[i + j * strideS1] + Scalar(x * x) * dataS2[i + j * strideS2] + Scalar(1.0 / (x)) * dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add20(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x * x) * dataS1[i + j * strideS1] + dataS2[i + j * strideS2] + Scalar(1.0 / (x)) * dataS3[i + j * strideS3];
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1] + dataS2[i + j * strideS2] + Scalar(-(1.0 / (x))) * dataS3[i + j * strideS3];
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
            dataC[i + j * strideC] = Scalar(-(x * x)) * dataS1[i + j * strideS1] + Scalar(1.0 / (x)) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add23(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1] -dataS2[i + j * strideS2] + Scalar(-(1.0 / (x))) * dataS3[i + j * strideS3] + Scalar(x * x) * dataS4[i + j * strideS4];
        }
    }
}

template <typename Scalar>
void S_Add24(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(x)) * dataS1[i + j * strideS1] + Scalar(x) * dataS2[i + j * strideS2] + Scalar(-(1.0 / (x))) * dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add25(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = dataS1[i + j * strideS1] + Scalar(-(1.0 / (x))) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add26(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& S5, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x * x) * dataS1[i + j * strideS1] + Scalar(1.0 / (x)) * dataS2[i + j * strideS2] + Scalar(-(1.0 / (x))) * dataS3[i + j * strideS3] + Scalar(-(x * x)) * dataS4[i + j * strideS4] + Scalar(x) * dataS5[i + j * strideS5];
        }
    }
}

template <typename Scalar>
void S_Add27(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x * x) * dataS1[i + j * strideS1] + Scalar(x) * dataS2[i + j * strideS2] + Scalar(1.0 / (x)) * dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add28(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x) * dataS1[i + j * strideS1] + Scalar(x) * dataS2[i + j * strideS2] + Scalar(1.0 / (x)) * dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add29(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x * x) * dataS1[i + j * strideS1] + Scalar(1.0 / (x)) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add30(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataS1[i + j * strideS1] + Scalar(x) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add31(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1];
        }
    }
}

template <typename Scalar>
void S_Add32(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1];
        }
    }
}

template <typename Scalar>
void S_Add33(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1];
        }
    }
}

template <typename Scalar>
void S_Add34(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& S5, Matrix<Scalar>& S6, Matrix<Scalar>& S7, Matrix<Scalar>& S8, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideS4 = S4.stride();
    const int strideS5 = S5.stride();
    const int strideS6 = S6.stride();
    const int strideS7 = S7.stride();
    const int strideS8 = S8.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    const Scalar *dataS4 = S4.data();
    const Scalar *dataS5 = S5.data();
    const Scalar *dataS6 = S6.data();
    const Scalar *dataS7 = S7.data();
    const Scalar *dataS8 = S8.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(x * x) * dataS1[i + j * strideS1] + dataS2[i + j * strideS2] + Scalar(-(1.0 / (x))) * dataS3[i + j * strideS3] + Scalar(1.0 / (x)) * dataS4[i + j * strideS4] + dataS5[i + j * strideS5] + Scalar(1.0 / (x)) * dataS6[i + j * strideS6] + Scalar(x * x) * dataS7[i + j * strideS7] + Scalar(-(x * x)) * dataS8[i + j * strideS8];
        }
    }
}

template <typename Scalar>
void S_Add35(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataS1[i + j * strideS1];
        }
    }
}

template <typename Scalar>
void S_Add36(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1] + Scalar(-(x * x)) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add37(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& S5, Matrix<Scalar>& S6, Matrix<Scalar>& S7, Matrix<Scalar>& S8, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideS4 = S4.stride();
    const int strideS5 = S5.stride();
    const int strideS6 = S6.stride();
    const int strideS7 = S7.stride();
    const int strideS8 = S8.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    const Scalar *dataS4 = S4.data();
    const Scalar *dataS5 = S5.data();
    const Scalar *dataS6 = S6.data();
    const Scalar *dataS7 = S7.data();
    const Scalar *dataS8 = S8.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(x * x) * dataS1[i + j * strideS1] + Scalar(x * x) * dataS2[i + j * strideS2] + Scalar(x) * dataS3[i + j * strideS3] + Scalar(-(x * x)) * dataS4[i + j * strideS4] + Scalar(-(x)) * dataS5[i + j * strideS5] + Scalar(1.0 / (x)) * dataS6[i + j * strideS6] + Scalar(1.0 / (x)) * dataS7[i + j * strideS7] + Scalar(1.0 / (x)) * dataS8[i + j * strideS8];
        }
    }
}

template <typename Scalar>
void S_Add38(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1];
        }
    }
}

template <typename Scalar>
void S_Add39(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataS1[i + j * strideS1];
        }
    }
}

template <typename Scalar>
void S_Add40(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& S5, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1] + Scalar(-(x * x)) * dataS2[i + j * strideS2] + Scalar(x * x) * dataS3[i + j * strideS3] + Scalar(1.0 / (x)) * dataS4[i + j * strideS4] + Scalar(1.0 / (x)) * dataS5[i + j * strideS5];
        }
    }
}

template <typename Scalar>
void S_Add41(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1];
        }
    }
}

template <typename Scalar>
void S_Add42(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1] + Scalar(x * x) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add43(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1] + Scalar(x * x) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add44(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1];
        }
    }
}

template <typename Scalar>
void S_Add45(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1];
        }
    }
}

template <typename Scalar>
void S_Add46(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataS1[i + j * strideS1];
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + dataT2[i + j * strideT2] + dataT3[i + j * strideT3] + Scalar(x) * dataT4[i + j * strideT4];
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
            dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataT1[i + j * strideT1] + Scalar(1.0 / (x)) * dataT2[i + j * strideT2] + dataT3[i + j * strideT3];
        }
    }
}

template <typename Scalar>
void T_Add3(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + Scalar(x) * dataT2[i + j * strideT2] + Scalar(x * x) * dataT3[i + j * strideT3];
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
            dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataT1[i + j * strideT1] + dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add5(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + dataT2[i + j * strideT2];
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
            dataC[i + j * strideC] = Scalar(x * x) * dataT1[i + j * strideT1] + Scalar(1.0 / (x)) * dataT2[i + j * strideT2];
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + Scalar(x * x) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add8(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataT1[i + j * strideT1] + Scalar(x * x) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add9(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] -dataT2[i + j * strideT2] + Scalar(-(x * x)) * dataT3[i + j * strideT3];
        }
    }
}

template <typename Scalar>
void T_Add10(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = -dataT1[i + j * strideT1] + Scalar(x * x) * dataT2[i + j * strideT2];
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add12(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x) * dataT1[i + j * strideT1] + dataT2[i + j * strideT2] + Scalar(x * x) * dataT3[i + j * strideT3];
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
            dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataT1[i + j * strideT1] + Scalar(x * x) * dataT2[i + j * strideT2] + Scalar(x * x) * dataT3[i + j * strideT3];
        }
    }
}

template <typename Scalar>
void T_Add14(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataT1[i + j * strideT1];
        }
    }
}

template <typename Scalar>
void T_Add15(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + Scalar(-(1.0 / (x))) * dataT2[i + j * strideT2] + Scalar(x) * dataT3[i + j * strideT3] + Scalar(x * x) * dataT4[i + j * strideT4];
        }
    }
}

template <typename Scalar>
void T_Add16(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = -dataT1[i + j * strideT1] + Scalar(x * x) * dataT2[i + j * strideT2] + Scalar(x) * dataT3[i + j * strideT3];
        }
    }
}

template <typename Scalar>
void T_Add17(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + dataT2[i + j * strideT2] + Scalar(1.0 / (x)) * dataT3[i + j * strideT3];
        }
    }
}

template <typename Scalar>
void T_Add18(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x * x) * dataT1[i + j * strideT1] + dataT2[i + j * strideT2] + Scalar(-(x * x)) * dataT3[i + j * strideT3];
        }
    }
}

template <typename Scalar>
void T_Add19(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + Scalar(x * x) * dataT2[i + j * strideT2];
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
            dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataT1[i + j * strideT1] + Scalar(1.0 / (x)) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add21(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + Scalar(-(x * x)) * dataT2[i + j * strideT2];
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
            dataC[i + j * strideC] = Scalar(x * x) * dataT1[i + j * strideT1] + Scalar(1.0 / (x)) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add23(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + dataT2[i + j * strideT2] + Scalar(-(x)) * dataT3[i + j * strideT3];
        }
    }
}

template <typename Scalar>
void T_Add24(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + dataT2[i + j * strideT2] + Scalar(-(1.0 / (x))) * dataT3[i + j * strideT3] + Scalar(x) * dataT4[i + j * strideT4];
        }
    }
}

template <typename Scalar>
void T_Add25(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& T5, Matrix<Scalar>& T6, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + Scalar(1.0 / (x)) * dataT2[i + j * strideT2] + Scalar(1.0 / (x)) * dataT3[i + j * strideT3] + Scalar(-(x * x)) * dataT4[i + j * strideT4] + dataT5[i + j * strideT5] -dataT6[i + j * strideT6];
        }
    }
}

template <typename Scalar>
void T_Add26(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + Scalar(-(x * x)) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add27(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] -dataT2[i + j * strideT2] + Scalar(x * x) * dataT3[i + j * strideT3];
        }
    }
}

template <typename Scalar>
void T_Add28(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataT1[i + j * strideT1];
        }
    }
}

template <typename Scalar>
void T_Add29(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + Scalar(x * x) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add30(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(x)) * dataT1[i + j * strideT1] + Scalar(1.0 / (x)) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add31(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = -dataT1[i + j * strideT1] + Scalar(1.0 / (x)) * dataT2[i + j * strideT2] + Scalar(-(x * x)) * dataT3[i + j * strideT3] + dataT4[i + j * strideT4];
        }
    }
}

template <typename Scalar>
void T_Add32(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& T5, Matrix<Scalar>& T6, Matrix<Scalar>& T7, Matrix<Scalar>& T8, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideT3 = T3.stride();
    const int strideT4 = T4.stride();
    const int strideT5 = T5.stride();
    const int strideT6 = T6.stride();
    const int strideT7 = T7.stride();
    const int strideT8 = T8.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    const Scalar *dataT3 = T3.data();
    const Scalar *dataT4 = T4.data();
    const Scalar *dataT5 = T5.data();
    const Scalar *dataT6 = T6.data();
    const Scalar *dataT7 = T7.data();
    const Scalar *dataT8 = T8.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataT1[i + j * strideT1] + dataT2[i + j * strideT2] + Scalar(1.0 / (x)) * dataT3[i + j * strideT3] + dataT4[i + j * strideT4] + Scalar(-(1.0 / (x))) * dataT5[i + j * strideT5] + Scalar(-(x * x)) * dataT6[i + j * strideT6] + Scalar(-(x)) * dataT7[i + j * strideT7] + Scalar(x * x) * dataT8[i + j * strideT8];
        }
    }
}

template <typename Scalar>
void T_Add33(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x) * dataT1[i + j * strideT1] + Scalar(-(1.0 / (x))) * dataT2[i + j * strideT2] -dataT3[i + j * strideT3] + Scalar(x) * dataT4[i + j * strideT4];
        }
    }
}

template <typename Scalar>
void T_Add34(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1];
        }
    }
}

template <typename Scalar>
void T_Add35(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& T5, Matrix<Scalar>& T6, Matrix<Scalar>& T7, Matrix<Scalar>& T8, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideT3 = T3.stride();
    const int strideT4 = T4.stride();
    const int strideT5 = T5.stride();
    const int strideT6 = T6.stride();
    const int strideT7 = T7.stride();
    const int strideT8 = T8.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    const Scalar *dataT3 = T3.data();
    const Scalar *dataT4 = T4.data();
    const Scalar *dataT5 = T5.data();
    const Scalar *dataT6 = T6.data();
    const Scalar *dataT7 = T7.data();
    const Scalar *dataT8 = T8.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataT1[i + j * strideT1] + Scalar(-(x)) * dataT2[i + j * strideT2] + Scalar(x) * dataT3[i + j * strideT3] + Scalar(1.0 / (x)) * dataT4[i + j * strideT4] + dataT5[i + j * strideT5] + Scalar(x * x) * dataT6[i + j * strideT6] + Scalar(-(x * x)) * dataT7[i + j * strideT7] + Scalar(x * x) * dataT8[i + j * strideT8];
        }
    }
}

template <typename Scalar>
void T_Add36(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + Scalar(-(x)) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add37(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataT1[i + j * strideT1];
        }
    }
}

template <typename Scalar>
void T_Add38(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& T5, Matrix<Scalar>& T6, Matrix<Scalar>& T7, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideT3 = T3.stride();
    const int strideT4 = T4.stride();
    const int strideT5 = T5.stride();
    const int strideT6 = T6.stride();
    const int strideT7 = T7.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    const Scalar *dataT3 = T3.data();
    const Scalar *dataT4 = T4.data();
    const Scalar *dataT5 = T5.data();
    const Scalar *dataT6 = T6.data();
    const Scalar *dataT7 = T7.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataT1[i + j * strideT1] + Scalar(x * x) * dataT2[i + j * strideT2] + Scalar(1.0 / (x)) * dataT3[i + j * strideT3] + Scalar(-(x * x)) * dataT4[i + j * strideT4] + Scalar(x * x) * dataT5[i + j * strideT5] + dataT6[i + j * strideT6] + Scalar(-(x * x)) * dataT7[i + j * strideT7];
        }
    }
}

template <typename Scalar>
void T_Add39(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1];
        }
    }
}

template <typename Scalar>
void T_Add40(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataT1[i + j * strideT1];
        }
    }
}

template <typename Scalar>
void T_Add41(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& T5, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideT3 = T3.stride();
    const int strideT4 = T4.stride();
    const int strideT5 = T5.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    const Scalar *dataT3 = T3.data();
    const Scalar *dataT4 = T4.data();
    const Scalar *dataT5 = T5.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + Scalar(x * x) * dataT2[i + j * strideT2] + Scalar(1.0 / (x)) * dataT3[i + j * strideT3] + Scalar(-(x)) * dataT4[i + j * strideT4] + Scalar(-(1.0 / (x))) * dataT5[i + j * strideT5];
        }
    }
}

template <typename Scalar>
void T_Add42(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1];
        }
    }
}

template <typename Scalar>
void T_Add43(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1];
        }
    }
}

template <typename Scalar>
void T_Add44(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + Scalar(x * x) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add45(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + Scalar(x * x) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add46(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1];
        }
    }
}

template <typename Scalar>
void M_Add1(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
    const int strideM1 = M1.stride();
    const int strideM2 = M2.stride();
    const int strideM3 = M3.stride();
    const int strideM4 = M4.stride();
    const int strideC = C.stride();
    const Scalar *dataM1 = M1.data();
    const Scalar *dataM2 = M2.data();
    const Scalar *dataM3 = M3.data();
    const Scalar *dataM4 = M4.data();
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + Scalar(x * x) * dataM4[i + j * strideM4] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + Scalar(x * x) * dataM4[i + j * strideM4];
            }
        }
    }
}

template <typename Scalar>
void M_Add2(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
    const int strideM1 = M1.stride();
    const int strideM2 = M2.stride();
    const int strideM3 = M3.stride();
    const int strideM4 = M4.stride();
    const int strideM5 = M5.stride();
    const int strideM6 = M6.stride();
    const int strideC = C.stride();
    const Scalar *dataM1 = M1.data();
    const Scalar *dataM2 = M2.data();
    const Scalar *dataM3 = M3.data();
    const Scalar *dataM4 = M4.data();
    const Scalar *dataM5 = M5.data();
    const Scalar *dataM6 = M6.data();
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(x) * dataM1[i + j * strideM1] + dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(x * x) * dataM5[i + j * strideM5] + Scalar(-(1.0 / (x))) * dataM6[i + j * strideM6] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(x) * dataM1[i + j * strideM1] + dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(x * x) * dataM5[i + j * strideM5] + Scalar(-(1.0 / (x))) * dataM6[i + j * strideM6];
            }
        }
    }
}

template <typename Scalar>
void M_Add3(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
    const int strideM1 = M1.stride();
    const int strideM2 = M2.stride();
    const int strideM3 = M3.stride();
    const int strideM4 = M4.stride();
    const int strideM5 = M5.stride();
    const int strideM6 = M6.stride();
    const int strideM7 = M7.stride();
    const int strideC = C.stride();
    const Scalar *dataM1 = M1.data();
    const Scalar *dataM2 = M2.data();
    const Scalar *dataM3 = M3.data();
    const Scalar *dataM4 = M4.data();
    const Scalar *dataM5 = M5.data();
    const Scalar *dataM6 = M6.data();
    const Scalar *dataM7 = M7.data();
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = dataM1[i + j * strideM1] + Scalar(x) * dataM2[i + j * strideM2] + dataM3[i + j * strideM3] -dataM4[i + j * strideM4] + dataM5[i + j * strideM5] + Scalar(x * x) * dataM6[i + j * strideM6] -dataM7[i + j * strideM7] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = dataM1[i + j * strideM1] + Scalar(x) * dataM2[i + j * strideM2] + dataM3[i + j * strideM3] -dataM4[i + j * strideM4] + dataM5[i + j * strideM5] + Scalar(x * x) * dataM6[i + j * strideM6] -dataM7[i + j * strideM7];
            }
        }
    }
}

template <typename Scalar>
void M_Add4(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
    const int strideM1 = M1.stride();
    const int strideM2 = M2.stride();
    const int strideM3 = M3.stride();
    const int strideM4 = M4.stride();
    const int strideM5 = M5.stride();
    const int strideM6 = M6.stride();
    const int strideC = C.stride();
    const Scalar *dataM1 = M1.data();
    const Scalar *dataM2 = M2.data();
    const Scalar *dataM3 = M3.data();
    const Scalar *dataM4 = M4.data();
    const Scalar *dataM5 = M5.data();
    const Scalar *dataM6 = M6.data();
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(-(1.0 / (x))) * dataM3[i + j * strideM3] + dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(1.0 / (x)) * dataM6[i + j * strideM6] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(-(1.0 / (x))) * dataM3[i + j * strideM3] + dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(1.0 / (x)) * dataM6[i + j * strideM6];
            }
        }
    }
}

template <typename Scalar>
void M_Add5(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& M8, Matrix<Scalar>& M9, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
    const int strideM1 = M1.stride();
    const int strideM2 = M2.stride();
    const int strideM3 = M3.stride();
    const int strideM4 = M4.stride();
    const int strideM5 = M5.stride();
    const int strideM6 = M6.stride();
    const int strideM7 = M7.stride();
    const int strideM8 = M8.stride();
    const int strideM9 = M9.stride();
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
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(x) * dataM1[i + j * strideM1] + Scalar(x * x) * dataM2[i + j * strideM2] + Scalar(x) * dataM3[i + j * strideM3] + Scalar(-(x)) * dataM4[i + j * strideM4] -dataM5[i + j * strideM5] + Scalar(x) * dataM6[i + j * strideM6] + dataM7[i + j * strideM7] + dataM8[i + j * strideM8] -dataM9[i + j * strideM9] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(x) * dataM1[i + j * strideM1] + Scalar(x * x) * dataM2[i + j * strideM2] + Scalar(x) * dataM3[i + j * strideM3] + Scalar(-(x)) * dataM4[i + j * strideM4] -dataM5[i + j * strideM5] + Scalar(x) * dataM6[i + j * strideM6] + dataM7[i + j * strideM7] + dataM8[i + j * strideM8] -dataM9[i + j * strideM9];
            }
        }
    }
}

template <typename Scalar>
void M_Add6(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
    const int strideM1 = M1.stride();
    const int strideM2 = M2.stride();
    const int strideM3 = M3.stride();
    const int strideM4 = M4.stride();
    const int strideM5 = M5.stride();
    const int strideM6 = M6.stride();
    const int strideC = C.stride();
    const Scalar *dataM1 = M1.data();
    const Scalar *dataM2 = M2.data();
    const Scalar *dataM3 = M3.data();
    const Scalar *dataM4 = M4.data();
    const Scalar *dataM5 = M5.data();
    const Scalar *dataM6 = M6.data();
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(x) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(x) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(-(1.0 / (x))) * dataM5[i + j * strideM5] + Scalar(1.0 / (x)) * dataM6[i + j * strideM6] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(x) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(x) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(-(1.0 / (x))) * dataM5[i + j * strideM5] + Scalar(1.0 / (x)) * dataM6[i + j * strideM6];
            }
        }
    }
}

template <typename Scalar>
void M_Add7(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
    const int strideM1 = M1.stride();
    const int strideM2 = M2.stride();
    const int strideM3 = M3.stride();
    const int strideM4 = M4.stride();
    const int strideM5 = M5.stride();
    const int strideC = C.stride();
    const Scalar *dataM1 = M1.data();
    const Scalar *dataM2 = M2.data();
    const Scalar *dataM3 = M3.data();
    const Scalar *dataM4 = M4.data();
    const Scalar *dataM5 = M5.data();
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(x) * dataM1[i + j * strideM1] + Scalar(x * x) * dataM2[i + j * strideM2] + Scalar(-(x * x)) * dataM3[i + j * strideM3] + Scalar(x * x) * dataM4[i + j * strideM4] + Scalar((x + x * x)) * dataM5[i + j * strideM5] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(x) * dataM1[i + j * strideM1] + Scalar(x * x) * dataM2[i + j * strideM2] + Scalar(-(x * x)) * dataM3[i + j * strideM3] + Scalar(x * x) * dataM4[i + j * strideM4] + Scalar((x + x * x)) * dataM5[i + j * strideM5];
            }
        }
    }
}

template <typename Scalar>
void M_Add8(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
    const int strideM1 = M1.stride();
    const int strideM2 = M2.stride();
    const int strideM3 = M3.stride();
    const int strideM4 = M4.stride();
    const int strideM5 = M5.stride();
    const int strideC = C.stride();
    const Scalar *dataM1 = M1.data();
    const Scalar *dataM2 = M2.data();
    const Scalar *dataM3 = M3.data();
    const Scalar *dataM4 = M4.data();
    const Scalar *dataM5 = M5.data();
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5];
            }
        }
    }
}

template <typename Scalar>
void M_Add9(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
    const int strideM1 = M1.stride();
    const int strideM2 = M2.stride();
    const int strideM3 = M3.stride();
    const int strideM4 = M4.stride();
    const int strideM5 = M5.stride();
    const int strideM6 = M6.stride();
    const int strideM7 = M7.stride();
    const int strideC = C.stride();
    const Scalar *dataM1 = M1.data();
    const Scalar *dataM2 = M2.data();
    const Scalar *dataM3 = M3.data();
    const Scalar *dataM4 = M4.data();
    const Scalar *dataM5 = M5.data();
    const Scalar *dataM6 = M6.data();
    const Scalar *dataM7 = M7.data();
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataM1[i + j * strideM1] + dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(x) * dataM6[i + j * strideM6] + Scalar(1.0 / (x)) * dataM7[i + j * strideM7] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataM1[i + j * strideM1] + dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(x) * dataM6[i + j * strideM6] + Scalar(1.0 / (x)) * dataM7[i + j * strideM7];
            }
        }
    }
}

template <typename Scalar>
void M_Add10(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
    const int strideM1 = M1.stride();
    const int strideM2 = M2.stride();
    const int strideM3 = M3.stride();
    const int strideM4 = M4.stride();
    const int strideM5 = M5.stride();
    const int strideC = C.stride();
    const Scalar *dataM1 = M1.data();
    const Scalar *dataM2 = M2.data();
    const Scalar *dataM3 = M3.data();
    const Scalar *dataM4 = M4.data();
    const Scalar *dataM5 = M5.data();
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataM1[i + j * strideM1] + Scalar(-(x)) * dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataM1[i + j * strideM1] + Scalar(-(x)) * dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5];
            }
        }
    }
}

template <typename Scalar>
void M_Add11(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
    const int strideM1 = M1.stride();
    const int strideM2 = M2.stride();
    const int strideM3 = M3.stride();
    const int strideM4 = M4.stride();
    const int strideM5 = M5.stride();
    const int strideC = C.stride();
    const Scalar *dataM1 = M1.data();
    const Scalar *dataM2 = M2.data();
    const Scalar *dataM3 = M3.data();
    const Scalar *dataM4 = M4.data();
    const Scalar *dataM5 = M5.data();
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(-(1.0 / (x))) * dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(-(1.0 / (x))) * dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5];
            }
        }
    }
}

template <typename Scalar>
void M_Add12(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& M8, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
    const int strideM1 = M1.stride();
    const int strideM2 = M2.stride();
    const int strideM3 = M3.stride();
    const int strideM4 = M4.stride();
    const int strideM5 = M5.stride();
    const int strideM6 = M6.stride();
    const int strideM7 = M7.stride();
    const int strideM8 = M8.stride();
    const int strideC = C.stride();
    const Scalar *dataM1 = M1.data();
    const Scalar *dataM2 = M2.data();
    const Scalar *dataM3 = M3.data();
    const Scalar *dataM4 = M4.data();
    const Scalar *dataM5 = M5.data();
    const Scalar *dataM6 = M6.data();
    const Scalar *dataM7 = M7.data();
    const Scalar *dataM8 = M8.data();
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(-(1.0 / (x))) * dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + Scalar(-(1.0 / (x))) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(-(1.0 / (x))) * dataM6[i + j * strideM6] + Scalar(1.0 / (x)) * dataM7[i + j * strideM7] + Scalar(-(1.0 / (x))) * dataM8[i + j * strideM8] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(-(1.0 / (x))) * dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + Scalar(-(1.0 / (x))) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(-(1.0 / (x))) * dataM6[i + j * strideM6] + Scalar(1.0 / (x)) * dataM7[i + j * strideM7] + Scalar(-(1.0 / (x))) * dataM8[i + j * strideM8];
            }
        }
    }
}

template <typename Scalar>
void M_Add13(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
    const int strideM1 = M1.stride();
    const int strideM2 = M2.stride();
    const int strideM3 = M3.stride();
    const int strideM4 = M4.stride();
    const int strideM5 = M5.stride();
    const int strideM6 = M6.stride();
    const int strideC = C.stride();
    const Scalar *dataM1 = M1.data();
    const Scalar *dataM2 = M2.data();
    const Scalar *dataM3 = M3.data();
    const Scalar *dataM4 = M4.data();
    const Scalar *dataM5 = M5.data();
    const Scalar *dataM6 = M6.data();
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = -dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(1.0 / (x)) * dataM6[i + j * strideM6] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = -dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(1.0 / (x)) * dataM6[i + j * strideM6];
            }
        }
    }
}

template <typename Scalar>
void M_Add14(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
    const int strideM1 = M1.stride();
    const int strideM2 = M2.stride();
    const int strideM3 = M3.stride();
    const int strideM4 = M4.stride();
    const int strideM5 = M5.stride();
    const int strideM6 = M6.stride();
    const int strideM7 = M7.stride();
    const int strideC = C.stride();
    const Scalar *dataM1 = M1.data();
    const Scalar *dataM2 = M2.data();
    const Scalar *dataM3 = M3.data();
    const Scalar *dataM4 = M4.data();
    const Scalar *dataM5 = M5.data();
    const Scalar *dataM6 = M6.data();
    const Scalar *dataM7 = M7.data();
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(-(x)) * dataM1[i + j * strideM1] + Scalar(x) * dataM2[i + j * strideM2] + Scalar(x) * dataM3[i + j * strideM3] + Scalar(x) * dataM4[i + j * strideM4] + Scalar(-(x)) * dataM5[i + j * strideM5] + Scalar(x * x) * dataM6[i + j * strideM6] + Scalar(x) * dataM7[i + j * strideM7] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(-(x)) * dataM1[i + j * strideM1] + Scalar(x) * dataM2[i + j * strideM2] + Scalar(x) * dataM3[i + j * strideM3] + Scalar(x) * dataM4[i + j * strideM4] + Scalar(-(x)) * dataM5[i + j * strideM5] + Scalar(x * x) * dataM6[i + j * strideM6] + Scalar(x) * dataM7[i + j * strideM7];
            }
        }
    }
}

template <typename Scalar>
void M_Add15(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
    const int strideM1 = M1.stride();
    const int strideM2 = M2.stride();
    const int strideM3 = M3.stride();
    const int strideM4 = M4.stride();
    const int strideM5 = M5.stride();
    const int strideM6 = M6.stride();
    const int strideC = C.stride();
    const Scalar *dataM1 = M1.data();
    const Scalar *dataM2 = M2.data();
    const Scalar *dataM3 = M3.data();
    const Scalar *dataM4 = M4.data();
    const Scalar *dataM5 = M5.data();
    const Scalar *dataM6 = M6.data();
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(x) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(x) * dataM5[i + j * strideM5] + Scalar(1.0 / (x)) * dataM6[i + j * strideM6] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(x) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(x) * dataM5[i + j * strideM5] + Scalar(1.0 / (x)) * dataM6[i + j * strideM6];
            }
        }
    }
}

template <typename Scalar>
void M_Add16(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
    const int strideM1 = M1.stride();
    const int strideM2 = M2.stride();
    const int strideM3 = M3.stride();
    const int strideM4 = M4.stride();
    const int strideM5 = M5.stride();
    const int strideC = C.stride();
    const Scalar *dataM1 = M1.data();
    const Scalar *dataM2 = M2.data();
    const Scalar *dataM3 = M3.data();
    const Scalar *dataM4 = M4.data();
    const Scalar *dataM5 = M5.data();
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(-(x)) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + Scalar(-(1.0 / (x))) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(-(x)) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + Scalar(-(1.0 / (x))) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5];
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

    Matrix<Scalar> A11 = A.Subblock(4, 4, 1, 1);
    Matrix<Scalar> A12 = A.Subblock(4, 4, 1, 2);
    Matrix<Scalar> A13 = A.Subblock(4, 4, 1, 3);
    Matrix<Scalar> A14 = A.Subblock(4, 4, 1, 4);
    Matrix<Scalar> A21 = A.Subblock(4, 4, 2, 1);
    Matrix<Scalar> A22 = A.Subblock(4, 4, 2, 2);
    Matrix<Scalar> A23 = A.Subblock(4, 4, 2, 3);
    Matrix<Scalar> A24 = A.Subblock(4, 4, 2, 4);
    Matrix<Scalar> A31 = A.Subblock(4, 4, 3, 1);
    Matrix<Scalar> A32 = A.Subblock(4, 4, 3, 2);
    Matrix<Scalar> A33 = A.Subblock(4, 4, 3, 3);
    Matrix<Scalar> A34 = A.Subblock(4, 4, 3, 4);
    Matrix<Scalar> A41 = A.Subblock(4, 4, 4, 1);
    Matrix<Scalar> A42 = A.Subblock(4, 4, 4, 2);
    Matrix<Scalar> A43 = A.Subblock(4, 4, 4, 3);
    Matrix<Scalar> A44 = A.Subblock(4, 4, 4, 4);
    Matrix<Scalar> B11 = B.Subblock(4, 4, 1, 1);
    Matrix<Scalar> B12 = B.Subblock(4, 4, 1, 2);
    Matrix<Scalar> B13 = B.Subblock(4, 4, 1, 3);
    Matrix<Scalar> B14 = B.Subblock(4, 4, 1, 4);
    Matrix<Scalar> B21 = B.Subblock(4, 4, 2, 1);
    Matrix<Scalar> B22 = B.Subblock(4, 4, 2, 2);
    Matrix<Scalar> B23 = B.Subblock(4, 4, 2, 3);
    Matrix<Scalar> B24 = B.Subblock(4, 4, 2, 4);
    Matrix<Scalar> B31 = B.Subblock(4, 4, 3, 1);
    Matrix<Scalar> B32 = B.Subblock(4, 4, 3, 2);
    Matrix<Scalar> B33 = B.Subblock(4, 4, 3, 3);
    Matrix<Scalar> B34 = B.Subblock(4, 4, 3, 4);
    Matrix<Scalar> B41 = B.Subblock(4, 4, 4, 1);
    Matrix<Scalar> B42 = B.Subblock(4, 4, 4, 2);
    Matrix<Scalar> B43 = B.Subblock(4, 4, 4, 3);
    Matrix<Scalar> B44 = B.Subblock(4, 4, 4, 4);
    Matrix<Scalar> C11 = C.Subblock(4, 4, 1, 1);
    Matrix<Scalar> C12 = C.Subblock(4, 4, 1, 2);
    Matrix<Scalar> C13 = C.Subblock(4, 4, 1, 3);
    Matrix<Scalar> C14 = C.Subblock(4, 4, 1, 4);
    Matrix<Scalar> C21 = C.Subblock(4, 4, 2, 1);
    Matrix<Scalar> C22 = C.Subblock(4, 4, 2, 2);
    Matrix<Scalar> C23 = C.Subblock(4, 4, 2, 3);
    Matrix<Scalar> C24 = C.Subblock(4, 4, 2, 4);
    Matrix<Scalar> C31 = C.Subblock(4, 4, 3, 1);
    Matrix<Scalar> C32 = C.Subblock(4, 4, 3, 2);
    Matrix<Scalar> C33 = C.Subblock(4, 4, 3, 3);
    Matrix<Scalar> C34 = C.Subblock(4, 4, 3, 4);
    Matrix<Scalar> C41 = C.Subblock(4, 4, 4, 1);
    Matrix<Scalar> C42 = C.Subblock(4, 4, 4, 2);
    Matrix<Scalar> C43 = C.Subblock(4, 4, 4, 3);
    Matrix<Scalar> C44 = C.Subblock(4, 4, 4, 4);


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
    Matrix<Scalar> M23(mem_mngr.GetMem(start_index, 23, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M24(mem_mngr.GetMem(start_index, 24, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M25(mem_mngr.GetMem(start_index, 25, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M26(mem_mngr.GetMem(start_index, 26, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M27(mem_mngr.GetMem(start_index, 27, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M28(mem_mngr.GetMem(start_index, 28, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M29(mem_mngr.GetMem(start_index, 29, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M30(mem_mngr.GetMem(start_index, 30, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M31(mem_mngr.GetMem(start_index, 31, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M32(mem_mngr.GetMem(start_index, 32, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M33(mem_mngr.GetMem(start_index, 33, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M34(mem_mngr.GetMem(start_index, 34, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M35(mem_mngr.GetMem(start_index, 35, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M36(mem_mngr.GetMem(start_index, 36, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M37(mem_mngr.GetMem(start_index, 37, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M38(mem_mngr.GetMem(start_index, 38, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M39(mem_mngr.GetMem(start_index, 39, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M40(mem_mngr.GetMem(start_index, 40, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M41(mem_mngr.GetMem(start_index, 41, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M42(mem_mngr.GetMem(start_index, 42, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M43(mem_mngr.GetMem(start_index, 43, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M44(mem_mngr.GetMem(start_index, 44, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M45(mem_mngr.GetMem(start_index, 45, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M46(mem_mngr.GetMem(start_index, 46, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
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
    Matrix<Scalar> M23(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M24(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M25(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M26(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M27(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M28(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M29(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M30(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M31(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M32(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M33(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M34(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M35(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M36(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M37(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M38(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M39(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M40(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M41(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M42(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M43(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M44(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M45(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M46(C11.m(), C11.n(), C.multiplier());
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
    bool sequential1 = should_launch_task(46, total_steps, steps_left, start_index, 1, num_threads);
    bool sequential2 = should_launch_task(46, total_steps, steps_left, start_index, 2, num_threads);
    bool sequential3 = should_launch_task(46, total_steps, steps_left, start_index, 3, num_threads);
    bool sequential4 = should_launch_task(46, total_steps, steps_left, start_index, 4, num_threads);
    bool sequential5 = should_launch_task(46, total_steps, steps_left, start_index, 5, num_threads);
    bool sequential6 = should_launch_task(46, total_steps, steps_left, start_index, 6, num_threads);
    bool sequential7 = should_launch_task(46, total_steps, steps_left, start_index, 7, num_threads);
    bool sequential8 = should_launch_task(46, total_steps, steps_left, start_index, 8, num_threads);
    bool sequential9 = should_launch_task(46, total_steps, steps_left, start_index, 9, num_threads);
    bool sequential10 = should_launch_task(46, total_steps, steps_left, start_index, 10, num_threads);
    bool sequential11 = should_launch_task(46, total_steps, steps_left, start_index, 11, num_threads);
    bool sequential12 = should_launch_task(46, total_steps, steps_left, start_index, 12, num_threads);
    bool sequential13 = should_launch_task(46, total_steps, steps_left, start_index, 13, num_threads);
    bool sequential14 = should_launch_task(46, total_steps, steps_left, start_index, 14, num_threads);
    bool sequential15 = should_launch_task(46, total_steps, steps_left, start_index, 15, num_threads);
    bool sequential16 = should_launch_task(46, total_steps, steps_left, start_index, 16, num_threads);
    bool sequential17 = should_launch_task(46, total_steps, steps_left, start_index, 17, num_threads);
    bool sequential18 = should_launch_task(46, total_steps, steps_left, start_index, 18, num_threads);
    bool sequential19 = should_launch_task(46, total_steps, steps_left, start_index, 19, num_threads);
    bool sequential20 = should_launch_task(46, total_steps, steps_left, start_index, 20, num_threads);
    bool sequential21 = should_launch_task(46, total_steps, steps_left, start_index, 21, num_threads);
    bool sequential22 = should_launch_task(46, total_steps, steps_left, start_index, 22, num_threads);
    bool sequential23 = should_launch_task(46, total_steps, steps_left, start_index, 23, num_threads);
    bool sequential24 = should_launch_task(46, total_steps, steps_left, start_index, 24, num_threads);
    bool sequential25 = should_launch_task(46, total_steps, steps_left, start_index, 25, num_threads);
    bool sequential26 = should_launch_task(46, total_steps, steps_left, start_index, 26, num_threads);
    bool sequential27 = should_launch_task(46, total_steps, steps_left, start_index, 27, num_threads);
    bool sequential28 = should_launch_task(46, total_steps, steps_left, start_index, 28, num_threads);
    bool sequential29 = should_launch_task(46, total_steps, steps_left, start_index, 29, num_threads);
    bool sequential30 = should_launch_task(46, total_steps, steps_left, start_index, 30, num_threads);
    bool sequential31 = should_launch_task(46, total_steps, steps_left, start_index, 31, num_threads);
    bool sequential32 = should_launch_task(46, total_steps, steps_left, start_index, 32, num_threads);
    bool sequential33 = should_launch_task(46, total_steps, steps_left, start_index, 33, num_threads);
    bool sequential34 = should_launch_task(46, total_steps, steps_left, start_index, 34, num_threads);
    bool sequential35 = should_launch_task(46, total_steps, steps_left, start_index, 35, num_threads);
    bool sequential36 = should_launch_task(46, total_steps, steps_left, start_index, 36, num_threads);
    bool sequential37 = should_launch_task(46, total_steps, steps_left, start_index, 37, num_threads);
    bool sequential38 = should_launch_task(46, total_steps, steps_left, start_index, 38, num_threads);
    bool sequential39 = should_launch_task(46, total_steps, steps_left, start_index, 39, num_threads);
    bool sequential40 = should_launch_task(46, total_steps, steps_left, start_index, 40, num_threads);
    bool sequential41 = should_launch_task(46, total_steps, steps_left, start_index, 41, num_threads);
    bool sequential42 = should_launch_task(46, total_steps, steps_left, start_index, 42, num_threads);
    bool sequential43 = should_launch_task(46, total_steps, steps_left, start_index, 43, num_threads);
    bool sequential44 = should_launch_task(46, total_steps, steps_left, start_index, 44, num_threads);
    bool sequential45 = should_launch_task(46, total_steps, steps_left, start_index, 45, num_threads);
    bool sequential46 = should_launch_task(46, total_steps, steps_left, start_index, 46, num_threads);
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
    bool sequential23 = false;
    bool sequential24 = false;
    bool sequential25 = false;
    bool sequential26 = false;
    bool sequential27 = false;
    bool sequential28 = false;
    bool sequential29 = false;
    bool sequential30 = false;
    bool sequential31 = false;
    bool sequential32 = false;
    bool sequential33 = false;
    bool sequential34 = false;
    bool sequential35 = false;
    bool sequential36 = false;
    bool sequential37 = false;
    bool sequential38 = false;
    bool sequential39 = false;
    bool sequential40 = false;
    bool sequential41 = false;
    bool sequential42 = false;
    bool sequential43 = false;
    bool sequential44 = false;
    bool sequential45 = false;
    bool sequential46 = false;
#endif



    // M1 = (x * A11 + x * A12 + -(x) * A13 + 1.0 / (x) * A14 + 1.0 / (x) * A22) * (1.0 / (x) * B13 + 1.0 * B21 + 1.0 * B31 + x * B43)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential1) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S1(mem_mngr.GetMem(start_index, 1, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S1(A11.m(), A11.n());
#endif
    S_Add1(A11, A12, A13, A14, A22, S1, x, sequential1);
#ifdef _PARALLEL_
    Matrix<Scalar> T1(mem_mngr.GetMem(start_index, 1, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T1(B11.m(), B11.n());
#endif
    T_Add1(B13, B21, B31, B43, T1, x, sequential1);
    FastMatmulRecursive(locker, mem_mngr, S1, T1, M1, total_steps, steps_left - 1, (start_index + 1 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S1.deallocate();
    T1.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 1, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M2 = (1.0 * A23 + 1.0 / (x) * A24 + -(x) * A33) * (-(1.0 / (x)) * B12 + 1.0 / (x) * B14 + 1.0 * B42)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential2) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S2(mem_mngr.GetMem(start_index, 2, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S2(A11.m(), A11.n());
#endif
    S_Add2(A23, A24, A33, S2, x, sequential2);
#ifdef _PARALLEL_
    Matrix<Scalar> T2(mem_mngr.GetMem(start_index, 2, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T2(B11.m(), B11.n());
#endif
    T_Add2(B12, B14, B42, T2, x, sequential2);
    FastMatmulRecursive(locker, mem_mngr, S2, T2, M2, total_steps, steps_left - 1, (start_index + 2 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S2.deallocate();
    T2.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 2, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M3 = (-(x * x) * A31 + -(x) * A34 + 1.0 / (x) * A44) * (1.0 / (x) * B12 + x * B13 + x * x * B43)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential3) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S3(mem_mngr.GetMem(start_index, 3, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S3(A11.m(), A11.n());
#endif
    S_Add3(A31, A34, A44, S3, x, sequential3);
#ifdef _PARALLEL_
    Matrix<Scalar> T3(mem_mngr.GetMem(start_index, 3, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T3(B11.m(), B11.n());
#endif
    T_Add3(B12, B13, B43, T3, x, sequential3);
    FastMatmulRecursive(locker, mem_mngr, S3, T3, M3, total_steps, steps_left - 1, (start_index + 3 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S3.deallocate();
    T3.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 3, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M4 = (1.0 * A23 + 1.0 / (x) * A24 + -(x) * A33 + 1.0 / (x) * A34 + -1.0 * A41 + 1.0 / (x) * A44) * (-(1.0 / (x)) * B12 + 1.0 * B42)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential4) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S4(mem_mngr.GetMem(start_index, 4, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S4(A11.m(), A11.n());
#endif
    S_Add4(A23, A24, A33, A34, A41, A44, S4, x, sequential4);
#ifdef _PARALLEL_
    Matrix<Scalar> T4(mem_mngr.GetMem(start_index, 4, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T4(B11.m(), B11.n());
#endif
    T_Add4(B12, B42, T4, x, sequential4);
    FastMatmulRecursive(locker, mem_mngr, S4, T4, M4, total_steps, steps_left - 1, (start_index + 4 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S4.deallocate();
    T4.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 4, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M5 = (x * x * A31 + x * A33 + -1.0 * A41 + 1.0 / (x) * A44) * (1.0 / (x) * B14 + 1.0 * B34)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential5) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S5(mem_mngr.GetMem(start_index, 5, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S5(A11.m(), A11.n());
#endif
    S_Add5(A31, A33, A41, A44, S5, x, sequential5);
#ifdef _PARALLEL_
    Matrix<Scalar> T5(mem_mngr.GetMem(start_index, 5, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T5(B11.m(), B11.n());
#endif
    T_Add5(B14, B34, T5, x, sequential5);
    FastMatmulRecursive(locker, mem_mngr, S5, T5, M5, total_steps, steps_left - 1, (start_index + 5 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S5.deallocate();
    T5.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 5, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M6 = (x * x * A21 + 1.0 / (x) * A22 + 1.0 * A23 + x * x * A33) * (x * x * B24 + 1.0 / (x) * B33)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential6) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S6(mem_mngr.GetMem(start_index, 6, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S6(A11.m(), A11.n());
#endif
    S_Add6(A21, A22, A23, A33, S6, x, sequential6);
#ifdef _PARALLEL_
    Matrix<Scalar> T6(mem_mngr.GetMem(start_index, 6, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T6(B11.m(), B11.n());
#endif
    T_Add6(B24, B33, T6, x, sequential6);
    FastMatmulRecursive(locker, mem_mngr, S6, T6, M6, total_steps, steps_left - 1, (start_index + 6 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S6.deallocate();
    T6.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 6, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M7 = (x * x * A21 + x * x * A23 + x * A24 + 1.0 / (x) * A42) * (1.0 / (x) * B12 + x * x * B23)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential7) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S7(mem_mngr.GetMem(start_index, 7, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S7(A11.m(), A11.n());
#endif
    S_Add7(A21, A23, A24, A42, S7, x, sequential7);
#ifdef _PARALLEL_
    Matrix<Scalar> T7(mem_mngr.GetMem(start_index, 7, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T7(B11.m(), B11.n());
#endif
    T_Add7(B12, B23, T7, x, sequential7);
    FastMatmulRecursive(locker, mem_mngr, S7, T7, M7, total_steps, steps_left - 1, (start_index + 7 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S7.deallocate();
    T7.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 7, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M8 = (1.0 / (x) * A14 + x * x * A33) * (-(1.0 / (x)) * B31 + x * x * B44)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential8) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S8(mem_mngr.GetMem(start_index, 8, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S8(A11.m(), A11.n());
#endif
    S_Add8(A14, A33, S8, x, sequential8);
#ifdef _PARALLEL_
    Matrix<Scalar> T8(mem_mngr.GetMem(start_index, 8, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T8(B11.m(), B11.n());
#endif
    T_Add8(B31, B44, T8, x, sequential8);
    FastMatmulRecursive(locker, mem_mngr, S8, T8, M8, total_steps, steps_left - 1, (start_index + 8 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S8.deallocate();
    T8.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 8, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M9 = (1.0 / (x) * A34) * (1.0 / (x) * B12 + -1.0 * B42 + -(x * x) * B44)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential9) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> T9(mem_mngr.GetMem(start_index, 9, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T9(B11.m(), B11.n());
#endif
    T_Add9(B12, B42, B44, T9, x, sequential9);
    M9.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, A34, T9, M9, total_steps, steps_left - 1, (start_index + 9 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    T9.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 9, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M10 = (x * A13 + x * A14 + -(x) * A23 + x * A33 + 1.0 / (x) * A42 + 1.0 / (x) * A43 + 1.0 / (x) * A44) * (-1.0 * B32 + x * x * B34)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential10) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S10(mem_mngr.GetMem(start_index, 10, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S10(A11.m(), A11.n());
#endif
    S_Add10(A13, A14, A23, A33, A42, A43, A44, S10, x, sequential10);
#ifdef _PARALLEL_
    Matrix<Scalar> T10(mem_mngr.GetMem(start_index, 10, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T10(B11.m(), B11.n());
#endif
    T_Add10(B32, B34, T10, x, sequential10);
    FastMatmulRecursive(locker, mem_mngr, S10, T10, M10, total_steps, steps_left - 1, (start_index + 10 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S10.deallocate();
    T10.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 10, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M11 = (-(1.0 / (x)) * A14 + x * A32) * (1.0 / (x) * B13 + 1.0 * B21)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential11) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S11(mem_mngr.GetMem(start_index, 11, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S11(A11.m(), A11.n());
#endif
    S_Add11(A14, A32, S11, x, sequential11);
#ifdef _PARALLEL_
    Matrix<Scalar> T11(mem_mngr.GetMem(start_index, 11, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T11(B11.m(), B11.n());
#endif
    T_Add11(B13, B21, T11, x, sequential11);
    FastMatmulRecursive(locker, mem_mngr, S11, T11, M11, total_steps, steps_left - 1, (start_index + 11 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S11.deallocate();
    T11.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 11, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M12 = (x * A33 + 1.0 / (x) * A44) * (x * B14 + 1.0 * B32 + x * x * B44)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential12) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S12(mem_mngr.GetMem(start_index, 12, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S12(A11.m(), A11.n());
#endif
    S_Add12(A33, A44, S12, x, sequential12);
#ifdef _PARALLEL_
    Matrix<Scalar> T12(mem_mngr.GetMem(start_index, 12, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T12(B11.m(), B11.n());
#endif
    T_Add12(B14, B32, B44, T12, x, sequential12);
    FastMatmulRecursive(locker, mem_mngr, S12, T12, M12, total_steps, steps_left - 1, (start_index + 12 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S12.deallocate();
    T12.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 12, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M13 = (1.0 / (x) * A12 + x * A41) * (-(1.0 / (x)) * B11 + x * x * B21 + x * x * B41)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential13) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S13(mem_mngr.GetMem(start_index, 13, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S13(A11.m(), A11.n());
#endif
    S_Add13(A12, A41, S13, x, sequential13);
#ifdef _PARALLEL_
    Matrix<Scalar> T13(mem_mngr.GetMem(start_index, 13, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T13(B11.m(), B11.n());
#endif
    T_Add13(B11, B21, B41, T13, x, sequential13);
    FastMatmulRecursive(locker, mem_mngr, S13, T13, M13, total_steps, steps_left - 1, (start_index + 13 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S13.deallocate();
    T13.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 13, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M14 = (1.0 * A23 + 1.0 / (x) * A24) * (-(1.0 / (x)) * B14)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential14) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S14(mem_mngr.GetMem(start_index, 14, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S14(A11.m(), A11.n());
#endif
    S_Add14(A23, A24, S14, x, sequential14);
    M14.UpdateMultiplier(Scalar(-(1.0 / (x))));
    FastMatmulRecursive(locker, mem_mngr, S14, B14, M14, total_steps, steps_left - 1, (start_index + 14 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S14.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 14, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M15 = (x * x * A21 + 1.0 * A23 + 1.0 / (x) * A24) * (1.0 / (x) * B14 + -(1.0 / (x)) * B33 + x * B34 + x * x * B44)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential15) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S15(mem_mngr.GetMem(start_index, 15, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S15(A11.m(), A11.n());
#endif
    S_Add15(A21, A23, A24, S15, x, sequential15);
#ifdef _PARALLEL_
    Matrix<Scalar> T15(mem_mngr.GetMem(start_index, 15, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T15(B11.m(), B11.n());
#endif
    T_Add15(B14, B33, B34, B44, T15, x, sequential15);
    FastMatmulRecursive(locker, mem_mngr, S15, T15, M15, total_steps, steps_left - 1, (start_index + 15 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S15.deallocate();
    T15.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 15, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M16 = (1.0 / (x) * A14 + 1.0 / (x) * A34) * (-1.0 * B12 + x * x * B41 + x * B42)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential16) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S16(mem_mngr.GetMem(start_index, 16, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S16(A11.m(), A11.n());
#endif
    S_Add16(A14, A34, S16, x, sequential16);
#ifdef _PARALLEL_
    Matrix<Scalar> T16(mem_mngr.GetMem(start_index, 16, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T16(B11.m(), B11.n());
#endif
    T_Add16(B12, B41, B42, T16, x, sequential16);
    FastMatmulRecursive(locker, mem_mngr, S16, T16, M16, total_steps, steps_left - 1, (start_index + 16 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S16.deallocate();
    T16.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 16, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M17 = (-(x) * A13 + 1.0 / (x) * A22) * (1.0 / (x) * B23 + 1.0 * B31 + 1.0 / (x) * B33)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential17) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S17(mem_mngr.GetMem(start_index, 17, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S17(A11.m(), A11.n());
#endif
    S_Add17(A13, A22, S17, x, sequential17);
#ifdef _PARALLEL_
    Matrix<Scalar> T17(mem_mngr.GetMem(start_index, 17, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T17(B11.m(), B11.n());
#endif
    T_Add17(B23, B31, B33, T17, x, sequential17);
    FastMatmulRecursive(locker, mem_mngr, S17, T17, M17, total_steps, steps_left - 1, (start_index + 17 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S17.deallocate();
    T17.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 17, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M18 = (x * A23 + -(1.0 / (x)) * A42) * (x * x * B24 + 1.0 * B32 + -(x * x) * B34)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential18) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S18(mem_mngr.GetMem(start_index, 18, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S18(A11.m(), A11.n());
#endif
    S_Add18(A23, A42, S18, x, sequential18);
#ifdef _PARALLEL_
    Matrix<Scalar> T18(mem_mngr.GetMem(start_index, 18, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T18(B11.m(), B11.n());
#endif
    T_Add18(B24, B32, B34, T18, x, sequential18);
    FastMatmulRecursive(locker, mem_mngr, S18, T18, M18, total_steps, steps_left - 1, (start_index + 18 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S18.deallocate();
    T18.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 18, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M19 = (-(1.0 / (x)) * A24 + x * x * A31 + 1.0 / (x) * A32) * (1.0 / (x) * B13 + x * x * B23)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential19) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S19(mem_mngr.GetMem(start_index, 19, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S19(A11.m(), A11.n());
#endif
    S_Add19(A24, A31, A32, S19, x, sequential19);
#ifdef _PARALLEL_
    Matrix<Scalar> T19(mem_mngr.GetMem(start_index, 19, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T19(B11.m(), B11.n());
#endif
    T_Add19(B13, B23, T19, x, sequential19);
    FastMatmulRecursive(locker, mem_mngr, S19, T19, M19, total_steps, steps_left - 1, (start_index + 19 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S19.deallocate();
    T19.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 19, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M20 = (x * x * A11 + 1.0 * A23 + 1.0 / (x) * A24) * (-(1.0 / (x)) * B14 + 1.0 / (x) * B31)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential20) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S20(mem_mngr.GetMem(start_index, 20, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S20(A11.m(), A11.n());
#endif
    S_Add20(A11, A23, A24, S20, x, sequential20);
#ifdef _PARALLEL_
    Matrix<Scalar> T20(mem_mngr.GetMem(start_index, 20, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T20(B11.m(), B11.n());
#endif
    T_Add20(B14, B31, T20, x, sequential20);
    FastMatmulRecursive(locker, mem_mngr, S20, T20, M20, total_steps, steps_left - 1, (start_index + 20 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S20.deallocate();
    T20.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 20, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M21 = (1.0 / (x) * A32 + 1.0 * A41 + -(1.0 / (x)) * A44) * (1.0 / (x) * B13 + -(x * x) * B24)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential21) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S21(mem_mngr.GetMem(start_index, 21, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S21(A11.m(), A11.n());
#endif
    S_Add21(A32, A41, A44, S21, x, sequential21);
#ifdef _PARALLEL_
    Matrix<Scalar> T21(mem_mngr.GetMem(start_index, 21, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T21(B11.m(), B11.n());
#endif
    T_Add21(B13, B24, T21, x, sequential21);
    FastMatmulRecursive(locker, mem_mngr, S21, T21, M21, total_steps, steps_left - 1, (start_index + 21 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S21.deallocate();
    T21.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 21, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M22 = (-(x * x) * A22 + 1.0 / (x) * A42) * (x * x * B21 + 1.0 / (x) * B22)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential22) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S22(mem_mngr.GetMem(start_index, 22, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S22(A11.m(), A11.n());
#endif
    S_Add22(A22, A42, S22, x, sequential22);
#ifdef _PARALLEL_
    Matrix<Scalar> T22(mem_mngr.GetMem(start_index, 22, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T22(B11.m(), B11.n());
#endif
    T_Add22(B21, B22, T22, x, sequential22);
    FastMatmulRecursive(locker, mem_mngr, S22, T22, M22, total_steps, steps_left - 1, (start_index + 22 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S22.deallocate();
    T22.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 22, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M23 = (1.0 / (x) * A21 + -1.0 * A22 + -(1.0 / (x)) * A24 + x * x * A31) * (1.0 / (x) * B13 + 1.0 * B31 + -(x) * B41)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential23) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S23(mem_mngr.GetMem(start_index, 23, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S23(A11.m(), A11.n());
#endif
    S_Add23(A21, A22, A24, A31, S23, x, sequential23);
#ifdef _PARALLEL_
    Matrix<Scalar> T23(mem_mngr.GetMem(start_index, 23, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T23(B11.m(), B11.n());
#endif
    T_Add23(B13, B31, B41, T23, x, sequential23);
    FastMatmulRecursive(locker, mem_mngr, S23, T23, M23, total_steps, steps_left - 1, (start_index + 23 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S23.deallocate();
    T23.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 23, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M24 = (-(x) * A12 + x * A13 + -(1.0 / (x)) * A22) * (1.0 / (x) * B13 + 1.0 * B21 + -(1.0 / (x)) * B23 + x * B43)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential24) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S24(mem_mngr.GetMem(start_index, 24, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S24(A11.m(), A11.n());
#endif
    S_Add24(A12, A13, A22, S24, x, sequential24);
#ifdef _PARALLEL_
    Matrix<Scalar> T24(mem_mngr.GetMem(start_index, 24, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T24(B11.m(), B11.n());
#endif
    T_Add24(B13, B21, B23, B43, T24, x, sequential24);
    FastMatmulRecursive(locker, mem_mngr, S24, T24, M24, total_steps, steps_left - 1, (start_index + 24 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S24.deallocate();
    T24.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 24, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M25 = (1.0 * A41 + -(1.0 / (x)) * A44) * (1.0 / (x) * B12 + 1.0 / (x) * B13 + 1.0 / (x) * B14 + -(x * x) * B24 + 1.0 * B34 + -1.0 * B42)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential25) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S25(mem_mngr.GetMem(start_index, 25, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S25(A11.m(), A11.n());
#endif
    S_Add25(A41, A44, S25, x, sequential25);
#ifdef _PARALLEL_
    Matrix<Scalar> T25(mem_mngr.GetMem(start_index, 25, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T25(B11.m(), B11.n());
#endif
    T_Add25(B12, B13, B14, B24, B34, B42, T25, x, sequential25);
    FastMatmulRecursive(locker, mem_mngr, S25, T25, M25, total_steps, steps_left - 1, (start_index + 25 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S25.deallocate();
    T25.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 25, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M26 = (x * x * A11 + 1.0 / (x) * A12 + -(1.0 / (x)) * A14 + -(x * x) * A31 + x * A41) * (1.0 / (x) * B11 + -(x * x) * B41)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential26) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S26(mem_mngr.GetMem(start_index, 26, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S26(A11.m(), A11.n());
#endif
    S_Add26(A11, A12, A14, A31, A41, S26, x, sequential26);
#ifdef _PARALLEL_
    Matrix<Scalar> T26(mem_mngr.GetMem(start_index, 26, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T26(B11.m(), B11.n());
#endif
    T_Add26(B11, B41, T26, x, sequential26);
    FastMatmulRecursive(locker, mem_mngr, S26, T26, M26, total_steps, steps_left - 1, (start_index + 26 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S26.deallocate();
    T26.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 26, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M27 = (x * x * A11 + x * A14 + 1.0 / (x) * A43) * (1.0 / (x) * B12 + -1.0 * B32 + x * x * B33)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential27) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S27(mem_mngr.GetMem(start_index, 27, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S27(A11.m(), A11.n());
#endif
    S_Add27(A11, A14, A43, S27, x, sequential27);
#ifdef _PARALLEL_
    Matrix<Scalar> T27(mem_mngr.GetMem(start_index, 27, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T27(B11.m(), B11.n());
#endif
    T_Add27(B12, B32, B33, T27, x, sequential27);
    FastMatmulRecursive(locker, mem_mngr, S27, T27, M27, total_steps, steps_left - 1, (start_index + 27 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S27.deallocate();
    T27.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 27, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M28 = (x * A13 + x * A14 + 1.0 / (x) * A43) * (1.0 * B32)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential28) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S28(mem_mngr.GetMem(start_index, 28, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S28(A11.m(), A11.n());
#endif
    S_Add28(A13, A14, A43, S28, x, sequential28);
    M28.UpdateMultiplier(Scalar(1.0));
    FastMatmulRecursive(locker, mem_mngr, S28, B32, M28, total_steps, steps_left - 1, (start_index + 28 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S28.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 28, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M29 = (x * x * A32 + 1.0 / (x) * A44) * (1.0 / (x) * B22 + x * x * B41)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential29) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S29(mem_mngr.GetMem(start_index, 29, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S29(A11.m(), A11.n());
#endif
    S_Add29(A32, A44, S29, x, sequential29);
#ifdef _PARALLEL_
    Matrix<Scalar> T29(mem_mngr.GetMem(start_index, 29, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T29(B11.m(), B11.n());
#endif
    T_Add29(B22, B41, T29, x, sequential29);
    FastMatmulRecursive(locker, mem_mngr, S29, T29, M29, total_steps, steps_left - 1, (start_index + 29 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S29.deallocate();
    T29.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 29, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M30 = (-(1.0 / (x)) * A12 + x * A43) * (-(x) * B24 + 1.0 / (x) * B31)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential30) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S30(mem_mngr.GetMem(start_index, 30, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S30(A11.m(), A11.n());
#endif
    S_Add30(A12, A43, S30, x, sequential30);
#ifdef _PARALLEL_
    Matrix<Scalar> T30(mem_mngr.GetMem(start_index, 30, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T30(B11.m(), B11.n());
#endif
    T_Add30(B24, B31, T30, x, sequential30);
    FastMatmulRecursive(locker, mem_mngr, S30, T30, M30, total_steps, steps_left - 1, (start_index + 30 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S30.deallocate();
    T30.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 30, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M31 = (1.0 / (x) * A12) * (-1.0 * B11 + 1.0 / (x) * B22 + -(x * x) * B24 + 1.0 * B31)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential31) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> T31(mem_mngr.GetMem(start_index, 31, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T31(B11.m(), B11.n());
#endif
    T_Add31(B11, B22, B24, B31, T31, x, sequential31);
    M31.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, A12, T31, M31, total_steps, steps_left - 1, (start_index + 31 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    T31.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 31, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M32 = (1.0 / (x) * A14) * (-(1.0 / (x)) * B11 + 1.0 * B12 + 1.0 / (x) * B13 + 1.0 * B21 + -(1.0 / (x)) * B31 + -(x * x) * B41 + -(x) * B42 + x * x * B44)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential32) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> T32(mem_mngr.GetMem(start_index, 32, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T32(B11.m(), B11.n());
#endif
    T_Add32(B11, B12, B13, B21, B31, B41, B42, B44, T32, x, sequential32);
    M32.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, A14, T32, M32, total_steps, steps_left - 1, (start_index + 32 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    T32.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 32, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M33 = (1.0 / (x) * A21) * (x * B11 + -(1.0 / (x)) * B13 + -1.0 * B31 + x * B41)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential33) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> T33(mem_mngr.GetMem(start_index, 33, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T33(B11.m(), B11.n());
#endif
    T_Add33(B11, B13, B31, B41, T33, x, sequential33);
    M33.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, A21, T33, M33, total_steps, steps_left - 1, (start_index + 33 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    T33.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 33, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M34 = (x * x * A11 + 1.0 * A12 + -(1.0 / (x)) * A13 + 1.0 / (x) * A14 + 1.0 * A23 + 1.0 / (x) * A24 + x * x * A33 + -(x * x) * A43) * (1.0 / (x) * B31)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential34) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S34(mem_mngr.GetMem(start_index, 34, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S34(A11.m(), A11.n());
#endif
    S_Add34(A11, A12, A13, A14, A23, A24, A33, A43, S34, x, sequential34);
    M34.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, S34, B31, M34, total_steps, steps_left - 1, (start_index + 34 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S34.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 34, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M35 = (-(1.0 / (x)) * A44) * (-(1.0 / (x)) * B12 + -(x) * B13 + x * B14 + 1.0 / (x) * B22 + 1.0 * B32 + x * x * B41 + -(x * x) * B43 + x * x * B44)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential35) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> T35(mem_mngr.GetMem(start_index, 35, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T35(B11.m(), B11.n());
#endif
    T_Add35(B12, B13, B14, B22, B32, B41, B43, B44, T35, x, sequential35);
    M35.UpdateMultiplier(Scalar(-(1.0 / (x))));
    FastMatmulRecursive(locker, mem_mngr, A44, T35, M35, total_steps, steps_left - 1, (start_index + 35 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    T35.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 35, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M36 = (1.0 / (x) * A24 + -(x * x) * A33) * (1.0 / (x) * B33 + -(x) * B34)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential36) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S36(mem_mngr.GetMem(start_index, 36, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S36(A11.m(), A11.n());
#endif
    S_Add36(A24, A33, S36, x, sequential36);
#ifdef _PARALLEL_
    Matrix<Scalar> T36(mem_mngr.GetMem(start_index, 36, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T36(B11.m(), B11.n());
#endif
    T_Add36(B33, B34, T36, x, sequential36);
    FastMatmulRecursive(locker, mem_mngr, S36, T36, M36, total_steps, steps_left - 1, (start_index + 36 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S36.deallocate();
    T36.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 36, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M37 = (x * x * A21 + x * x * A23 + x * A24 + -(x * x) * A31 + -(x) * A34 + 1.0 / (x) * A42 + 1.0 / (x) * A43 + 1.0 / (x) * A44) * (-(1.0 / (x)) * B12)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential37) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S37(mem_mngr.GetMem(start_index, 37, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S37(A11.m(), A11.n());
#endif
    S_Add37(A21, A23, A24, A31, A34, A42, A43, A44, S37, x, sequential37);
    M37.UpdateMultiplier(Scalar(-(1.0 / (x))));
    FastMatmulRecursive(locker, mem_mngr, S37, B12, M37, total_steps, steps_left - 1, (start_index + 37 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S37.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 37, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M38 = (1.0 / (x) * A42) * (-(1.0 / (x)) * B12 + x * x * B21 + 1.0 / (x) * B22 + -(x * x) * B23 + x * x * B24 + 1.0 * B32 + -(x * x) * B34)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential38) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> T38(mem_mngr.GetMem(start_index, 38, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T38(B11.m(), B11.n());
#endif
    T_Add38(B12, B21, B22, B23, B24, B32, B34, T38, x, sequential38);
    M38.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, A42, T38, M38, total_steps, steps_left - 1, (start_index + 38 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    T38.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 38, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M39 = (-(1.0 / (x)) * A32) * (1.0 / (x) * B13)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential39) shared(mem_mngr, locker) untied default(shared)
    {
#endif
    M39.UpdateMultiplier(Scalar(-(1.0 / (x))));
    M39.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, A32, B13, M39, total_steps, steps_left - 1, (start_index + 39 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 39, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M40 = (1.0 / (x) * A12 + -(x * x) * A22 + x * x * A32 + 1.0 / (x) * A42 + 1.0 / (x) * A44) * (-(1.0 / (x)) * B22)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential40) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S40(mem_mngr.GetMem(start_index, 40, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S40(A11.m(), A11.n());
#endif
    S_Add40(A12, A22, A32, A42, A44, S40, x, sequential40);
    M40.UpdateMultiplier(Scalar(-(1.0 / (x))));
    FastMatmulRecursive(locker, mem_mngr, S40, B22, M40, total_steps, steps_left - 1, (start_index + 40 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S40.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 40, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M41 = (1.0 / (x) * A24) * (1.0 / (x) * B13 + x * x * B23 + 1.0 / (x) * B33 + -(x) * B34 + -(1.0 / (x)) * B43)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential41) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> T41(mem_mngr.GetMem(start_index, 41, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T41(B11.m(), B11.n());
#endif
    T_Add41(B13, B23, B33, B34, B43, T41, x, sequential41);
    M41.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, A24, T41, M41, total_steps, steps_left - 1, (start_index + 41 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    T41.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 41, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M42 = (1.0 / (x) * A24 + x * x * A34) * (1.0 / (x) * B43)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential42) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S42(mem_mngr.GetMem(start_index, 42, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S42(A11.m(), A11.n());
#endif
    S_Add42(A24, A34, S42, x, sequential42);
    M42.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, S42, B43, M42, total_steps, steps_left - 1, (start_index + 42 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S42.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 42, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M43 = (1.0 / (x) * A14 + x * x * A31) * (1.0 / (x) * B11)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential43) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S43(mem_mngr.GetMem(start_index, 43, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S43(A11.m(), A11.n());
#endif
    S_Add43(A14, A31, S43, x, sequential43);
    M43.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, S43, B11, M43, total_steps, steps_left - 1, (start_index + 43 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S43.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 43, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M44 = (1.0 / (x) * A43) * (1.0 / (x) * B12 + x * x * B33)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential44) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> T44(mem_mngr.GetMem(start_index, 44, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T44(B11.m(), B11.n());
#endif
    T_Add44(B12, B33, T44, x, sequential44);
    M44.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, A43, T44, M44, total_steps, steps_left - 1, (start_index + 44 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    T44.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 44, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M45 = (1.0 / (x) * A13) * (1.0 / (x) * B31 + x * x * B34)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential45) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> T45(mem_mngr.GetMem(start_index, 45, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T45(B11.m(), B11.n());
#endif
    T_Add45(B31, B34, T45, x, sequential45);
    M45.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, A13, T45, M45, total_steps, steps_left - 1, (start_index + 45 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    T45.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 45, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M46 = (-(1.0 / (x)) * A22) * (1.0 / (x) * B33)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential46) shared(mem_mngr, locker) untied default(shared)
    {
#endif
    M46.UpdateMultiplier(Scalar(-(1.0 / (x))));
    M46.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, A22, B33, M46, total_steps, steps_left - 1, (start_index + 46 - 1) * 46, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(46, total_steps, steps_left, start_index, 46, num_threads)) {
# pragma omp taskwait
    }
#endif

    M_Add1(M13, M26, M43, M45, C11, x, false, beta);
    M_Add2(M9, M16, M27, M28, M31, M44, C12, x, false, beta);
    M_Add3(M1, M8, M11, M17, M24, M39, M46, C13, x, false, beta);
    M_Add4(M8, M14, M20, M30, M34, M45, C14, x, false, beta);
    M_Add5(M1, M8, M11, M14, M19, M20, M23, M33, M39, C21, x, false, beta);
    M_Add6(M2, M7, M14, M18, M22, M38, C22, x, false, beta);
    M_Add7(M6, M17, M33, M42, M46, C23, x, false, beta);
    M_Add8(M6, M14, M15, M36, M46, C24, x, false, beta);
    M_Add9(M8, M9, M11, M16, M32, M39, M43, C31, x, false, beta);
    M_Add10(M3, M9, M12, M29, M35, C32, x, false, beta);
    M_Add11(M19, M36, M39, M41, M42, C33, x, false, beta);
    M_Add12(M2, M4, M5, M9, M14, M21, M25, M39, C34, x, false, beta);
    M_Add13(M13, M22, M29, M30, M31, M40, C41, x, false, beta);
    M_Add14(M2, M3, M4, M9, M14, M22, M28, C42, x, false, beta);
    M_Add15(M3, M7, M21, M37, M39, M44, C43, x, false, beta);
    M_Add16(M5, M10, M12, M18, M28, C44, x, false, beta);

    // Handle edge cases with dynamic peeling
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
    if (total_steps == steps_left) {
        mkl_set_num_threads_local(num_threads);
        mkl_set_dynamic(0);
    }
#endif
    DynamicPeeling(A, B, C, 4, 4, 4, beta);
}

// C := alpha * A * B + beta * C
template <typename Scalar>
double FastMatmul(Matrix<Scalar>& A, Matrix<Scalar>& B, Matrix<Scalar>& C,
    int num_steps, double x=1e-8, int num_threads=12, Scalar alpha=Scalar(1.0), Scalar beta=Scalar(0.0)) {
    MemoryManager<Scalar> mem_mngr;
#ifdef _PARALLEL_
    mem_mngr.Allocate(4, 4, 4, 46, num_steps, A.m(), A.n(), B.n());
#endif
    A.set_multiplier(alpha);
    int num_multiplies_per_step = 46;
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

}  // namespace smirnov444_46_352_approx

#endif  // _smirnov444_46_352_approx_HPP_
