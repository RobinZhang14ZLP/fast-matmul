#ifndef _smirnov555_90_710_approx_HPP_
#define _smirnov555_90_710_approx_HPP_

// This is an automatically generated file from gen.py.
#include "common.hpp"

namespace smirnov555_90_710_approx {

template <typename Scalar>
void S_Add1(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(4.0 * (x * x)) * dataS1[i + j * strideS1] + Scalar(-(1.0 / (x))) * dataS2[i + j * strideS2] + Scalar(-(x * x)) * dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add2(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& S5, Matrix<Scalar>& S6, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataS1[i + j * strideS1] + Scalar(-(x * x)) * dataS2[i + j * strideS2] + Scalar(-(0.2 * (x))) * dataS3[i + j * strideS3] + Scalar(x * x) * dataS4[i + j * strideS4] + Scalar(x * x) * dataS5[i + j * strideS5] + Scalar(-(1.0 / (x))) * dataS6[i + j * strideS6];
        }
    }
}

template <typename Scalar>
void S_Add3(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& S5, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataS1[i + j * strideS1] + Scalar(-(x * x)) * dataS2[i + j * strideS2] + Scalar(-(0.2 * (x))) * dataS3[i + j * strideS3] + Scalar(x * x) * dataS4[i + j * strideS4] + Scalar(x * x) * dataS5[i + j * strideS5];
        }
    }
}

template <typename Scalar>
void S_Add4(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& S5, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(x * x)) * dataS1[i + j * strideS1] + Scalar(x * x) * dataS2[i + j * strideS2] + Scalar(x * x) * dataS3[i + j * strideS3] + Scalar(-(x * x)) * dataS4[i + j * strideS4] + Scalar(1.0 / (x)) * dataS5[i + j * strideS5];
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
            dataC[i + j * strideC] = Scalar(-0.2) * dataS1[i + j * strideS1] + Scalar(-(1.0 / (x))) * dataS2[i + j * strideS2] + Scalar(-5.0) * dataS3[i + j * strideS3] + Scalar(x * x) * dataS4[i + j * strideS4];
        }
    }
}

template <typename Scalar>
void S_Add6(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(0.2) * dataS1[i + j * strideS1] + Scalar(1.0 / (x)) * dataS2[i + j * strideS2] + Scalar(0.8 * (x)) * dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add7(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-5.0) * dataS1[i + j * strideS1] + Scalar(x * x) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add8(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& S5, Matrix<Scalar>& S6, Matrix<Scalar>& S7, Matrix<Scalar>& S8, Matrix<Scalar>& S9, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideS4 = S4.stride();
    const int strideS5 = S5.stride();
    const int strideS6 = S6.stride();
    const int strideS7 = S7.stride();
    const int strideS8 = S8.stride();
    const int strideS9 = S9.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    const Scalar *dataS4 = S4.data();
    const Scalar *dataS5 = S5.data();
    const Scalar *dataS6 = S6.data();
    const Scalar *dataS7 = S7.data();
    const Scalar *dataS8 = S8.data();
    const Scalar *dataS9 = S9.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(0.5) * dataS1[i + j * strideS1] + Scalar(-(0.5 * (x * x))) * dataS2[i + j * strideS2] + Scalar(0.05 * (x)) * dataS3[i + j * strideS3] + Scalar(0.25 * (1.0 / (x))) * dataS4[i + j * strideS4] + Scalar(x * x) * dataS5[i + j * strideS5] + Scalar(-0.125) * dataS6[i + j * strideS6] + Scalar(-0.125) * dataS7[i + j * strideS7] + Scalar(1.0 / (x)) * dataS8[i + j * strideS8] + Scalar(0.025 * (x * x)) * dataS9[i + j * strideS9];
        }
    }
}

template <typename Scalar>
void S_Add9(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataS1[i + j * strideS1] + Scalar(0.25 * (1.0 / (x))) * dataS2[i + j * strideS2] + Scalar(0.25 * (x * x)) * dataS3[i + j * strideS3];
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
            dataC[i + j * strideC] = Scalar(0.5) * dataS1[i + j * strideS1] + Scalar(-0.2) * dataS2[i + j * strideS2] + Scalar(1.0 / (x)) * dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add11(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& S5, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(0.2 * (x))) * dataS1[i + j * strideS1] + Scalar(1.0 / (x)) * dataS2[i + j * strideS2] + Scalar(-(0.2 * (x))) * dataS3[i + j * strideS3] + Scalar(0.2 * (x)) * dataS4[i + j * strideS4] + Scalar(0.25) * dataS5[i + j * strideS5];
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
            dataC[i + j * strideC] = Scalar(0.2) * dataS1[i + j * strideS1] + Scalar(1.0 / (x)) * dataS2[i + j * strideS2];
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
            dataC[i + j * strideC] = Scalar(-(0.2 * (x))) * dataS1[i + j * strideS1] + Scalar(1.0 / (x)) * dataS2[i + j * strideS2] + Scalar(0.2 * (x)) * dataS3[i + j * strideS3] + Scalar(0.25) * dataS4[i + j * strideS4];
        }
    }
}

template <typename Scalar>
void S_Add14(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& S5, Matrix<Scalar>& S6, Matrix<Scalar>& S7, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(x * x)) * dataS1[i + j * strideS1] + Scalar(x * x) * dataS2[i + j * strideS2] + Scalar(-(1.0 / (x))) * dataS3[i + j * strideS3] + Scalar(x * x) * dataS4[i + j * strideS4] + Scalar(1.0 / (x)) * dataS5[i + j * strideS5] + Scalar(-(x * x)) * dataS6[i + j * strideS6] + Scalar(x * x) * dataS7[i + j * strideS7];
        }
    }
}

template <typename Scalar>
void S_Add15(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
void S_Add16(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& S5, Matrix<Scalar>& S6, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(0.2 * (x))) * dataS1[i + j * strideS1] + Scalar(0.6 * (x)) * dataS2[i + j * strideS2] -dataS3[i + j * strideS3] + Scalar(1.0 / (x)) * dataS4[i + j * strideS4] + Scalar(-(0.8 * (x * x))) * dataS5[i + j * strideS5] + Scalar(0.2 * (x)) * dataS6[i + j * strideS6];
        }
    }
}

template <typename Scalar>
void S_Add17(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(0.2 * (x))) * dataS1[i + j * strideS1] + Scalar(1.0 / (x)) * dataS2[i + j * strideS2] + Scalar(0.2 * (x)) * dataS3[i + j * strideS3];
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
            dataC[i + j * strideC] = Scalar(5.0) * dataS1[i + j * strideS1] + Scalar(x * x) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add19(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataS1[i + j * strideS1] + Scalar(x * x) * dataS2[i + j * strideS2] + Scalar(1.0 / (x)) * dataS3[i + j * strideS3] + Scalar(x * x) * dataS4[i + j * strideS4];
        }
    }
}

template <typename Scalar>
void S_Add20(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(0.2 * (x))) * dataS1[i + j * strideS1] + Scalar(1.0 / (x)) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add21(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = -dataS1[i + j * strideS1] + Scalar(x * x) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add23(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
void S_Add24(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
void S_Add25(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x * x) * dataS1[i + j * strideS1] + Scalar(1.0 / (x)) * dataS2[i + j * strideS2] + Scalar(1.0 / (x)) * dataS3[i + j * strideS3] + Scalar(1.0 / (x)) * dataS4[i + j * strideS4];
        }
    }
}

template <typename Scalar>
void S_Add26(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(x * x)) * dataS1[i + j * strideS1] + Scalar(1.0 / (x)) * dataS2[i + j * strideS2] + Scalar(-(x * x)) * dataS3[i + j * strideS3] + Scalar(-(x * x)) * dataS4[i + j * strideS4];
        }
    }
}

template <typename Scalar>
void S_Add27(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
void S_Add28(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(5.0) * dataS1[i + j * strideS1];
        }
    }
}

template <typename Scalar>
void S_Add29(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1] + Scalar(-0.5) * dataS2[i + j * strideS2] + Scalar(0.2) * dataS3[i + j * strideS3] + Scalar(-(1.0 / (x))) * dataS4[i + j * strideS4];
        }
    }
}

template <typename Scalar>
void S_Add30(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& S5, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(0.2 * (x))) * dataS1[i + j * strideS1] + Scalar(-(1.0 / (x))) * dataS2[i + j * strideS2] + Scalar(0.5) * dataS3[i + j * strideS3] + Scalar(0.5) * dataS4[i + j * strideS4] + Scalar(-(0.1 * (x * x))) * dataS5[i + j * strideS5];
        }
    }
}

template <typename Scalar>
void S_Add31(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x * x) * dataS1[i + j * strideS1] + Scalar(-(1.0 / (x))) * dataS2[i + j * strideS2] + Scalar(x * x) * dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add32(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-0.2) * dataS1[i + j * strideS1] + Scalar(1.0 / (x)) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add33(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataS1[i + j * strideS1] + Scalar(1.0 / (x)) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add34(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x * x) * dataS1[i + j * strideS1] + Scalar(-(x * x)) * dataS2[i + j * strideS2] + Scalar(x * x) * dataS3[i + j * strideS3] + Scalar(-(1.0 / (x))) * dataS4[i + j * strideS4];
        }
    }
}

template <typename Scalar>
void S_Add35(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
void S_Add36(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1] + Scalar(-0.2) * dataS2[i + j * strideS2] + Scalar(-(1.0 / (x))) * dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add37(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataS1[i + j * strideS1] + Scalar(-0.2) * dataS2[i + j * strideS2] + Scalar(1.0 / (x)) * dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add38(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
void S_Add39(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1] + Scalar(5.0) * dataS2[i + j * strideS2] + Scalar(1.0 / (x)) * dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add40(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
void S_Add41(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1] + Scalar(x * x) * dataS2[i + j * strideS2] + Scalar(0.2 * (x)) * dataS3[i + j * strideS3] + Scalar(-(x * x)) * dataS4[i + j * strideS4];
        }
    }
}

template <typename Scalar>
void S_Add42(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataS1[i + j * strideS1] + Scalar(1.0 / (x)) * dataS2[i + j * strideS2] + Scalar(x * x) * dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add43(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x, bool sequential) {
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
void S_Add44(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(0.5) * dataS1[i + j * strideS1] + Scalar(-(0.5 * (x * x))) * dataS2[i + j * strideS2] + Scalar(1.0 / (x)) * dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add45(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(x * x)) * dataS1[i + j * strideS1] + Scalar(x * x) * dataS2[i + j * strideS2] + Scalar(1.0 / (x)) * dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add46(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& S5, Matrix<Scalar>& S6, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar((0.2 * (x) + -(x * x))) * dataS1[i + j * strideS1] + Scalar(-0.2) * dataS2[i + j * strideS2] + Scalar(-(1.0 / (x))) * dataS3[i + j * strideS3] + Scalar(-(1.0 / (x))) * dataS4[i + j * strideS4] + Scalar(-(0.8 * (x))) * dataS5[i + j * strideS5] + Scalar(1.0 / (x)) * dataS6[i + j * strideS6];
        }
    }
}

template <typename Scalar>
void S_Add47(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x, bool sequential) {
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
void S_Add48(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1] + Scalar(-(1.0 / (x))) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add49(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x, bool sequential) {
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
void S_Add50(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1] + Scalar(-(x * x)) * dataS2[i + j * strideS2] + Scalar(x * x) * dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add51(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(2.0 * (x)) * dataS1[i + j * strideS1] + Scalar(-5.0) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add52(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1] + Scalar(-(10.0 * (x))) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add53(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
void S_Add54(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = dataS1[i + j * strideS1] + Scalar(-(x * x)) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add55(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x, bool sequential) {
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
void S_Add56(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
void S_Add57(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = dataS1[i + j * strideS1] + Scalar(10.0) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add58(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x * x) * dataS1[i + j * strideS1] + Scalar(-(1.0 / (x))) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add59(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1] + Scalar(x * x) * dataS2[i + j * strideS2] + Scalar(-(x * x)) * dataS3[i + j * strideS3] + Scalar(x * x) * dataS4[i + j * strideS4];
        }
    }
}

template <typename Scalar>
void S_Add60(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
void S_Add61(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
void S_Add62(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = dataS1[i + j * strideS1] -dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add63(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(0.2) * dataS1[i + j * strideS1] + Scalar(1.0 / (x)) * dataS2[i + j * strideS2] + Scalar(-(10.0 * (x))) * dataS3[i + j * strideS3] + Scalar(x * x) * dataS4[i + j * strideS4];
        }
    }
}

template <typename Scalar>
void S_Add64(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
void S_Add65(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& S5, Matrix<Scalar>& S6, Matrix<Scalar>& S7, Matrix<Scalar>& S8, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1] + Scalar(x * x) * dataS2[i + j * strideS2] + Scalar(0.2 * (x)) * dataS3[i + j * strideS3] + Scalar(-(x * x)) * dataS4[i + j * strideS4] + Scalar(1.0 / (x)) * dataS5[i + j * strideS5] + Scalar(-(x * x)) * dataS6[i + j * strideS6] + Scalar(-(1.0 / (x))) * dataS7[i + j * strideS7] + Scalar(x * x) * dataS8[i + j * strideS8];
        }
    }
}

template <typename Scalar>
void S_Add66(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1] + Scalar(-(0.2 * (x))) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add67(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(10.0) * dataS1[i + j * strideS1];
        }
    }
}

template <typename Scalar>
void S_Add68(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x, bool sequential) {
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
void S_Add69(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataS1[i + j * strideS1];
        }
    }
}

template <typename Scalar>
void S_Add70(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& S5, Matrix<Scalar>& S6, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1] + Scalar(-(1.0 / (x))) * dataS2[i + j * strideS2] + Scalar(-(1.0 / (x))) * dataS3[i + j * strideS3] + Scalar(0.25 * (1.0 / (x))) * dataS4[i + j * strideS4] + Scalar(0.25 * (x * x)) * dataS5[i + j * strideS5] + Scalar(x * x) * dataS6[i + j * strideS6];
        }
    }
}

template <typename Scalar>
void S_Add71(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(x * x)) * dataS1[i + j * strideS1] + Scalar(1.0 / (x)) * dataS2[i + j * strideS2] + Scalar(x * x) * dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add72(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x, bool sequential) {
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
void S_Add73(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataS1[i + j * strideS1] + Scalar(x * x) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add74(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x, bool sequential) {
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
void S_Add75(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x * x) * dataS1[i + j * strideS1] + Scalar(-(1.0 / (x))) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add76(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x * x) * dataS1[i + j * strideS1] + Scalar(-(x * x)) * dataS2[i + j * strideS2] + Scalar(10.0) * dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add77(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x, bool sequential) {
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
void S_Add78(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x, bool sequential) {
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
void S_Add79(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1] + Scalar(0.2 * (x)) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add80(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x, bool sequential) {
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
void S_Add81(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
void S_Add82(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x, bool sequential) {
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
void S_Add83(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x, bool sequential) {
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
void S_Add84(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
void S_Add85(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataS1[i + j * strideS1] + Scalar(1.0 / (x)) * dataS2[i + j * strideS2] + Scalar(1.0 / (x)) * dataS3[i + j * strideS3] + Scalar(-(0.2 * (x))) * dataS4[i + j * strideS4];
        }
    }
}

template <typename Scalar>
void S_Add86(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = dataS1[i + j * strideS1] + Scalar(-(x * x)) * dataS2[i + j * strideS2] + Scalar(1.0 / (x)) * dataS3[i + j * strideS3] + Scalar(1.0 / (x)) * dataS4[i + j * strideS4];
        }
    }
}

template <typename Scalar>
void S_Add87(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x * x) * dataS1[i + j * strideS1] + Scalar(1.0 / (x)) * dataS2[i + j * strideS2] + Scalar(-(1.0 / (x))) * dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add88(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x * x) * dataS1[i + j * strideS1] + Scalar(-(1.0 / (x))) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add89(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataS1[i + j * strideS1];
        }
    }
}

template <typename Scalar>
void S_Add90(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideS1 = S1.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataS1[i + j * strideS1];
        }
    }
}

template <typename Scalar>
void T_Add1(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + Scalar(-(1.0 / (x))) * dataT2[i + j * strideT2] + Scalar(x * x) * dataT3[i + j * strideT3];
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
            dataC[i + j * strideC] = Scalar(0.5 * (1.0 / (x))) * dataT1[i + j * strideT1] + Scalar(x * x) * dataT2[i + j * strideT2] + Scalar(0.5 * (x * x)) * dataT3[i + j * strideT3];
        }
    }
}

template <typename Scalar>
void T_Add3(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& T5, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(0.5 * (1.0 / (x))) * dataT1[i + j * strideT1] + Scalar(1.0 / (x)) * dataT2[i + j * strideT2] + Scalar(1.0 / (x)) * dataT3[i + j * strideT3] + Scalar(x * x) * dataT4[i + j * strideT4] + Scalar(0.5 * (x * x)) * dataT5[i + j * strideT5];
        }
    }
}

template <typename Scalar>
void T_Add4(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& T5, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x * x) * dataT1[i + j * strideT1] + Scalar(1.0 / (x)) * dataT2[i + j * strideT2] + Scalar(-(1.0 / (x))) * dataT3[i + j * strideT3] + Scalar(x * x) * dataT4[i + j * strideT4] + Scalar(-(x * x)) * dataT5[i + j * strideT5];
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
            dataC[i + j * strideC] = Scalar(0.25) * dataT1[i + j * strideT1] + Scalar(0.2) * dataT2[i + j * strideT2] + Scalar(0.25 * (x * x)) * dataT3[i + j * strideT3] + Scalar(1.0 / (x)) * dataT4[i + j * strideT4];
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
            dataC[i + j * strideC] = Scalar(5.0) * dataT1[i + j * strideT1] + Scalar(x * x) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add7(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& T5, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(0.25) * dataT1[i + j * strideT1] + Scalar(-0.3) * dataT2[i + j * strideT2] + Scalar(0.25 * (x * x)) * dataT3[i + j * strideT3] + Scalar(1.0 / (x)) * dataT4[i + j * strideT4] + Scalar(0.2 * (x)) * dataT5[i + j * strideT5];
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + Scalar(x * x) * dataT2[i + j * strideT2];
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
            dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataT1[i + j * strideT1] + Scalar(x * x) * dataT2[i + j * strideT2];
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
            dataC[i + j * strideC] = Scalar(2.0 * (x)) * dataT1[i + j * strideT1] + Scalar(1.0 / (x)) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add11(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + Scalar(-5.0) * dataT2[i + j * strideT2] + Scalar(1.0 / (x)) * dataT3[i + j * strideT3];
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
            dataC[i + j * strideC] = -dataT1[i + j * strideT1] + Scalar(0.2) * dataT2[i + j * strideT2] + Scalar(1.0 / (x)) * dataT3[i + j * strideT3];
        }
    }
}

template <typename Scalar>
void T_Add13(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& T5, Matrix<Scalar>& T6, Matrix<Scalar>& T7, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataT1[i + j * strideT1] + Scalar(5.0) * dataT2[i + j * strideT2] + Scalar(-(0.5 * (x * x))) * dataT3[i + j * strideT3] + Scalar(4.0 * (x)) * dataT4[i + j * strideT4] + Scalar(-(1.0 / (x))) * dataT5[i + j * strideT5] + Scalar(x * x) * dataT6[i + j * strideT6] + Scalar(2.4 * (x * x)) * dataT7[i + j * strideT7];
        }
    }
}

template <typename Scalar>
void T_Add14(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(x * x)) * dataT1[i + j * strideT1] + Scalar(1.0 / (x)) * dataT2[i + j * strideT2];
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + Scalar(x * x) * dataT2[i + j * strideT2];
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
            dataC[i + j * strideC] = Scalar(5.0) * dataT1[i + j * strideT1] + Scalar(0.5 * (x * x)) * dataT2[i + j * strideT2] + Scalar(-(4.0 * (x))) * dataT3[i + j * strideT3] + Scalar(1.6 * (x * x)) * dataT4[i + j * strideT4];
        }
    }
}

template <typename Scalar>
void T_Add17(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(5.0) * dataT1[i + j * strideT1] + Scalar(0.5 * (x * x)) * dataT2[i + j * strideT2] + Scalar(-(4.0 * (x))) * dataT3[i + j * strideT3] + Scalar(-(2.4 * (x * x))) * dataT4[i + j * strideT4];
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + Scalar(0.2 * (x)) * dataT2[i + j * strideT2] + Scalar(-(x * x)) * dataT3[i + j * strideT3];
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + Scalar(-(x * x)) * dataT2[i + j * strideT2] + Scalar(0.1 * (x * x * x)) * dataT3[i + j * strideT3];
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
            dataC[i + j * strideC] = Scalar(5.0) * dataT1[i + j * strideT1] + Scalar(x * x) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add21(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& T5, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + Scalar(1.0 / (x)) * dataT2[i + j * strideT2] + Scalar(-(1.0 / (x))) * dataT3[i + j * strideT3] + Scalar(x * x) * dataT4[i + j * strideT4] + Scalar(x * x) * dataT5[i + j * strideT5];
        }
    }
}

template <typename Scalar>
void T_Add22(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(x)) * dataT1[i + j * strideT1] + Scalar(0.5 * (1.0 / (x))) * dataT2[i + j * strideT2] + Scalar(1.0 / (x)) * dataT3[i + j * strideT3];
        }
    }
}

template <typename Scalar>
void T_Add23(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
void T_Add24(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
void T_Add25(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
void T_Add26(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
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
void T_Add27(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
void T_Add28(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& T5, Matrix<Scalar>& T6, Matrix<Scalar>& T7, Matrix<Scalar>& T8, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataT1[i + j * strideT1] + Scalar(0.25) * dataT2[i + j * strideT2] + Scalar(-0.3) * dataT3[i + j * strideT3] + Scalar(-(0.2 * (x))) * dataT4[i + j * strideT4] + Scalar(0.25 * (x * x)) * dataT5[i + j * strideT5] + Scalar(1.0 / (x)) * dataT6[i + j * strideT6] + Scalar(0.2 * (x)) * dataT7[i + j * strideT7] + Scalar(x * x) * dataT8[i + j * strideT8];
        }
    }
}

template <typename Scalar>
void T_Add29(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
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
void T_Add30(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(0.25 * (1.0 / (x))) * dataT1[i + j * strideT1] + Scalar(5.0) * dataT2[i + j * strideT2] + Scalar(0.25 * (x * x)) * dataT3[i + j * strideT3];
        }
    }
}

template <typename Scalar>
void T_Add31(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& T5, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(x * x)) * dataT1[i + j * strideT1] + Scalar(1.0 / (x)) * dataT2[i + j * strideT2] + Scalar(-(1.0 / (x))) * dataT3[i + j * strideT3] + Scalar(-(1.0 / (x))) * dataT4[i + j * strideT4] + Scalar(x * x) * dataT5[i + j * strideT5];
        }
    }
}

template <typename Scalar>
void T_Add32(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(2.0 * (x)) * dataT1[i + j * strideT1] + Scalar(1.0 / (x)) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add33(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + Scalar(1.0 / (x)) * dataT2[i + j * strideT2];
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
void T_Add35(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x * x) * dataT1[i + j * strideT1] + Scalar(1.0 / (x)) * dataT2[i + j * strideT2] + Scalar(-(1.0 / (x))) * dataT3[i + j * strideT3];
        }
    }
}

template <typename Scalar>
void T_Add36(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = -dataT1[i + j * strideT1] + Scalar(0.2) * dataT2[i + j * strideT2] + Scalar(1.0 / (x)) * dataT3[i + j * strideT3] + Scalar(1.0 / (x)) * dataT4[i + j * strideT4];
        }
    }
}

template <typename Scalar>
void T_Add37(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
void T_Add38(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + Scalar(x * x) * dataT2[i + j * strideT2] + Scalar(-(x * x)) * dataT3[i + j * strideT3];
        }
    }
}

template <typename Scalar>
void T_Add39(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1];
        }
    }
}

template <typename Scalar>
void T_Add41(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
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
void T_Add42(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x * x) * dataT1[i + j * strideT1] + Scalar(-(1.0 / (x))) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add43(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& T5, Matrix<Scalar>& T6, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-0.1) * dataT1[i + j * strideT1] -dataT2[i + j * strideT2] + Scalar(0.2) * dataT3[i + j * strideT3] + Scalar(1.0 / (x)) * dataT4[i + j * strideT4] + Scalar(x * x) * dataT5[i + j * strideT5] + Scalar(1.0 / (x)) * dataT6[i + j * strideT6];
        }
    }
}

template <typename Scalar>
void T_Add44(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(5.0) * dataT1[i + j * strideT1];
        }
    }
}

template <typename Scalar>
void T_Add47(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& T5, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + Scalar(5.0) * dataT2[i + j * strideT2] + Scalar(-(1.0 / (x))) * dataT3[i + j * strideT3] + Scalar(x * x) * dataT4[i + j * strideT4] + Scalar(x * x) * dataT5[i + j * strideT5];
        }
    }
}

template <typename Scalar>
void T_Add48(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& T5, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x * x) * dataT1[i + j * strideT1] + Scalar(-(1.0 / (x))) * dataT2[i + j * strideT2] + Scalar(x * x) * dataT3[i + j * strideT3] + Scalar(1.0 / (x)) * dataT4[i + j * strideT4] + Scalar(-(x * x)) * dataT5[i + j * strideT5];
        }
    }
}

template <typename Scalar>
void T_Add49(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& T5, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + Scalar(1.0 / (x)) * dataT2[i + j * strideT2] + Scalar(-(x * x)) * dataT3[i + j * strideT3] + Scalar(x * x) * dataT4[i + j * strideT4] + Scalar(-(x * x)) * dataT5[i + j * strideT5];
        }
    }
}

template <typename Scalar>
void T_Add50(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
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
void T_Add51(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(0.2 * (x))) * dataT1[i + j * strideT1] + Scalar(0.5) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add52(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(0.1) * dataT1[i + j * strideT1] + Scalar(1.0 / (x)) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add53(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
void T_Add54(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + Scalar(-(1.0 / (x))) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add55(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(x * x)) * dataT1[i + j * strideT1] + Scalar(1.0 / (x)) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add56(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x) * dataT1[i + j * strideT1] + Scalar(1.0 / (x)) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add57(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x) * dataT1[i + j * strideT1] + Scalar(1.0 / (x)) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add58(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
void T_Add59(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
void T_Add60(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
void T_Add61(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
void T_Add62(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = dataT1[i + j * strideT1] + dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add63(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
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
void T_Add64(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
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
void T_Add65(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
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
void T_Add66(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + Scalar(5.0) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add67(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& T5, Matrix<Scalar>& T6, Matrix<Scalar>& T7, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(x)) * dataT1[i + j * strideT1] + Scalar(-(1.0 / (x))) * dataT2[i + j * strideT2] + Scalar(0.1 * (x)) * dataT3[i + j * strideT3] + Scalar(-(1.0 / (x))) * dataT4[i + j * strideT4] + Scalar(-(0.1 * (x))) * dataT5[i + j * strideT5] + Scalar(0.1 * (x)) * dataT6[i + j * strideT6] + Scalar(0.25) * dataT7[i + j * strideT7];
        }
    }
}

template <typename Scalar>
void T_Add68(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& T5, Matrix<Scalar>& T6, Matrix<Scalar>& T7, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x) * dataT1[i + j * strideT1] + Scalar(1.0 / (x)) * dataT2[i + j * strideT2] + Scalar(1.0 / (x)) * dataT3[i + j * strideT3] + Scalar(-(1.0 / (x))) * dataT4[i + j * strideT4] + Scalar(1.0 / (x)) * dataT5[i + j * strideT5] + Scalar(x * x) * dataT6[i + j * strideT6] + Scalar(1.0 / (x)) * dataT7[i + j * strideT7];
        }
    }
}

template <typename Scalar>
void T_Add69(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
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
void T_Add70(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
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
void T_Add71(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
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
void T_Add72(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(x * x)) * dataT1[i + j * strideT1] + Scalar(1.0 / (x)) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add73(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + Scalar(-(x * x)) * dataT2[i + j * strideT2] + Scalar(1.0 / (x)) * dataT3[i + j * strideT3];
        }
    }
}

template <typename Scalar>
void T_Add74(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + Scalar(-(1.0 / (x))) * dataT2[i + j * strideT2] + Scalar(1.0 / (x)) * dataT3[i + j * strideT3];
        }
    }
}

template <typename Scalar>
void T_Add75(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + Scalar(1.0 / (x)) * dataT2[i + j * strideT2] + Scalar(5.0) * dataT3[i + j * strideT3];
        }
    }
}

template <typename Scalar>
void T_Add76(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
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
void T_Add77(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x * x) * dataT1[i + j * strideT1] + Scalar(0.5 * (x * x)) * dataT2[i + j * strideT2] + Scalar(-(1.0 / (x))) * dataT3[i + j * strideT3];
        }
    }
}

template <typename Scalar>
void T_Add78(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataT1[i + j * strideT1] + Scalar(-(1.0 / (x))) * dataT2[i + j * strideT2] + Scalar(1.0 / (x)) * dataT3[i + j * strideT3];
        }
    }
}

template <typename Scalar>
void T_Add79(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(x * x) * dataT1[i + j * strideT1] + Scalar(-(1.0 / (x))) * dataT2[i + j * strideT2] + Scalar(5.0) * dataT3[i + j * strideT3];
        }
    }
}

template <typename Scalar>
void T_Add80(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& T5, Matrix<Scalar>& C, double x, bool sequential) {
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
            dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataT1[i + j * strideT1] + Scalar(-(1.0 / (x))) * dataT2[i + j * strideT2] + Scalar(x * x) * dataT3[i + j * strideT3] + Scalar(1.0 / (x)) * dataT4[i + j * strideT4] + Scalar(-(1.0 / (x))) * dataT5[i + j * strideT5];
        }
    }
}

template <typename Scalar>
void T_Add81(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
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
void T_Add82(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
    const int strideT1 = T1.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(5.0) * dataT1[i + j * strideT1];
        }
    }
}

template <typename Scalar>
void T_Add83(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
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
void T_Add84(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
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
void T_Add85(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
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
void T_Add86(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
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
void T_Add87(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
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
void T_Add88(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
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
void T_Add89(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
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
void T_Add90(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x, bool sequential) {
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
void M_Add1(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& M8, Matrix<Scalar>& M9, Matrix<Scalar>& M10, Matrix<Scalar>& M11, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
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
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(0.25 * (x * x)) * dataM1[i + j * strideM1] + Scalar(x * x) * dataM2[i + j * strideM2] + Scalar(0.25 * (x * x)) * dataM3[i + j * strideM3] + Scalar(x * x) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(1.0 / (x)) * dataM6[i + j * strideM6] + Scalar(x * x) * dataM7[i + j * strideM7] + dataM8[i + j * strideM8] + Scalar(1.0 / (x)) * dataM9[i + j * strideM9] + Scalar(1.0 / (x)) * dataM10[i + j * strideM10] + Scalar(1.0 / (x)) * dataM11[i + j * strideM11] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(0.25 * (x * x)) * dataM1[i + j * strideM1] + Scalar(x * x) * dataM2[i + j * strideM2] + Scalar(0.25 * (x * x)) * dataM3[i + j * strideM3] + Scalar(x * x) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(1.0 / (x)) * dataM6[i + j * strideM6] + Scalar(x * x) * dataM7[i + j * strideM7] + dataM8[i + j * strideM8] + Scalar(1.0 / (x)) * dataM9[i + j * strideM9] + Scalar(1.0 / (x)) * dataM10[i + j * strideM10] + Scalar(1.0 / (x)) * dataM11[i + j * strideM11];
            }
        }
    }
}

template <typename Scalar>
void M_Add2(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& M8, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
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
                dataC[i + j * strideC] = Scalar(0.25 * (x * x)) * dataM1[i + j * strideM1] + Scalar(x * x) * dataM2[i + j * strideM2] + Scalar(0.25 * (x * x)) * dataM3[i + j * strideM3] + Scalar(x * x) * dataM4[i + j * strideM4] + Scalar(x) * dataM5[i + j * strideM5] + Scalar(1.0 / (x)) * dataM6[i + j * strideM6] + Scalar(-(x)) * dataM7[i + j * strideM7] + Scalar(1.0 / (x)) * dataM8[i + j * strideM8] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(0.25 * (x * x)) * dataM1[i + j * strideM1] + Scalar(x * x) * dataM2[i + j * strideM2] + Scalar(0.25 * (x * x)) * dataM3[i + j * strideM3] + Scalar(x * x) * dataM4[i + j * strideM4] + Scalar(x) * dataM5[i + j * strideM5] + Scalar(1.0 / (x)) * dataM6[i + j * strideM6] + Scalar(-(x)) * dataM7[i + j * strideM7] + Scalar(1.0 / (x)) * dataM8[i + j * strideM8];
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
                dataC[i + j * strideC] = Scalar(0.25 * (1.0 / (x))) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(-(0.25 * (x * x))) * dataM3[i + j * strideM3] + Scalar(0.25 * (1.0 / (x))) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(x * x) * dataM6[i + j * strideM6] + Scalar(-(0.25 * (x * x))) * dataM7[i + j * strideM7] + Scalar(1.0 / (x)) * dataM8[i + j * strideM8] + Scalar(1.0 / (x)) * dataM9[i + j * strideM9] + Scalar(1.0 / (x)) * dataM10[i + j * strideM10] + Scalar(1.0 / (x)) * dataM11[i + j * strideM11] + Scalar(1.0 / (x)) * dataM12[i + j * strideM12] + Scalar(-(1.0 / (x))) * dataM13[i + j * strideM13] + Scalar(-(1.0 / (x))) * dataM14[i + j * strideM14] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(0.25 * (1.0 / (x))) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(-(0.25 * (x * x))) * dataM3[i + j * strideM3] + Scalar(0.25 * (1.0 / (x))) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(x * x) * dataM6[i + j * strideM6] + Scalar(-(0.25 * (x * x))) * dataM7[i + j * strideM7] + Scalar(1.0 / (x)) * dataM8[i + j * strideM8] + Scalar(1.0 / (x)) * dataM9[i + j * strideM9] + Scalar(1.0 / (x)) * dataM10[i + j * strideM10] + Scalar(1.0 / (x)) * dataM11[i + j * strideM11] + Scalar(1.0 / (x)) * dataM12[i + j * strideM12] + Scalar(-(1.0 / (x))) * dataM13[i + j * strideM13] + Scalar(-(1.0 / (x))) * dataM14[i + j * strideM14];
            }
        }
    }
}

template <typename Scalar>
void M_Add4(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& M8, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
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
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(-(x * x)) * dataM2[i + j * strideM2] + Scalar(-(x * x)) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(-(1.0 / (x))) * dataM6[i + j * strideM6] + Scalar(1.0 / (x)) * dataM7[i + j * strideM7] + Scalar(1.0 / (x)) * dataM8[i + j * strideM8] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(-(x * x)) * dataM2[i + j * strideM2] + Scalar(-(x * x)) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(-(1.0 / (x))) * dataM6[i + j * strideM6] + Scalar(1.0 / (x)) * dataM7[i + j * strideM7] + Scalar(1.0 / (x)) * dataM8[i + j * strideM8];
            }
        }
    }
}

template <typename Scalar>
void M_Add5(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& M8, Matrix<Scalar>& M9, Matrix<Scalar>& M10, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
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
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataM1[i + j * strideM1] + Scalar(0.25 * (x * x)) * dataM2[i + j * strideM2] + Scalar(-(1.0 / (x))) * dataM3[i + j * strideM3] + Scalar(-(2.0 * (x))) * dataM4[i + j * strideM4] + Scalar(-(1.0 / (x))) * dataM5[i + j * strideM5] + Scalar(1.0 / (x)) * dataM6[i + j * strideM6] + Scalar(0.25 * (x * x)) * dataM7[i + j * strideM7] + Scalar(-(1.0 / (x))) * dataM8[i + j * strideM8] + Scalar(2.0 * (x)) * dataM9[i + j * strideM9] + Scalar(-(1.0 / (x))) * dataM10[i + j * strideM10] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataM1[i + j * strideM1] + Scalar(0.25 * (x * x)) * dataM2[i + j * strideM2] + Scalar(-(1.0 / (x))) * dataM3[i + j * strideM3] + Scalar(-(2.0 * (x))) * dataM4[i + j * strideM4] + Scalar(-(1.0 / (x))) * dataM5[i + j * strideM5] + Scalar(1.0 / (x)) * dataM6[i + j * strideM6] + Scalar(0.25 * (x * x)) * dataM7[i + j * strideM7] + Scalar(-(1.0 / (x))) * dataM8[i + j * strideM8] + Scalar(2.0 * (x)) * dataM9[i + j * strideM9] + Scalar(-(1.0 / (x))) * dataM10[i + j * strideM10];
            }
        }
    }
}

template <typename Scalar>
void M_Add6(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
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
                dataC[i + j * strideC] = Scalar(x * x) * dataM1[i + j * strideM1] + Scalar(5.0 * (x * x * x)) * dataM2[i + j * strideM2] + Scalar(x * x) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(x * x) * dataM5[i + j * strideM5] + Scalar(1.0 / (x)) * dataM6[i + j * strideM6] + Scalar(-(1.0 / (x))) * dataM7[i + j * strideM7] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(x * x) * dataM1[i + j * strideM1] + Scalar(5.0 * (x * x * x)) * dataM2[i + j * strideM2] + Scalar(x * x) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(x * x) * dataM5[i + j * strideM5] + Scalar(1.0 / (x)) * dataM6[i + j * strideM6] + Scalar(-(1.0 / (x))) * dataM7[i + j * strideM7];
            }
        }
    }
}

template <typename Scalar>
void M_Add7(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& M8, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
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
                dataC[i + j * strideC] = Scalar(-(2.0 * (x * x))) * dataM1[i + j * strideM1] + Scalar(x * x) * dataM2[i + j * strideM2] + Scalar(x * x) * dataM3[i + j * strideM3] + Scalar(x * x) * dataM4[i + j * strideM4] + Scalar(x * x) * dataM5[i + j * strideM5] + Scalar(x * x) * dataM6[i + j * strideM6] + Scalar(x * x) * dataM7[i + j * strideM7] + Scalar(x * x) * dataM8[i + j * strideM8] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(-(2.0 * (x * x))) * dataM1[i + j * strideM1] + Scalar(x * x) * dataM2[i + j * strideM2] + Scalar(x * x) * dataM3[i + j * strideM3] + Scalar(x * x) * dataM4[i + j * strideM4] + Scalar(x * x) * dataM5[i + j * strideM5] + Scalar(x * x) * dataM6[i + j * strideM6] + Scalar(x * x) * dataM7[i + j * strideM7] + Scalar(x * x) * dataM8[i + j * strideM8];
            }
        }
    }
}

template <typename Scalar>
void M_Add8(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& M8, Matrix<Scalar>& M9, Matrix<Scalar>& M10, Matrix<Scalar>& M11, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
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
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(x * x) * dataM1[i + j * strideM1] + Scalar(-(x * x)) * dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(1.0 / (x)) * dataM6[i + j * strideM6] + Scalar(0.5 * (x * x)) * dataM7[i + j * strideM7] + Scalar(-(x * x)) * dataM8[i + j * strideM8] + Scalar(1.0 / (x)) * dataM9[i + j * strideM9] + Scalar(1.0 / (x)) * dataM10[i + j * strideM10] + Scalar(1.0 / (x)) * dataM11[i + j * strideM11] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(x * x) * dataM1[i + j * strideM1] + Scalar(-(x * x)) * dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(1.0 / (x)) * dataM6[i + j * strideM6] + Scalar(0.5 * (x * x)) * dataM7[i + j * strideM7] + Scalar(-(x * x)) * dataM8[i + j * strideM8] + Scalar(1.0 / (x)) * dataM9[i + j * strideM9] + Scalar(1.0 / (x)) * dataM10[i + j * strideM10] + Scalar(1.0 / (x)) * dataM11[i + j * strideM11];
            }
        }
    }
}

template <typename Scalar>
void M_Add9(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& M8, Matrix<Scalar>& M9, Matrix<Scalar>& M10, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
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
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(0.2 * (x * x)) * dataM1[i + j * strideM1] + Scalar(x * x) * dataM2[i + j * strideM2] + Scalar(-(x * x)) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(x * x) * dataM5[i + j * strideM5] + Scalar(-(x * x)) * dataM6[i + j * strideM6] + Scalar(x * x) * dataM7[i + j * strideM7] + Scalar(-(1.0 / (x))) * dataM8[i + j * strideM8] + Scalar(-(x * x)) * dataM9[i + j * strideM9] + Scalar(1.0 / (x)) * dataM10[i + j * strideM10] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(0.2 * (x * x)) * dataM1[i + j * strideM1] + Scalar(x * x) * dataM2[i + j * strideM2] + Scalar(-(x * x)) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(x * x) * dataM5[i + j * strideM5] + Scalar(-(x * x)) * dataM6[i + j * strideM6] + Scalar(x * x) * dataM7[i + j * strideM7] + Scalar(-(1.0 / (x))) * dataM8[i + j * strideM8] + Scalar(-(x * x)) * dataM9[i + j * strideM9] + Scalar(1.0 / (x)) * dataM10[i + j * strideM10];
            }
        }
    }
}

template <typename Scalar>
void M_Add10(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
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
                dataC[i + j * strideC] = Scalar(0.25) * dataM1[i + j * strideM1] + Scalar(4.0 * (x * x)) * dataM2[i + j * strideM2] + Scalar(0.25 * (1.0 / (x))) * dataM3[i + j * strideM3] + Scalar(-(0.25 * (1.0 / (x)))) * dataM4[i + j * strideM4] + Scalar(-(4.0 * (x * x))) * dataM5[i + j * strideM5] + Scalar(x * x) * dataM6[i + j * strideM6] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(0.25) * dataM1[i + j * strideM1] + Scalar(4.0 * (x * x)) * dataM2[i + j * strideM2] + Scalar(0.25 * (1.0 / (x))) * dataM3[i + j * strideM3] + Scalar(-(0.25 * (1.0 / (x)))) * dataM4[i + j * strideM4] + Scalar(-(4.0 * (x * x))) * dataM5[i + j * strideM5] + Scalar(x * x) * dataM6[i + j * strideM6];
            }
        }
    }
}

template <typename Scalar>
void M_Add11(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& M8, Matrix<Scalar>& M9, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
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
                dataC[i + j * strideC] = Scalar(-0.4) * dataM1[i + j * strideM1] + Scalar(-0.1) * dataM2[i + j * strideM2] + Scalar(0.4) * dataM3[i + j * strideM3] + Scalar(-0.4) * dataM4[i + j * strideM4] + Scalar(x * x) * dataM5[i + j * strideM5] + Scalar(x * x) * dataM6[i + j * strideM6] + Scalar(-(1.0 / (x))) * dataM7[i + j * strideM7] + Scalar(1.0 / (x)) * dataM8[i + j * strideM8] + Scalar(0.1 * (x)) * dataM9[i + j * strideM9] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(-0.4) * dataM1[i + j * strideM1] + Scalar(-0.1) * dataM2[i + j * strideM2] + Scalar(0.4) * dataM3[i + j * strideM3] + Scalar(-0.4) * dataM4[i + j * strideM4] + Scalar(x * x) * dataM5[i + j * strideM5] + Scalar(x * x) * dataM6[i + j * strideM6] + Scalar(-(1.0 / (x))) * dataM7[i + j * strideM7] + Scalar(1.0 / (x)) * dataM8[i + j * strideM8] + Scalar(0.1 * (x)) * dataM9[i + j * strideM9];
            }
        }
    }
}

template <typename Scalar>
void M_Add12(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& M8, Matrix<Scalar>& M9, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
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
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(-(1.0 / (x))) * dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + Scalar(-(1.0 / (x))) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(1.0 / (x)) * dataM6[i + j * strideM6] -dataM7[i + j * strideM7] + Scalar(1.0 / (x)) * dataM8[i + j * strideM8] + Scalar(x * x) * dataM9[i + j * strideM9] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(-(1.0 / (x))) * dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + Scalar(-(1.0 / (x))) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(1.0 / (x)) * dataM6[i + j * strideM6] -dataM7[i + j * strideM7] + Scalar(1.0 / (x)) * dataM8[i + j * strideM8] + Scalar(x * x) * dataM9[i + j * strideM9];
            }
        }
    }
}

template <typename Scalar>
void M_Add13(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& M8, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
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
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(-(1.0 / (x))) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(1.0 / (x)) * dataM6[i + j * strideM6] + Scalar(1.0 / (x)) * dataM7[i + j * strideM7] + Scalar(1.0 / (x)) * dataM8[i + j * strideM8] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(-(1.0 / (x))) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(1.0 / (x)) * dataM6[i + j * strideM6] + Scalar(1.0 / (x)) * dataM7[i + j * strideM7] + Scalar(1.0 / (x)) * dataM8[i + j * strideM8];
            }
        }
    }
}

template <typename Scalar>
void M_Add14(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& M8, Matrix<Scalar>& M9, Matrix<Scalar>& M10, Matrix<Scalar>& M11, Matrix<Scalar>& M12, Matrix<Scalar>& M13, Matrix<Scalar>& M14, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
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
                dataC[i + j * strideC] = Scalar(-(0.2 * (x))) * dataM1[i + j * strideM1] + Scalar(-(0.05 * (x))) * dataM2[i + j * strideM2] + Scalar(0.2 * (x)) * dataM3[i + j * strideM3] + Scalar(-(0.2 * (x))) * dataM4[i + j * strideM4] + Scalar(x * x) * dataM5[i + j * strideM5] + Scalar(x * x) * dataM6[i + j * strideM6] + Scalar(0.2 * (x)) * dataM7[i + j * strideM7] + Scalar(x * x) * dataM8[i + j * strideM8] + Scalar(x * x) * dataM9[i + j * strideM9] + Scalar(-(x * x)) * dataM10[i + j * strideM10] + Scalar(0.1 * (x)) * dataM11[i + j * strideM11] + Scalar((0.1 * (x) + -(x * x))) * dataM12[i + j * strideM12] + Scalar(x * x) * dataM13[i + j * strideM13] + Scalar(x * x) * dataM14[i + j * strideM14] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(-(0.2 * (x))) * dataM1[i + j * strideM1] + Scalar(-(0.05 * (x))) * dataM2[i + j * strideM2] + Scalar(0.2 * (x)) * dataM3[i + j * strideM3] + Scalar(-(0.2 * (x))) * dataM4[i + j * strideM4] + Scalar(x * x) * dataM5[i + j * strideM5] + Scalar(x * x) * dataM6[i + j * strideM6] + Scalar(0.2 * (x)) * dataM7[i + j * strideM7] + Scalar(x * x) * dataM8[i + j * strideM8] + Scalar(x * x) * dataM9[i + j * strideM9] + Scalar(-(x * x)) * dataM10[i + j * strideM10] + Scalar(0.1 * (x)) * dataM11[i + j * strideM11] + Scalar((0.1 * (x) + -(x * x))) * dataM12[i + j * strideM12] + Scalar(x * x) * dataM13[i + j * strideM13] + Scalar(x * x) * dataM14[i + j * strideM14];
            }
        }
    }
}

template <typename Scalar>
void M_Add15(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& M8, Matrix<Scalar>& M9, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
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
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(0.25 * (1.0 / (x))) * dataM2[i + j * strideM2] + Scalar(-(1.0 / (x))) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(1.0 / (x)) * dataM6[i + j * strideM6] + Scalar(-(1.0 / (x))) * dataM7[i + j * strideM7] + Scalar(1.0 / (x)) * dataM8[i + j * strideM8] + Scalar(1.0 / (x)) * dataM9[i + j * strideM9] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(0.25 * (1.0 / (x))) * dataM2[i + j * strideM2] + Scalar(-(1.0 / (x))) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(1.0 / (x)) * dataM6[i + j * strideM6] + Scalar(-(1.0 / (x))) * dataM7[i + j * strideM7] + Scalar(1.0 / (x)) * dataM8[i + j * strideM8] + Scalar(1.0 / (x)) * dataM9[i + j * strideM9];
            }
        }
    }
}

template <typename Scalar>
void M_Add16(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& M8, Matrix<Scalar>& M9, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
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
                dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(-(x * x)) * dataM3[i + j * strideM3] + Scalar(-(0.5 * (1.0 / (x)))) * dataM4[i + j * strideM4] + Scalar(-(x * x)) * dataM5[i + j * strideM5] + Scalar(1.0 / (x)) * dataM6[i + j * strideM6] + Scalar(1.0 / (x)) * dataM7[i + j * strideM7] + Scalar(1.0 / (x)) * dataM8[i + j * strideM8] + Scalar(1.0 / (x)) * dataM9[i + j * strideM9] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(-(x * x)) * dataM3[i + j * strideM3] + Scalar(-(0.5 * (1.0 / (x)))) * dataM4[i + j * strideM4] + Scalar(-(x * x)) * dataM5[i + j * strideM5] + Scalar(1.0 / (x)) * dataM6[i + j * strideM6] + Scalar(1.0 / (x)) * dataM7[i + j * strideM7] + Scalar(1.0 / (x)) * dataM8[i + j * strideM8] + Scalar(1.0 / (x)) * dataM9[i + j * strideM9];
            }
        }
    }
}

template <typename Scalar>
void M_Add17(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
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
                dataC[i + j * strideC] = Scalar(x * x) * dataM1[i + j * strideM1] + Scalar(x * x) * dataM2[i + j * strideM2] + Scalar(x * x) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(1.0 / (x)) * dataM6[i + j * strideM6] + Scalar(x * x) * dataM7[i + j * strideM7] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(x * x) * dataM1[i + j * strideM1] + Scalar(x * x) * dataM2[i + j * strideM2] + Scalar(x * x) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(1.0 / (x)) * dataM6[i + j * strideM6] + Scalar(x * x) * dataM7[i + j * strideM7];
            }
        }
    }
}

template <typename Scalar>
void M_Add18(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
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
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + Scalar(-(1.0 / (x))) * dataM4[i + j * strideM4] + Scalar(-(x * x)) * dataM5[i + j * strideM5] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + Scalar(-(1.0 / (x))) * dataM4[i + j * strideM4] + Scalar(-(x * x)) * dataM5[i + j * strideM5];
            }
        }
    }
}

template <typename Scalar>
void M_Add19(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& M8, Matrix<Scalar>& M9, Matrix<Scalar>& M10, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
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
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(-(0.2 * (x))) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(x * x) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(-(x * x)) * dataM5[i + j * strideM5] + Scalar(x * x) * dataM6[i + j * strideM6] + Scalar(1.0 / (x)) * dataM7[i + j * strideM7] + Scalar(1.0 / (x)) * dataM8[i + j * strideM8] + Scalar(-(1.0 / (x))) * dataM9[i + j * strideM9] + Scalar(1.0 / (x)) * dataM10[i + j * strideM10] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(-(0.2 * (x))) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(x * x) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(-(x * x)) * dataM5[i + j * strideM5] + Scalar(x * x) * dataM6[i + j * strideM6] + Scalar(1.0 / (x)) * dataM7[i + j * strideM7] + Scalar(1.0 / (x)) * dataM8[i + j * strideM8] + Scalar(-(1.0 / (x))) * dataM9[i + j * strideM9] + Scalar(1.0 / (x)) * dataM10[i + j * strideM10];
            }
        }
    }
}

template <typename Scalar>
void M_Add20(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& M8, Matrix<Scalar>& M9, Matrix<Scalar>& M10, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
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
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(0.2 * (x)) * dataM1[i + j * strideM1] + Scalar(x * x) * dataM2[i + j * strideM2] + Scalar(0.2 * (x)) * dataM3[i + j * strideM3] + Scalar(x * x) * dataM4[i + j * strideM4] + Scalar(0.2 * (x)) * dataM5[i + j * strideM5] + Scalar(-(1.0 / (x))) * dataM6[i + j * strideM6] + Scalar(1.0 / (x)) * dataM7[i + j * strideM7] + Scalar(1.0 / (x)) * dataM8[i + j * strideM8] + Scalar(1.0 / (x)) * dataM9[i + j * strideM9] + dataM10[i + j * strideM10] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(0.2 * (x)) * dataM1[i + j * strideM1] + Scalar(x * x) * dataM2[i + j * strideM2] + Scalar(0.2 * (x)) * dataM3[i + j * strideM3] + Scalar(x * x) * dataM4[i + j * strideM4] + Scalar(0.2 * (x)) * dataM5[i + j * strideM5] + Scalar(-(1.0 / (x))) * dataM6[i + j * strideM6] + Scalar(1.0 / (x)) * dataM7[i + j * strideM7] + Scalar(1.0 / (x)) * dataM8[i + j * strideM8] + Scalar(1.0 / (x)) * dataM9[i + j * strideM9] + dataM10[i + j * strideM10];
            }
        }
    }
}

template <typename Scalar>
void M_Add21(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& M8, Matrix<Scalar>& M9, Matrix<Scalar>& M10, Matrix<Scalar>& M11, Matrix<Scalar>& M12, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
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
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(0.25 * (x * x)) * dataM1[i + j * strideM1] + Scalar(x * x) * dataM2[i + j * strideM2] + Scalar(0.25 * (x * x)) * dataM3[i + j * strideM3] + Scalar(x * x) * dataM4[i + j * strideM4] + Scalar(x * x) * dataM5[i + j * strideM5] + Scalar(-(1.0 / (x))) * dataM6[i + j * strideM6] + Scalar(1.0 / (x)) * dataM7[i + j * strideM7] + Scalar(-(1.0 / (x))) * dataM8[i + j * strideM8] + Scalar(-(1.0 / (x))) * dataM9[i + j * strideM9] + Scalar(-(1.0 / (x))) * dataM10[i + j * strideM10] + Scalar(1.0 / (x)) * dataM11[i + j * strideM11] + Scalar(1.0 / (x)) * dataM12[i + j * strideM12] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(0.25 * (x * x)) * dataM1[i + j * strideM1] + Scalar(x * x) * dataM2[i + j * strideM2] + Scalar(0.25 * (x * x)) * dataM3[i + j * strideM3] + Scalar(x * x) * dataM4[i + j * strideM4] + Scalar(x * x) * dataM5[i + j * strideM5] + Scalar(-(1.0 / (x))) * dataM6[i + j * strideM6] + Scalar(1.0 / (x)) * dataM7[i + j * strideM7] + Scalar(-(1.0 / (x))) * dataM8[i + j * strideM8] + Scalar(-(1.0 / (x))) * dataM9[i + j * strideM9] + Scalar(-(1.0 / (x))) * dataM10[i + j * strideM10] + Scalar(1.0 / (x)) * dataM11[i + j * strideM11] + Scalar(1.0 / (x)) * dataM12[i + j * strideM12];
            }
        }
    }
}

template <typename Scalar>
void M_Add22(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& M8, Matrix<Scalar>& M9, Matrix<Scalar>& M10, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
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
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(x * x) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(-(1.0 / (x))) * dataM6[i + j * strideM6] + Scalar(x * x) * dataM7[i + j * strideM7] + Scalar(x * x) * dataM8[i + j * strideM8] + Scalar(1.0 / (x)) * dataM9[i + j * strideM9] + Scalar(-(x * x)) * dataM10[i + j * strideM10] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(x * x) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(-(1.0 / (x))) * dataM6[i + j * strideM6] + Scalar(x * x) * dataM7[i + j * strideM7] + Scalar(x * x) * dataM8[i + j * strideM8] + Scalar(1.0 / (x)) * dataM9[i + j * strideM9] + Scalar(-(x * x)) * dataM10[i + j * strideM10];
            }
        }
    }
}

template <typename Scalar>
void M_Add23(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& M8, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
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
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + Scalar(-(1.0 / (x))) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(-(1.0 / (x))) * dataM6[i + j * strideM6] + Scalar(1.0 / (x)) * dataM7[i + j * strideM7] + Scalar(x * x) * dataM8[i + j * strideM8] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + Scalar(-(1.0 / (x))) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(-(1.0 / (x))) * dataM6[i + j * strideM6] + Scalar(1.0 / (x)) * dataM7[i + j * strideM7] + Scalar(x * x) * dataM8[i + j * strideM8];
            }
        }
    }
}

template <typename Scalar>
void M_Add24(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
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
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + Scalar(x * x) * dataM4[i + j * strideM4] + Scalar(x * x) * dataM5[i + j * strideM5] + Scalar(x * x) * dataM6[i + j * strideM6] + Scalar(x * x) * dataM7[i + j * strideM7] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + Scalar(x * x) * dataM4[i + j * strideM4] + Scalar(x * x) * dataM5[i + j * strideM5] + Scalar(x * x) * dataM6[i + j * strideM6] + Scalar(x * x) * dataM7[i + j * strideM7];
            }
        }
    }
}

template <typename Scalar>
void M_Add25(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& M8, Matrix<Scalar>& M9, Matrix<Scalar>& M10, Matrix<Scalar>& M11, Matrix<Scalar>& C, double x, bool sequential, Scalar beta) {
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
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(0.05 * (x)) * dataM1[i + j * strideM1] + Scalar(0.2 * (x)) * dataM2[i + j * strideM2] + Scalar(x * x) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(0.05 * (x)) * dataM6[i + j * strideM6] + Scalar(1.0 / (x)) * dataM7[i + j * strideM7] + Scalar(-(1.0 / (x))) * dataM8[i + j * strideM8] + Scalar(0.2 * (x)) * dataM9[i + j * strideM9] + Scalar(0.2 * (x)) * dataM10[i + j * strideM10] + dataM11[i + j * strideM11] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for if(!sequential)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(0.05 * (x)) * dataM1[i + j * strideM1] + Scalar(0.2 * (x)) * dataM2[i + j * strideM2] + Scalar(x * x) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + Scalar(1.0 / (x)) * dataM5[i + j * strideM5] + Scalar(0.05 * (x)) * dataM6[i + j * strideM6] + Scalar(1.0 / (x)) * dataM7[i + j * strideM7] + Scalar(-(1.0 / (x))) * dataM8[i + j * strideM8] + Scalar(0.2 * (x)) * dataM9[i + j * strideM9] + Scalar(0.2 * (x)) * dataM10[i + j * strideM10] + dataM11[i + j * strideM11];
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

    Matrix<Scalar> A11 = A.Subblock(5, 5, 1, 1);
    Matrix<Scalar> A12 = A.Subblock(5, 5, 1, 2);
    Matrix<Scalar> A13 = A.Subblock(5, 5, 1, 3);
    Matrix<Scalar> A14 = A.Subblock(5, 5, 1, 4);
    Matrix<Scalar> A15 = A.Subblock(5, 5, 1, 5);
    Matrix<Scalar> A21 = A.Subblock(5, 5, 2, 1);
    Matrix<Scalar> A22 = A.Subblock(5, 5, 2, 2);
    Matrix<Scalar> A23 = A.Subblock(5, 5, 2, 3);
    Matrix<Scalar> A24 = A.Subblock(5, 5, 2, 4);
    Matrix<Scalar> A25 = A.Subblock(5, 5, 2, 5);
    Matrix<Scalar> A31 = A.Subblock(5, 5, 3, 1);
    Matrix<Scalar> A32 = A.Subblock(5, 5, 3, 2);
    Matrix<Scalar> A33 = A.Subblock(5, 5, 3, 3);
    Matrix<Scalar> A34 = A.Subblock(5, 5, 3, 4);
    Matrix<Scalar> A35 = A.Subblock(5, 5, 3, 5);
    Matrix<Scalar> A41 = A.Subblock(5, 5, 4, 1);
    Matrix<Scalar> A42 = A.Subblock(5, 5, 4, 2);
    Matrix<Scalar> A43 = A.Subblock(5, 5, 4, 3);
    Matrix<Scalar> A44 = A.Subblock(5, 5, 4, 4);
    Matrix<Scalar> A45 = A.Subblock(5, 5, 4, 5);
    Matrix<Scalar> A51 = A.Subblock(5, 5, 5, 1);
    Matrix<Scalar> A52 = A.Subblock(5, 5, 5, 2);
    Matrix<Scalar> A53 = A.Subblock(5, 5, 5, 3);
    Matrix<Scalar> A54 = A.Subblock(5, 5, 5, 4);
    Matrix<Scalar> A55 = A.Subblock(5, 5, 5, 5);
    Matrix<Scalar> B11 = B.Subblock(5, 5, 1, 1);
    Matrix<Scalar> B12 = B.Subblock(5, 5, 1, 2);
    Matrix<Scalar> B13 = B.Subblock(5, 5, 1, 3);
    Matrix<Scalar> B14 = B.Subblock(5, 5, 1, 4);
    Matrix<Scalar> B15 = B.Subblock(5, 5, 1, 5);
    Matrix<Scalar> B21 = B.Subblock(5, 5, 2, 1);
    Matrix<Scalar> B22 = B.Subblock(5, 5, 2, 2);
    Matrix<Scalar> B23 = B.Subblock(5, 5, 2, 3);
    Matrix<Scalar> B24 = B.Subblock(5, 5, 2, 4);
    Matrix<Scalar> B25 = B.Subblock(5, 5, 2, 5);
    Matrix<Scalar> B31 = B.Subblock(5, 5, 3, 1);
    Matrix<Scalar> B32 = B.Subblock(5, 5, 3, 2);
    Matrix<Scalar> B33 = B.Subblock(5, 5, 3, 3);
    Matrix<Scalar> B34 = B.Subblock(5, 5, 3, 4);
    Matrix<Scalar> B35 = B.Subblock(5, 5, 3, 5);
    Matrix<Scalar> B41 = B.Subblock(5, 5, 4, 1);
    Matrix<Scalar> B42 = B.Subblock(5, 5, 4, 2);
    Matrix<Scalar> B43 = B.Subblock(5, 5, 4, 3);
    Matrix<Scalar> B44 = B.Subblock(5, 5, 4, 4);
    Matrix<Scalar> B45 = B.Subblock(5, 5, 4, 5);
    Matrix<Scalar> B51 = B.Subblock(5, 5, 5, 1);
    Matrix<Scalar> B52 = B.Subblock(5, 5, 5, 2);
    Matrix<Scalar> B53 = B.Subblock(5, 5, 5, 3);
    Matrix<Scalar> B54 = B.Subblock(5, 5, 5, 4);
    Matrix<Scalar> B55 = B.Subblock(5, 5, 5, 5);
    Matrix<Scalar> C11 = C.Subblock(5, 5, 1, 1);
    Matrix<Scalar> C12 = C.Subblock(5, 5, 1, 2);
    Matrix<Scalar> C13 = C.Subblock(5, 5, 1, 3);
    Matrix<Scalar> C14 = C.Subblock(5, 5, 1, 4);
    Matrix<Scalar> C15 = C.Subblock(5, 5, 1, 5);
    Matrix<Scalar> C21 = C.Subblock(5, 5, 2, 1);
    Matrix<Scalar> C22 = C.Subblock(5, 5, 2, 2);
    Matrix<Scalar> C23 = C.Subblock(5, 5, 2, 3);
    Matrix<Scalar> C24 = C.Subblock(5, 5, 2, 4);
    Matrix<Scalar> C25 = C.Subblock(5, 5, 2, 5);
    Matrix<Scalar> C31 = C.Subblock(5, 5, 3, 1);
    Matrix<Scalar> C32 = C.Subblock(5, 5, 3, 2);
    Matrix<Scalar> C33 = C.Subblock(5, 5, 3, 3);
    Matrix<Scalar> C34 = C.Subblock(5, 5, 3, 4);
    Matrix<Scalar> C35 = C.Subblock(5, 5, 3, 5);
    Matrix<Scalar> C41 = C.Subblock(5, 5, 4, 1);
    Matrix<Scalar> C42 = C.Subblock(5, 5, 4, 2);
    Matrix<Scalar> C43 = C.Subblock(5, 5, 4, 3);
    Matrix<Scalar> C44 = C.Subblock(5, 5, 4, 4);
    Matrix<Scalar> C45 = C.Subblock(5, 5, 4, 5);
    Matrix<Scalar> C51 = C.Subblock(5, 5, 5, 1);
    Matrix<Scalar> C52 = C.Subblock(5, 5, 5, 2);
    Matrix<Scalar> C53 = C.Subblock(5, 5, 5, 3);
    Matrix<Scalar> C54 = C.Subblock(5, 5, 5, 4);
    Matrix<Scalar> C55 = C.Subblock(5, 5, 5, 5);


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
    Matrix<Scalar> M47(mem_mngr.GetMem(start_index, 47, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M48(mem_mngr.GetMem(start_index, 48, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M49(mem_mngr.GetMem(start_index, 49, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M50(mem_mngr.GetMem(start_index, 50, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M51(mem_mngr.GetMem(start_index, 51, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M52(mem_mngr.GetMem(start_index, 52, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M53(mem_mngr.GetMem(start_index, 53, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M54(mem_mngr.GetMem(start_index, 54, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M55(mem_mngr.GetMem(start_index, 55, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M56(mem_mngr.GetMem(start_index, 56, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M57(mem_mngr.GetMem(start_index, 57, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M58(mem_mngr.GetMem(start_index, 58, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M59(mem_mngr.GetMem(start_index, 59, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M60(mem_mngr.GetMem(start_index, 60, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M61(mem_mngr.GetMem(start_index, 61, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M62(mem_mngr.GetMem(start_index, 62, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M63(mem_mngr.GetMem(start_index, 63, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M64(mem_mngr.GetMem(start_index, 64, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M65(mem_mngr.GetMem(start_index, 65, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M66(mem_mngr.GetMem(start_index, 66, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M67(mem_mngr.GetMem(start_index, 67, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M68(mem_mngr.GetMem(start_index, 68, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M69(mem_mngr.GetMem(start_index, 69, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M70(mem_mngr.GetMem(start_index, 70, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M71(mem_mngr.GetMem(start_index, 71, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M72(mem_mngr.GetMem(start_index, 72, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M73(mem_mngr.GetMem(start_index, 73, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M74(mem_mngr.GetMem(start_index, 74, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M75(mem_mngr.GetMem(start_index, 75, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M76(mem_mngr.GetMem(start_index, 76, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M77(mem_mngr.GetMem(start_index, 77, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M78(mem_mngr.GetMem(start_index, 78, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M79(mem_mngr.GetMem(start_index, 79, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M80(mem_mngr.GetMem(start_index, 80, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M81(mem_mngr.GetMem(start_index, 81, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M82(mem_mngr.GetMem(start_index, 82, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M83(mem_mngr.GetMem(start_index, 83, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M84(mem_mngr.GetMem(start_index, 84, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M85(mem_mngr.GetMem(start_index, 85, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M86(mem_mngr.GetMem(start_index, 86, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M87(mem_mngr.GetMem(start_index, 87, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M88(mem_mngr.GetMem(start_index, 88, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M89(mem_mngr.GetMem(start_index, 89, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M90(mem_mngr.GetMem(start_index, 90, total_steps - steps_left, M), C11.m(), C11.m(), C11.n(), C.multiplier());
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
    Matrix<Scalar> M47(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M48(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M49(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M50(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M51(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M52(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M53(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M54(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M55(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M56(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M57(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M58(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M59(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M60(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M61(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M62(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M63(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M64(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M65(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M66(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M67(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M68(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M69(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M70(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M71(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M72(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M73(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M74(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M75(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M76(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M77(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M78(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M79(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M80(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M81(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M82(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M83(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M84(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M85(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M86(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M87(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M88(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M89(C11.m(), C11.n(), C.multiplier());
    Matrix<Scalar> M90(C11.m(), C11.n(), C.multiplier());
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
    bool sequential1 = should_launch_task(90, total_steps, steps_left, start_index, 1, num_threads);
    bool sequential2 = should_launch_task(90, total_steps, steps_left, start_index, 2, num_threads);
    bool sequential3 = should_launch_task(90, total_steps, steps_left, start_index, 3, num_threads);
    bool sequential4 = should_launch_task(90, total_steps, steps_left, start_index, 4, num_threads);
    bool sequential5 = should_launch_task(90, total_steps, steps_left, start_index, 5, num_threads);
    bool sequential6 = should_launch_task(90, total_steps, steps_left, start_index, 6, num_threads);
    bool sequential7 = should_launch_task(90, total_steps, steps_left, start_index, 7, num_threads);
    bool sequential8 = should_launch_task(90, total_steps, steps_left, start_index, 8, num_threads);
    bool sequential9 = should_launch_task(90, total_steps, steps_left, start_index, 9, num_threads);
    bool sequential10 = should_launch_task(90, total_steps, steps_left, start_index, 10, num_threads);
    bool sequential11 = should_launch_task(90, total_steps, steps_left, start_index, 11, num_threads);
    bool sequential12 = should_launch_task(90, total_steps, steps_left, start_index, 12, num_threads);
    bool sequential13 = should_launch_task(90, total_steps, steps_left, start_index, 13, num_threads);
    bool sequential14 = should_launch_task(90, total_steps, steps_left, start_index, 14, num_threads);
    bool sequential15 = should_launch_task(90, total_steps, steps_left, start_index, 15, num_threads);
    bool sequential16 = should_launch_task(90, total_steps, steps_left, start_index, 16, num_threads);
    bool sequential17 = should_launch_task(90, total_steps, steps_left, start_index, 17, num_threads);
    bool sequential18 = should_launch_task(90, total_steps, steps_left, start_index, 18, num_threads);
    bool sequential19 = should_launch_task(90, total_steps, steps_left, start_index, 19, num_threads);
    bool sequential20 = should_launch_task(90, total_steps, steps_left, start_index, 20, num_threads);
    bool sequential21 = should_launch_task(90, total_steps, steps_left, start_index, 21, num_threads);
    bool sequential22 = should_launch_task(90, total_steps, steps_left, start_index, 22, num_threads);
    bool sequential23 = should_launch_task(90, total_steps, steps_left, start_index, 23, num_threads);
    bool sequential24 = should_launch_task(90, total_steps, steps_left, start_index, 24, num_threads);
    bool sequential25 = should_launch_task(90, total_steps, steps_left, start_index, 25, num_threads);
    bool sequential26 = should_launch_task(90, total_steps, steps_left, start_index, 26, num_threads);
    bool sequential27 = should_launch_task(90, total_steps, steps_left, start_index, 27, num_threads);
    bool sequential28 = should_launch_task(90, total_steps, steps_left, start_index, 28, num_threads);
    bool sequential29 = should_launch_task(90, total_steps, steps_left, start_index, 29, num_threads);
    bool sequential30 = should_launch_task(90, total_steps, steps_left, start_index, 30, num_threads);
    bool sequential31 = should_launch_task(90, total_steps, steps_left, start_index, 31, num_threads);
    bool sequential32 = should_launch_task(90, total_steps, steps_left, start_index, 32, num_threads);
    bool sequential33 = should_launch_task(90, total_steps, steps_left, start_index, 33, num_threads);
    bool sequential34 = should_launch_task(90, total_steps, steps_left, start_index, 34, num_threads);
    bool sequential35 = should_launch_task(90, total_steps, steps_left, start_index, 35, num_threads);
    bool sequential36 = should_launch_task(90, total_steps, steps_left, start_index, 36, num_threads);
    bool sequential37 = should_launch_task(90, total_steps, steps_left, start_index, 37, num_threads);
    bool sequential38 = should_launch_task(90, total_steps, steps_left, start_index, 38, num_threads);
    bool sequential39 = should_launch_task(90, total_steps, steps_left, start_index, 39, num_threads);
    bool sequential40 = should_launch_task(90, total_steps, steps_left, start_index, 40, num_threads);
    bool sequential41 = should_launch_task(90, total_steps, steps_left, start_index, 41, num_threads);
    bool sequential42 = should_launch_task(90, total_steps, steps_left, start_index, 42, num_threads);
    bool sequential43 = should_launch_task(90, total_steps, steps_left, start_index, 43, num_threads);
    bool sequential44 = should_launch_task(90, total_steps, steps_left, start_index, 44, num_threads);
    bool sequential45 = should_launch_task(90, total_steps, steps_left, start_index, 45, num_threads);
    bool sequential46 = should_launch_task(90, total_steps, steps_left, start_index, 46, num_threads);
    bool sequential47 = should_launch_task(90, total_steps, steps_left, start_index, 47, num_threads);
    bool sequential48 = should_launch_task(90, total_steps, steps_left, start_index, 48, num_threads);
    bool sequential49 = should_launch_task(90, total_steps, steps_left, start_index, 49, num_threads);
    bool sequential50 = should_launch_task(90, total_steps, steps_left, start_index, 50, num_threads);
    bool sequential51 = should_launch_task(90, total_steps, steps_left, start_index, 51, num_threads);
    bool sequential52 = should_launch_task(90, total_steps, steps_left, start_index, 52, num_threads);
    bool sequential53 = should_launch_task(90, total_steps, steps_left, start_index, 53, num_threads);
    bool sequential54 = should_launch_task(90, total_steps, steps_left, start_index, 54, num_threads);
    bool sequential55 = should_launch_task(90, total_steps, steps_left, start_index, 55, num_threads);
    bool sequential56 = should_launch_task(90, total_steps, steps_left, start_index, 56, num_threads);
    bool sequential57 = should_launch_task(90, total_steps, steps_left, start_index, 57, num_threads);
    bool sequential58 = should_launch_task(90, total_steps, steps_left, start_index, 58, num_threads);
    bool sequential59 = should_launch_task(90, total_steps, steps_left, start_index, 59, num_threads);
    bool sequential60 = should_launch_task(90, total_steps, steps_left, start_index, 60, num_threads);
    bool sequential61 = should_launch_task(90, total_steps, steps_left, start_index, 61, num_threads);
    bool sequential62 = should_launch_task(90, total_steps, steps_left, start_index, 62, num_threads);
    bool sequential63 = should_launch_task(90, total_steps, steps_left, start_index, 63, num_threads);
    bool sequential64 = should_launch_task(90, total_steps, steps_left, start_index, 64, num_threads);
    bool sequential65 = should_launch_task(90, total_steps, steps_left, start_index, 65, num_threads);
    bool sequential66 = should_launch_task(90, total_steps, steps_left, start_index, 66, num_threads);
    bool sequential67 = should_launch_task(90, total_steps, steps_left, start_index, 67, num_threads);
    bool sequential68 = should_launch_task(90, total_steps, steps_left, start_index, 68, num_threads);
    bool sequential69 = should_launch_task(90, total_steps, steps_left, start_index, 69, num_threads);
    bool sequential70 = should_launch_task(90, total_steps, steps_left, start_index, 70, num_threads);
    bool sequential71 = should_launch_task(90, total_steps, steps_left, start_index, 71, num_threads);
    bool sequential72 = should_launch_task(90, total_steps, steps_left, start_index, 72, num_threads);
    bool sequential73 = should_launch_task(90, total_steps, steps_left, start_index, 73, num_threads);
    bool sequential74 = should_launch_task(90, total_steps, steps_left, start_index, 74, num_threads);
    bool sequential75 = should_launch_task(90, total_steps, steps_left, start_index, 75, num_threads);
    bool sequential76 = should_launch_task(90, total_steps, steps_left, start_index, 76, num_threads);
    bool sequential77 = should_launch_task(90, total_steps, steps_left, start_index, 77, num_threads);
    bool sequential78 = should_launch_task(90, total_steps, steps_left, start_index, 78, num_threads);
    bool sequential79 = should_launch_task(90, total_steps, steps_left, start_index, 79, num_threads);
    bool sequential80 = should_launch_task(90, total_steps, steps_left, start_index, 80, num_threads);
    bool sequential81 = should_launch_task(90, total_steps, steps_left, start_index, 81, num_threads);
    bool sequential82 = should_launch_task(90, total_steps, steps_left, start_index, 82, num_threads);
    bool sequential83 = should_launch_task(90, total_steps, steps_left, start_index, 83, num_threads);
    bool sequential84 = should_launch_task(90, total_steps, steps_left, start_index, 84, num_threads);
    bool sequential85 = should_launch_task(90, total_steps, steps_left, start_index, 85, num_threads);
    bool sequential86 = should_launch_task(90, total_steps, steps_left, start_index, 86, num_threads);
    bool sequential87 = should_launch_task(90, total_steps, steps_left, start_index, 87, num_threads);
    bool sequential88 = should_launch_task(90, total_steps, steps_left, start_index, 88, num_threads);
    bool sequential89 = should_launch_task(90, total_steps, steps_left, start_index, 89, num_threads);
    bool sequential90 = should_launch_task(90, total_steps, steps_left, start_index, 90, num_threads);
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
    bool sequential47 = false;
    bool sequential48 = false;
    bool sequential49 = false;
    bool sequential50 = false;
    bool sequential51 = false;
    bool sequential52 = false;
    bool sequential53 = false;
    bool sequential54 = false;
    bool sequential55 = false;
    bool sequential56 = false;
    bool sequential57 = false;
    bool sequential58 = false;
    bool sequential59 = false;
    bool sequential60 = false;
    bool sequential61 = false;
    bool sequential62 = false;
    bool sequential63 = false;
    bool sequential64 = false;
    bool sequential65 = false;
    bool sequential66 = false;
    bool sequential67 = false;
    bool sequential68 = false;
    bool sequential69 = false;
    bool sequential70 = false;
    bool sequential71 = false;
    bool sequential72 = false;
    bool sequential73 = false;
    bool sequential74 = false;
    bool sequential75 = false;
    bool sequential76 = false;
    bool sequential77 = false;
    bool sequential78 = false;
    bool sequential79 = false;
    bool sequential80 = false;
    bool sequential81 = false;
    bool sequential82 = false;
    bool sequential83 = false;
    bool sequential84 = false;
    bool sequential85 = false;
    bool sequential86 = false;
    bool sequential87 = false;
    bool sequential88 = false;
    bool sequential89 = false;
    bool sequential90 = false;
#endif



    // M1 = (4.0 * (x * x) * A14 + -(1.0 / (x)) * A45 + -(x * x) * A54) * (1.0 / (x) * B43 + -(1.0 / (x)) * B52 + x * x * B53)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential1) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S1(mem_mngr.GetMem(start_index, 1, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S1(A11.m(), A11.n());
#endif
    S_Add1(A14, A45, A54, S1, x, sequential1);
#ifdef _PARALLEL_
    Matrix<Scalar> T1(mem_mngr.GetMem(start_index, 1, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T1(B11.m(), B11.n());
#endif
    T_Add1(B43, B52, B53, T1, x, sequential1);
    FastMatmulRecursive(locker, mem_mngr, S1, T1, M1, total_steps, steps_left - 1, (start_index + 1 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S1.deallocate();
    T1.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 1, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M2 = (-(1.0 / (x)) * A21 + -(x * x) * A22 + -(0.2 * (x)) * A24 + x * x * A25 + x * x * A42 + -(1.0 / (x)) * A44) * (0.5 * (1.0 / (x)) * B12 + x * x * B41 + 0.5 * (x * x) * B44)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential2) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S2(mem_mngr.GetMem(start_index, 2, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S2(A11.m(), A11.n());
#endif
    S_Add2(A21, A22, A24, A25, A42, A44, S2, x, sequential2);
#ifdef _PARALLEL_
    Matrix<Scalar> T2(mem_mngr.GetMem(start_index, 2, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T2(B11.m(), B11.n());
#endif
    T_Add2(B12, B41, B44, T2, x, sequential2);
    FastMatmulRecursive(locker, mem_mngr, S2, T2, M2, total_steps, steps_left - 1, (start_index + 2 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S2.deallocate();
    T2.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 2, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M3 = (-(1.0 / (x)) * A21 + -(x * x) * A22 + -(0.2 * (x)) * A24 + x * x * A25 + x * x * A42) * (0.5 * (1.0 / (x)) * B12 + 1.0 / (x) * B13 + 1.0 / (x) * B21 + x * x * B41 + 0.5 * (x * x) * B44)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential3) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S3(mem_mngr.GetMem(start_index, 3, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S3(A11.m(), A11.n());
#endif
    S_Add3(A21, A22, A24, A25, A42, S3, x, sequential3);
#ifdef _PARALLEL_
    Matrix<Scalar> T3(mem_mngr.GetMem(start_index, 3, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T3(B11.m(), B11.n());
#endif
    T_Add3(B12, B13, B21, B41, B44, T3, x, sequential3);
    FastMatmulRecursive(locker, mem_mngr, S3, T3, M3, total_steps, steps_left - 1, (start_index + 3 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S3.deallocate();
    T3.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 3, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M4 = (-(x * x) * A13 + x * x * A22 + x * x * A32 + -(x * x) * A42 + 1.0 / (x) * A43) * (x * x * B13 + 1.0 / (x) * B23 + -(1.0 / (x)) * B31 + x * x * B33 + -(x * x) * B43)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential4) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S4(mem_mngr.GetMem(start_index, 4, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S4(A11.m(), A11.n());
#endif
    S_Add4(A13, A22, A32, A42, A43, S4, x, sequential4);
#ifdef _PARALLEL_
    Matrix<Scalar> T4(mem_mngr.GetMem(start_index, 4, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T4(B11.m(), B11.n());
#endif
    T_Add4(B13, B23, B31, B33, B43, T4, x, sequential4);
    FastMatmulRecursive(locker, mem_mngr, S4, T4, M4, total_steps, steps_left - 1, (start_index + 4 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S4.deallocate();
    T4.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 4, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M5 = (-0.2 * A23 + -(1.0 / (x)) * A24 + -5.0 * A34 + x * x * A54) * (0.25 * B35 + 0.2 * B41 + 0.25 * (x * x) * B43 + 1.0 / (x) * B44)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential5) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S5(mem_mngr.GetMem(start_index, 5, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S5(A11.m(), A11.n());
#endif
    S_Add5(A23, A24, A34, A54, S5, x, sequential5);
#ifdef _PARALLEL_
    Matrix<Scalar> T5(mem_mngr.GetMem(start_index, 5, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T5(B11.m(), B11.n());
#endif
    T_Add5(B35, B41, B43, B44, T5, x, sequential5);
    FastMatmulRecursive(locker, mem_mngr, S5, T5, M5, total_steps, steps_left - 1, (start_index + 5 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S5.deallocate();
    T5.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 5, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M6 = (0.2 * A23 + 1.0 / (x) * A24 + 0.8 * (x) * A33) * (5.0 * B35 + x * x * B43)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential6) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S6(mem_mngr.GetMem(start_index, 6, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S6(A11.m(), A11.n());
#endif
    S_Add6(A23, A24, A33, S6, x, sequential6);
#ifdef _PARALLEL_
    Matrix<Scalar> T6(mem_mngr.GetMem(start_index, 6, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T6(B11.m(), B11.n());
#endif
    T_Add6(B35, B43, T6, x, sequential6);
    FastMatmulRecursive(locker, mem_mngr, S6, T6, M6, total_steps, steps_left - 1, (start_index + 6 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S6.deallocate();
    T6.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 6, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M7 = (-5.0 * A34 + x * x * A54) * (0.25 * B35 + -0.3 * B41 + 0.25 * (x * x) * B43 + 1.0 / (x) * B44 + 0.2 * (x) * B45)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential7) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S7(mem_mngr.GetMem(start_index, 7, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S7(A11.m(), A11.n());
#endif
    S_Add7(A34, A54, S7, x, sequential7);
#ifdef _PARALLEL_
    Matrix<Scalar> T7(mem_mngr.GetMem(start_index, 7, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T7(B11.m(), B11.n());
#endif
    T_Add7(B35, B41, B43, B44, B45, T7, x, sequential7);
    FastMatmulRecursive(locker, mem_mngr, S7, T7, M7, total_steps, steps_left - 1, (start_index + 7 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S7.deallocate();
    T7.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 7, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M8 = (0.5 * A11 + -(0.5 * (x * x)) * A12 + 0.05 * (x) * A14 + 0.25 * (1.0 / (x)) * A21 + x * x * A31 + -0.125 * A32 + -0.125 * A34 + 1.0 / (x) * A35 + 0.025 * (x * x) * A54) * (1.0 / (x) * B15 + x * x * B55)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential8) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S8(mem_mngr.GetMem(start_index, 8, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S8(A11.m(), A11.n());
#endif
    S_Add8(A11, A12, A14, A21, A31, A32, A34, A35, A54, S8, x, sequential8);
#ifdef _PARALLEL_
    Matrix<Scalar> T8(mem_mngr.GetMem(start_index, 8, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T8(B11.m(), B11.n());
#endif
    T_Add8(B15, B55, T8, x, sequential8);
    FastMatmulRecursive(locker, mem_mngr, S8, T8, M8, total_steps, steps_left - 1, (start_index + 8 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S8.deallocate();
    T8.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 8, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M9 = (-(1.0 / (x)) * A15 + 0.25 * (1.0 / (x)) * A45 + 0.25 * (x * x) * A54) * (-(1.0 / (x)) * B52 + x * x * B55)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential9) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S9(mem_mngr.GetMem(start_index, 9, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S9(A11.m(), A11.n());
#endif
    S_Add9(A15, A45, A54, S9, x, sequential9);
#ifdef _PARALLEL_
    Matrix<Scalar> T9(mem_mngr.GetMem(start_index, 9, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T9(B11.m(), B11.n());
#endif
    T_Add9(B52, B55, T9, x, sequential9);
    FastMatmulRecursive(locker, mem_mngr, S9, T9, M9, total_steps, steps_left - 1, (start_index + 9 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S9.deallocate();
    T9.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 9, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M10 = (0.5 * A53 + -0.2 * A54 + 1.0 / (x) * A55) * (2.0 * (x) * B33 + 1.0 / (x) * B55)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential10) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S10(mem_mngr.GetMem(start_index, 10, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S10(A11.m(), A11.n());
#endif
    S_Add10(A53, A54, A55, S10, x, sequential10);
#ifdef _PARALLEL_
    Matrix<Scalar> T10(mem_mngr.GetMem(start_index, 10, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T10(B11.m(), B11.n());
#endif
    T_Add10(B33, B55, T10, x, sequential10);
    FastMatmulRecursive(locker, mem_mngr, S10, T10, M10, total_steps, steps_left - 1, (start_index + 10 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S10.deallocate();
    T10.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 10, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M11 = (-(0.2 * (x)) * A13 + 1.0 / (x) * A25 + -(0.2 * (x)) * A52 + 0.2 * (x) * A53 + 0.25 * A54) * (1.0 / (x) * B21 + -5.0 * B25 + 1.0 / (x) * B51)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential11) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S11(mem_mngr.GetMem(start_index, 11, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S11(A11.m(), A11.n());
#endif
    S_Add11(A13, A25, A52, A53, A54, S11, x, sequential11);
#ifdef _PARALLEL_
    Matrix<Scalar> T11(mem_mngr.GetMem(start_index, 11, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T11(B11.m(), B11.n());
#endif
    T_Add11(B21, B25, B51, T11, x, sequential11);
    FastMatmulRecursive(locker, mem_mngr, S11, T11, M11, total_steps, steps_left - 1, (start_index + 11 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S11.deallocate();
    T11.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 11, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M12 = (0.2 * A23 + 1.0 / (x) * A24) * (-1.0 * B35 + 0.2 * B41 + 1.0 / (x) * B44)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential12) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S12(mem_mngr.GetMem(start_index, 12, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S12(A11.m(), A11.n());
#endif
    S_Add12(A23, A24, S12, x, sequential12);
#ifdef _PARALLEL_
    Matrix<Scalar> T12(mem_mngr.GetMem(start_index, 12, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T12(B11.m(), B11.n());
#endif
    T_Add12(B35, B41, B44, T12, x, sequential12);
    FastMatmulRecursive(locker, mem_mngr, S12, T12, M12, total_steps, steps_left - 1, (start_index + 12 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S12.deallocate();
    T12.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 12, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M13 = (-(0.2 * (x)) * A13 + 1.0 / (x) * A25 + 0.2 * (x) * A53 + 0.25 * A54) * (-(1.0 / (x)) * B21 + 5.0 * B25 + -(0.5 * (x * x)) * B43 + 4.0 * (x) * B45 + -(1.0 / (x)) * B51 + x * x * B53 + 2.4 * (x * x) * B55)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential13) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S13(mem_mngr.GetMem(start_index, 13, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S13(A11.m(), A11.n());
#endif
    S_Add13(A13, A25, A53, A54, S13, x, sequential13);
#ifdef _PARALLEL_
    Matrix<Scalar> T13(mem_mngr.GetMem(start_index, 13, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T13(B11.m(), B11.n());
#endif
    T_Add13(B21, B25, B43, B45, B51, B53, B55, T13, x, sequential13);
    FastMatmulRecursive(locker, mem_mngr, S13, T13, M13, total_steps, steps_left - 1, (start_index + 13 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S13.deallocate();
    T13.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 13, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M14 = (-(x * x) * A13 + x * x * A22 + -(1.0 / (x)) * A31 + x * x * A32 + 1.0 / (x) * A33 + -(x * x) * A42 + x * x * A53) * (-(x * x) * B13 + 1.0 / (x) * B31)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential14) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S14(mem_mngr.GetMem(start_index, 14, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S14(A11.m(), A11.n());
#endif
    S_Add14(A13, A22, A31, A32, A33, A42, A53, S14, x, sequential14);
#ifdef _PARALLEL_
    Matrix<Scalar> T14(mem_mngr.GetMem(start_index, 14, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T14(B11.m(), B11.n());
#endif
    T_Add14(B13, B31, T14, x, sequential14);
    FastMatmulRecursive(locker, mem_mngr, S14, T14, M14, total_steps, steps_left - 1, (start_index + 14 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S14.deallocate();
    T14.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 14, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M15 = (1.0 / (x) * A45 + x * x * A54) * (1.0 / (x) * B43 + x * x * B53)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential15) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S15(mem_mngr.GetMem(start_index, 15, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S15(A11.m(), A11.n());
#endif
    S_Add15(A45, A54, S15, x, sequential15);
#ifdef _PARALLEL_
    Matrix<Scalar> T15(mem_mngr.GetMem(start_index, 15, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T15(B11.m(), B11.n());
#endif
    T_Add15(B43, B53, T15, x, sequential15);
    FastMatmulRecursive(locker, mem_mngr, S15, T15, M15, total_steps, steps_left - 1, (start_index + 15 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S15.deallocate();
    T15.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 15, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M16 = (-(0.2 * (x)) * A13 + 0.6 * (x) * A23 + -1.0 * A24 + 1.0 / (x) * A25 + -(0.8 * (x * x)) * A33 + 0.2 * (x) * A53) * (5.0 * B35 + 0.5 * (x * x) * B43 + -(4.0 * (x)) * B45 + 1.6 * (x * x) * B55)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential16) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S16(mem_mngr.GetMem(start_index, 16, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S16(A11.m(), A11.n());
#endif
    S_Add16(A13, A23, A24, A25, A33, A53, S16, x, sequential16);
#ifdef _PARALLEL_
    Matrix<Scalar> T16(mem_mngr.GetMem(start_index, 16, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T16(B11.m(), B11.n());
#endif
    T_Add16(B35, B43, B45, B55, T16, x, sequential16);
    FastMatmulRecursive(locker, mem_mngr, S16, T16, M16, total_steps, steps_left - 1, (start_index + 16 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S16.deallocate();
    T16.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 16, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M17 = (-(0.2 * (x)) * A13 + 1.0 / (x) * A25 + 0.2 * (x) * A53) * (5.0 * B35 + 0.5 * (x * x) * B43 + -(4.0 * (x)) * B45 + -(2.4 * (x * x)) * B55)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential17) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S17(mem_mngr.GetMem(start_index, 17, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S17(A11.m(), A11.n());
#endif
    S_Add17(A13, A25, A53, S17, x, sequential17);
#ifdef _PARALLEL_
    Matrix<Scalar> T17(mem_mngr.GetMem(start_index, 17, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T17(B11.m(), B11.n());
#endif
    T_Add17(B35, B43, B45, B55, T17, x, sequential17);
    FastMatmulRecursive(locker, mem_mngr, S17, T17, M17, total_steps, steps_left - 1, (start_index + 17 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S17.deallocate();
    T17.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 17, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M18 = (5.0 * A34 + x * x * A53) * (1.0 / (x) * B34 + 0.2 * (x) * B42 + -(x * x) * B52)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential18) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S18(mem_mngr.GetMem(start_index, 18, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S18(A11.m(), A11.n());
#endif
    S_Add18(A34, A53, S18, x, sequential18);
#ifdef _PARALLEL_
    Matrix<Scalar> T18(mem_mngr.GetMem(start_index, 18, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T18(B11.m(), B11.n());
#endif
    T_Add18(B34, B42, B52, T18, x, sequential18);
    FastMatmulRecursive(locker, mem_mngr, S18, T18, M18, total_steps, steps_left - 1, (start_index + 18 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S18.deallocate();
    T18.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 18, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M19 = (-(1.0 / (x)) * A23 + x * x * A33 + 1.0 / (x) * A43 + x * x * A53) * (1.0 / (x) * B32 + -(x * x) * B33 + 0.1 * (x * x * x) * B43)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential19) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S19(mem_mngr.GetMem(start_index, 19, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S19(A11.m(), A11.n());
#endif
    S_Add19(A23, A33, A43, A53, S19, x, sequential19);
#ifdef _PARALLEL_
    Matrix<Scalar> T19(mem_mngr.GetMem(start_index, 19, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T19(B11.m(), B11.n());
#endif
    T_Add19(B32, B33, B43, T19, x, sequential19);
    FastMatmulRecursive(locker, mem_mngr, S19, T19, M19, total_steps, steps_left - 1, (start_index + 19 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S19.deallocate();
    T19.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 19, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M20 = (-(0.2 * (x)) * A13 + 1.0 / (x) * A25) * (5.0 * B35 + x * x * B53)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential20) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S20(mem_mngr.GetMem(start_index, 20, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S20(A11.m(), A11.n());
#endif
    S_Add20(A13, A25, S20, x, sequential20);
#ifdef _PARALLEL_
    Matrix<Scalar> T20(mem_mngr.GetMem(start_index, 20, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T20(B11.m(), B11.n());
#endif
    T_Add20(B35, B53, T20, x, sequential20);
    FastMatmulRecursive(locker, mem_mngr, S20, T20, M20, total_steps, steps_left - 1, (start_index + 20 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S20.deallocate();
    T20.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 20, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M21 = (1.0 / (x) * A15) * (1.0 / (x) * B21 + 1.0 / (x) * B51 + -(1.0 / (x)) * B52 + x * x * B53 + x * x * B55)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential21) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> T21(mem_mngr.GetMem(start_index, 21, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T21(B11.m(), B11.n());
#endif
    T_Add21(B21, B51, B52, B53, B55, T21, x, sequential21);
    M21.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, A15, T21, M21, total_steps, steps_left - 1, (start_index + 21 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    T21.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 21, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M22 = (-1.0 * A11 + x * x * A12) * (-(x) * B14 + 0.5 * (1.0 / (x)) * B15 + 1.0 / (x) * B24)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential22) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S22(mem_mngr.GetMem(start_index, 22, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S22(A11.m(), A11.n());
#endif
    S_Add22(A11, A12, S22, x, sequential22);
#ifdef _PARALLEL_
    Matrix<Scalar> T22(mem_mngr.GetMem(start_index, 22, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T22(B11.m(), B11.n());
#endif
    T_Add22(B14, B15, B24, T22, x, sequential22);
    FastMatmulRecursive(locker, mem_mngr, S22, T22, M22, total_steps, steps_left - 1, (start_index + 22 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S22.deallocate();
    T22.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 22, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M23 = (-(x * x) * A12 + 1.0 / (x) * A25) * (1.0 / (x) * B25 + x * x * B54)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential23) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S23(mem_mngr.GetMem(start_index, 23, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S23(A11.m(), A11.n());
#endif
    S_Add23(A12, A25, S23, x, sequential23);
#ifdef _PARALLEL_
    Matrix<Scalar> T23(mem_mngr.GetMem(start_index, 23, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T23(B11.m(), B11.n());
#endif
    T_Add23(B25, B54, T23, x, sequential23);
    FastMatmulRecursive(locker, mem_mngr, S23, T23, M23, total_steps, steps_left - 1, (start_index + 23 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S23.deallocate();
    T23.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 23, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M24 = (1.0 / (x) * A14 + x * x * A55) * (x * x * B45 + 1.0 / (x) * B52)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential24) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S24(mem_mngr.GetMem(start_index, 24, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S24(A11.m(), A11.n());
#endif
    S_Add24(A14, A55, S24, x, sequential24);
#ifdef _PARALLEL_
    Matrix<Scalar> T24(mem_mngr.GetMem(start_index, 24, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T24(B11.m(), B11.n());
#endif
    T_Add24(B45, B52, T24, x, sequential24);
    FastMatmulRecursive(locker, mem_mngr, S24, T24, M24, total_steps, steps_left - 1, (start_index + 24 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S24.deallocate();
    T24.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 24, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M25 = (x * x * A13 + 1.0 / (x) * A33 + 1.0 / (x) * A34 + 1.0 / (x) * A35) * (1.0 / (x) * B34 + -(x * x) * B53)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential25) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S25(mem_mngr.GetMem(start_index, 25, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S25(A11.m(), A11.n());
#endif
    S_Add25(A13, A33, A34, A35, S25, x, sequential25);
#ifdef _PARALLEL_
    Matrix<Scalar> T25(mem_mngr.GetMem(start_index, 25, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T25(B11.m(), B11.n());
#endif
    T_Add25(B34, B53, T25, x, sequential25);
    FastMatmulRecursive(locker, mem_mngr, S25, T25, M25, total_steps, steps_left - 1, (start_index + 25 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S25.deallocate();
    T25.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 25, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M26 = (-(x * x) * A13 + 1.0 / (x) * A23 + -(x * x) * A33 + -(x * x) * A53) * (1.0 / (x) * B32)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential26) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S26(mem_mngr.GetMem(start_index, 26, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S26(A11.m(), A11.n());
#endif
    S_Add26(A13, A23, A33, A53, S26, x, sequential26);
    M26.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, S26, B32, M26, total_steps, steps_left - 1, (start_index + 26 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S26.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 26, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M27 = (x * x * A41 + 1.0 / (x) * A44) * (1.0 / (x) * B12 + x * x * B44)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential27) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S27(mem_mngr.GetMem(start_index, 27, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S27(A11.m(), A11.n());
#endif
    S_Add27(A41, A44, S27, x, sequential27);
#ifdef _PARALLEL_
    Matrix<Scalar> T27(mem_mngr.GetMem(start_index, 27, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T27(B11.m(), B11.n());
#endif
    T_Add27(B12, B44, T27, x, sequential27);
    FastMatmulRecursive(locker, mem_mngr, S27, T27, M27, total_steps, steps_left - 1, (start_index + 27 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S27.deallocate();
    T27.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 27, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M28 = (5.0 * A34) * (-(1.0 / (x)) * B34 + 0.25 * B35 + -0.3 * B41 + -(0.2 * (x)) * B42 + 0.25 * (x * x) * B43 + 1.0 / (x) * B44 + 0.2 * (x) * B45 + x * x * B52)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential28) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> T28(mem_mngr.GetMem(start_index, 28, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T28(B11.m(), B11.n());
#endif
    T_Add28(B34, B35, B41, B42, B43, B44, B45, B52, T28, x, sequential28);
    M28.UpdateMultiplier(Scalar(5.0));
    FastMatmulRecursive(locker, mem_mngr, A34, T28, M28, total_steps, steps_left - 1, (start_index + 28 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    T28.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 28, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M29 = (1.0 / (x) * A45 + -0.5 * A53 + 0.2 * A54 + -(1.0 / (x)) * A55) * (1.0 / (x) * B55)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential29) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S29(mem_mngr.GetMem(start_index, 29, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S29(A11.m(), A11.n());
#endif
    S_Add29(A45, A53, A54, A55, S29, x, sequential29);
    M29.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, S29, B55, M29, total_steps, steps_left - 1, (start_index + 29 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S29.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 29, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M30 = (-(0.2 * (x)) * A14 + -(1.0 / (x)) * A21 + 0.5 * A32 + 0.5 * A34 + -(0.1 * (x * x)) * A54) * (0.25 * (1.0 / (x)) * B15 + 5.0 * B41 + 0.25 * (x * x) * B55)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential30) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S30(mem_mngr.GetMem(start_index, 30, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S30(A11.m(), A11.n());
#endif
    S_Add30(A14, A21, A32, A34, A54, S30, x, sequential30);
#ifdef _PARALLEL_
    Matrix<Scalar> T30(mem_mngr.GetMem(start_index, 30, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T30(B11.m(), B11.n());
#endif
    T_Add30(B15, B41, B55, T30, x, sequential30);
    FastMatmulRecursive(locker, mem_mngr, S30, T30, M30, total_steps, steps_left - 1, (start_index + 30 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S30.deallocate();
    T30.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 30, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M31 = (x * x * A13 + -(1.0 / (x)) * A23 + x * x * A53) * (-(x * x) * B12 + 1.0 / (x) * B31 + -(1.0 / (x)) * B32 + -(1.0 / (x)) * B34 + x * x * B52)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential31) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S31(mem_mngr.GetMem(start_index, 31, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S31(A11.m(), A11.n());
#endif
    S_Add31(A13, A23, A53, S31, x, sequential31);
#ifdef _PARALLEL_
    Matrix<Scalar> T31(mem_mngr.GetMem(start_index, 31, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T31(B11.m(), B11.n());
#endif
    T_Add31(B12, B31, B32, B34, B52, T31, x, sequential31);
    FastMatmulRecursive(locker, mem_mngr, S31, T31, M31, total_steps, steps_left - 1, (start_index + 31 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S31.deallocate();
    T31.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 31, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M32 = (-0.2 * A54 + 1.0 / (x) * A55) * (2.0 * (x) * B33 + 1.0 / (x) * B54)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential32) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S32(mem_mngr.GetMem(start_index, 32, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S32(A11.m(), A11.n());
#endif
    S_Add32(A54, A55, S32, x, sequential32);
#ifdef _PARALLEL_
    Matrix<Scalar> T32(mem_mngr.GetMem(start_index, 32, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T32(B11.m(), B11.n());
#endif
    T_Add32(B33, B54, T32, x, sequential32);
    FastMatmulRecursive(locker, mem_mngr, S32, T32, M32, total_steps, steps_left - 1, (start_index + 32 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S32.deallocate();
    T32.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 32, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M33 = (-(1.0 / (x)) * A15 + 1.0 / (x) * A52) * (1.0 / (x) * B21 + 1.0 / (x) * B51)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential33) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S33(mem_mngr.GetMem(start_index, 33, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S33(A11.m(), A11.n());
#endif
    S_Add33(A15, A52, S33, x, sequential33);
#ifdef _PARALLEL_
    Matrix<Scalar> T33(mem_mngr.GetMem(start_index, 33, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T33(B11.m(), B11.n());
#endif
    T_Add33(B21, B51, T33, x, sequential33);
    FastMatmulRecursive(locker, mem_mngr, S33, T33, M33, total_steps, steps_left - 1, (start_index + 33 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S33.deallocate();
    T33.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 33, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M34 = (x * x * A13 + -(x * x) * A22 + x * x * A42 + -(1.0 / (x)) * A43) * (1.0 / (x) * B23)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential34) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S34(mem_mngr.GetMem(start_index, 34, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S34(A11.m(), A11.n());
#endif
    S_Add34(A13, A22, A42, A43, S34, x, sequential34);
    M34.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, S34, B23, M34, total_steps, steps_left - 1, (start_index + 34 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S34.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 34, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M35 = (-(x * x) * A41 + 1.0 / (x) * A51) * (x * x * B11 + 1.0 / (x) * B12 + -(1.0 / (x)) * B14)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential35) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S35(mem_mngr.GetMem(start_index, 35, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S35(A11.m(), A11.n());
#endif
    S_Add35(A41, A51, S35, x, sequential35);
#ifdef _PARALLEL_
    Matrix<Scalar> T35(mem_mngr.GetMem(start_index, 35, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T35(B11.m(), B11.n());
#endif
    T_Add35(B11, B12, B14, T35, x, sequential35);
    FastMatmulRecursive(locker, mem_mngr, S35, T35, M35, total_steps, steps_left - 1, (start_index + 35 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S35.deallocate();
    T35.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 35, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M36 = (1.0 / (x) * A14 + -0.2 * A23 + -(1.0 / (x)) * A24) * (-1.0 * B35 + 0.2 * B41 + 1.0 / (x) * B42 + 1.0 / (x) * B44)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential36) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S36(mem_mngr.GetMem(start_index, 36, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S36(A11.m(), A11.n());
#endif
    S_Add36(A14, A23, A24, S36, x, sequential36);
#ifdef _PARALLEL_
    Matrix<Scalar> T36(mem_mngr.GetMem(start_index, 36, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T36(B11.m(), B11.n());
#endif
    T_Add36(B35, B41, B42, B44, T36, x, sequential36);
    FastMatmulRecursive(locker, mem_mngr, S36, T36, M36, total_steps, steps_left - 1, (start_index + 36 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S36.deallocate();
    T36.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 36, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M37 = (-(1.0 / (x)) * A45 + -0.2 * A54 + 1.0 / (x) * A55) * (x * x * B53 + 1.0 / (x) * B54)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential37) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S37(mem_mngr.GetMem(start_index, 37, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S37(A11.m(), A11.n());
#endif
    S_Add37(A45, A54, A55, S37, x, sequential37);
#ifdef _PARALLEL_
    Matrix<Scalar> T37(mem_mngr.GetMem(start_index, 37, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T37(B11.m(), B11.n());
#endif
    T_Add37(B53, B54, T37, x, sequential37);
    FastMatmulRecursive(locker, mem_mngr, S37, T37, M37, total_steps, steps_left - 1, (start_index + 37 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S37.deallocate();
    T37.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 37, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M38 = (1.0 / (x) * A34 + 1.0 / (x) * A43) * (1.0 / (x) * B34 + x * x * B43 + -(x * x) * B53)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential38) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S38(mem_mngr.GetMem(start_index, 38, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S38(A11.m(), A11.n());
#endif
    S_Add38(A34, A43, S38, x, sequential38);
#ifdef _PARALLEL_
    Matrix<Scalar> T38(mem_mngr.GetMem(start_index, 38, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T38(B11.m(), B11.n());
#endif
    T_Add38(B34, B43, B53, T38, x, sequential38);
    FastMatmulRecursive(locker, mem_mngr, S38, T38, M38, total_steps, steps_left - 1, (start_index + 38 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S38.deallocate();
    T38.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 38, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M39 = (1.0 / (x) * A23 + 5.0 * A34 + 1.0 / (x) * A35) * (1.0 / (x) * B34 + -(x * x) * B52)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential39) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S39(mem_mngr.GetMem(start_index, 39, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S39(A11.m(), A11.n());
#endif
    S_Add39(A23, A34, A35, S39, x, sequential39);
#ifdef _PARALLEL_
    Matrix<Scalar> T39(mem_mngr.GetMem(start_index, 39, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T39(B11.m(), B11.n());
#endif
    T_Add39(B34, B52, T39, x, sequential39);
    FastMatmulRecursive(locker, mem_mngr, S39, T39, M39, total_steps, steps_left - 1, (start_index + 39 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S39.deallocate();
    T39.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 39, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M40 = (1.0 / (x) * A31 + -(x * x) * A53) * (1.0 / (x) * B31)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential40) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S40(mem_mngr.GetMem(start_index, 40, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S40(A11.m(), A11.n());
#endif
    S_Add40(A31, A53, S40, x, sequential40);
    M40.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, S40, B31, M40, total_steps, steps_left - 1, (start_index + 40 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S40.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 40, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M41 = (1.0 / (x) * A21 + x * x * A22 + 0.2 * (x) * A24 + -(x * x) * A25) * (1.0 / (x) * B21)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential41) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S41(mem_mngr.GetMem(start_index, 41, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S41(A11.m(), A11.n());
#endif
    S_Add41(A21, A22, A24, A25, S41, x, sequential41);
    M41.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, S41, B21, M41, total_steps, steps_left - 1, (start_index + 41 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S41.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 41, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M42 = (-(1.0 / (x)) * A23 + 1.0 / (x) * A31 + x * x * A32) * (x * x * B12 + -(1.0 / (x)) * B31)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential42) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S42(mem_mngr.GetMem(start_index, 42, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S42(A11.m(), A11.n());
#endif
    S_Add42(A23, A31, A32, S42, x, sequential42);
#ifdef _PARALLEL_
    Matrix<Scalar> T42(mem_mngr.GetMem(start_index, 42, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T42(B11.m(), B11.n());
#endif
    T_Add42(B12, B31, T42, x, sequential42);
    FastMatmulRecursive(locker, mem_mngr, S42, T42, M42, total_steps, steps_left - 1, (start_index + 42 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S42.deallocate();
    T42.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 42, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M43 = (-(1.0 / (x)) * A14) * (-0.1 * B32 + -1.0 * B35 + 0.2 * B41 + 1.0 / (x) * B44 + x * x * B45 + 1.0 / (x) * B52)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential43) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> T43(mem_mngr.GetMem(start_index, 43, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T43(B11.m(), B11.n());
#endif
    T_Add43(B32, B35, B41, B44, B45, B52, T43, x, sequential43);
    M43.UpdateMultiplier(Scalar(-(1.0 / (x))));
    FastMatmulRecursive(locker, mem_mngr, A14, T43, M43, total_steps, steps_left - 1, (start_index + 43 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    T43.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 43, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M44 = (0.5 * A11 + -(0.5 * (x * x)) * A12 + 1.0 / (x) * A35) * (1.0 / (x) * B15)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential44) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S44(mem_mngr.GetMem(start_index, 44, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S44(A11.m(), A11.n());
#endif
    S_Add44(A11, A12, A35, S44, x, sequential44);
    M44.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, S44, B15, M44, total_steps, steps_left - 1, (start_index + 44 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S44.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 44, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M45 = (-(x * x) * A13 + x * x * A22 + 1.0 / (x) * A43) * (1.0 / (x) * B23 + x * x * B33)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential45) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S45(mem_mngr.GetMem(start_index, 45, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S45(A11.m(), A11.n());
#endif
    S_Add45(A13, A22, A43, S45, x, sequential45);
#ifdef _PARALLEL_
    Matrix<Scalar> T45(mem_mngr.GetMem(start_index, 45, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T45(B11.m(), B11.n());
#endif
    T_Add45(B23, B33, T45, x, sequential45);
    FastMatmulRecursive(locker, mem_mngr, S45, T45, M45, total_steps, steps_left - 1, (start_index + 45 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S45.deallocate();
    T45.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 45, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M46 = ((0.2 * (x) + -(x * x)) * A13 + -0.2 * A23 + -(1.0 / (x)) * A24 + -(1.0 / (x)) * A25 + -(0.8 * (x)) * A33 + 1.0 / (x) * A43) * (5.0 * B35)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential46) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S46(mem_mngr.GetMem(start_index, 46, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S46(A11.m(), A11.n());
#endif
    S_Add46(A13, A23, A24, A25, A33, A43, S46, x, sequential46);
    M46.UpdateMultiplier(Scalar(5.0));
    FastMatmulRecursive(locker, mem_mngr, S46, B35, M46, total_steps, steps_left - 1, (start_index + 46 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S46.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 46, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M47 = (-(1.0 / (x)) * A25) * (1.0 / (x) * B25 + 5.0 * B35 + -(1.0 / (x)) * B52 + x * x * B53 + x * x * B54)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential47) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> T47(mem_mngr.GetMem(start_index, 47, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T47(B11.m(), B11.n());
#endif
    T_Add47(B25, B35, B52, B53, B54, T47, x, sequential47);
    M47.UpdateMultiplier(Scalar(-(1.0 / (x))));
    FastMatmulRecursive(locker, mem_mngr, A25, T47, M47, total_steps, steps_left - 1, (start_index + 47 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    T47.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 47, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M48 = (1.0 / (x) * A33 + -(1.0 / (x)) * A43) * (x * x * B13 + -(1.0 / (x)) * B31 + x * x * B33 + 1.0 / (x) * B34 + -(x * x) * B53)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential48) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S48(mem_mngr.GetMem(start_index, 48, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S48(A11.m(), A11.n());
#endif
    S_Add48(A33, A43, S48, x, sequential48);
#ifdef _PARALLEL_
    Matrix<Scalar> T48(mem_mngr.GetMem(start_index, 48, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T48(B11.m(), B11.n());
#endif
    T_Add48(B13, B31, B33, B34, B53, T48, x, sequential48);
    FastMatmulRecursive(locker, mem_mngr, S48, T48, M48, total_steps, steps_left - 1, (start_index + 48 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S48.deallocate();
    T48.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 48, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M49 = (-(1.0 / (x)) * A44) * (1.0 / (x) * B12 + 1.0 / (x) * B22 + -(x * x) * B42 + x * x * B44 + -(x * x) * B45)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential49) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> T49(mem_mngr.GetMem(start_index, 49, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T49(B11.m(), B11.n());
#endif
    T_Add49(B12, B22, B42, B44, B45, T49, x, sequential49);
    M49.UpdateMultiplier(Scalar(-(1.0 / (x))));
    FastMatmulRecursive(locker, mem_mngr, A44, T49, M49, total_steps, steps_left - 1, (start_index + 49 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    T49.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 49, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M50 = (1.0 / (x) * A31 + -(x * x) * A52 + x * x * A55) * (1.0 / (x) * B51)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential50) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S50(mem_mngr.GetMem(start_index, 50, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S50(A11.m(), A11.n());
#endif
    S_Add50(A31, A52, A55, S50, x, sequential50);
    M50.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, S50, B51, M50, total_steps, steps_left - 1, (start_index + 50 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S50.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 50, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M51 = (2.0 * (x) * A14 + -5.0 * A32) * (-(0.2 * (x)) * B25 + 0.5 * B41)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential51) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S51(mem_mngr.GetMem(start_index, 51, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S51(A11.m(), A11.n());
#endif
    S_Add51(A14, A32, S51, x, sequential51);
#ifdef _PARALLEL_
    Matrix<Scalar> T51(mem_mngr.GetMem(start_index, 51, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T51(B11.m(), B11.n());
#endif
    T_Add51(B25, B41, T51, x, sequential51);
    FastMatmulRecursive(locker, mem_mngr, S51, T51, M51, total_steps, steps_left - 1, (start_index + 51 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S51.deallocate();
    T51.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 51, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M52 = (1.0 / (x) * A14 + -(10.0 * (x)) * A53) * (0.1 * B32 + 1.0 / (x) * B42)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential52) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S52(mem_mngr.GetMem(start_index, 52, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S52(A11.m(), A11.n());
#endif
    S_Add52(A14, A53, S52, x, sequential52);
#ifdef _PARALLEL_
    Matrix<Scalar> T52(mem_mngr.GetMem(start_index, 52, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T52(B11.m(), B11.n());
#endif
    T_Add52(B32, B42, T52, x, sequential52);
    FastMatmulRecursive(locker, mem_mngr, S52, T52, M52, total_steps, steps_left - 1, (start_index + 52 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S52.deallocate();
    T52.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 52, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M53 = (x * x * A42 + 1.0 / (x) * A44) * (1.0 / (x) * B22 + -(x * x) * B45)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential53) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S53(mem_mngr.GetMem(start_index, 53, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S53(A11.m(), A11.n());
#endif
    S_Add53(A42, A44, S53, x, sequential53);
#ifdef _PARALLEL_
    Matrix<Scalar> T53(mem_mngr.GetMem(start_index, 53, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T53(B11.m(), B11.n());
#endif
    T_Add53(B22, B45, T53, x, sequential53);
    FastMatmulRecursive(locker, mem_mngr, S53, T53, M53, total_steps, steps_left - 1, (start_index + 53 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S53.deallocate();
    T53.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 53, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M54 = (1.0 * A11 + -(x * x) * A42) * (1.0 / (x) * B12 + -(1.0 / (x)) * B24)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential54) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S54(mem_mngr.GetMem(start_index, 54, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S54(A11.m(), A11.n());
#endif
    S_Add54(A11, A42, S54, x, sequential54);
#ifdef _PARALLEL_
    Matrix<Scalar> T54(mem_mngr.GetMem(start_index, 54, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T54(B11.m(), B11.n());
#endif
    T_Add54(B12, B24, T54, x, sequential54);
    FastMatmulRecursive(locker, mem_mngr, S54, T54, M54, total_steps, steps_left - 1, (start_index + 54 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S54.deallocate();
    T54.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 54, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M55 = (1.0 / (x) * A51) * (-(x * x) * B11 + 1.0 / (x) * B14)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential55) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> T55(mem_mngr.GetMem(start_index, 55, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T55(B11.m(), B11.n());
#endif
    T_Add55(B11, B14, T55, x, sequential55);
    M55.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, A51, T55, M55, total_steps, steps_left - 1, (start_index + 55 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    T55.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 55, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M56 = (1.0 * A11 + 1.0 / (x) * A52) * (x * B13 + 1.0 / (x) * B24)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential56) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S56(mem_mngr.GetMem(start_index, 56, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S56(A11.m(), A11.n());
#endif
    S_Add56(A11, A52, S56, x, sequential56);
#ifdef _PARALLEL_
    Matrix<Scalar> T56(mem_mngr.GetMem(start_index, 56, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T56(B11.m(), B11.n());
#endif
    T_Add56(B13, B24, T56, x, sequential56);
    FastMatmulRecursive(locker, mem_mngr, S56, T56, M56, total_steps, steps_left - 1, (start_index + 56 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S56.deallocate();
    T56.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 56, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M57 = (1.0 * A11 + 10.0 * A32) * (x * B11 + 1.0 / (x) * B24)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential57) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S57(mem_mngr.GetMem(start_index, 57, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S57(A11.m(), A11.n());
#endif
    S_Add57(A11, A32, S57, x, sequential57);
#ifdef _PARALLEL_
    Matrix<Scalar> T57(mem_mngr.GetMem(start_index, 57, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T57(B11.m(), B11.n());
#endif
    T_Add57(B11, B24, T57, x, sequential57);
    FastMatmulRecursive(locker, mem_mngr, S57, T57, M57, total_steps, steps_left - 1, (start_index + 57 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S57.deallocate();
    T57.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 57, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M58 = (x * x * A12 + -(1.0 / (x)) * A13) * (1.0 / (x) * B22 + -(x * x) * B32)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential58) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S58(mem_mngr.GetMem(start_index, 58, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S58(A11.m(), A11.n());
#endif
    S_Add58(A12, A13, S58, x, sequential58);
#ifdef _PARALLEL_
    Matrix<Scalar> T58(mem_mngr.GetMem(start_index, 58, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T58(B11.m(), B11.n());
#endif
    T_Add58(B22, B32, T58, x, sequential58);
    FastMatmulRecursive(locker, mem_mngr, S58, T58, M58, total_steps, steps_left - 1, (start_index + 58 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S58.deallocate();
    T58.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 58, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M59 = (1.0 / (x) * A31 + x * x * A35 + -(x * x) * A52 + x * x * A55) * (x * x * B11 + 1.0 / (x) * B51)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential59) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S59(mem_mngr.GetMem(start_index, 59, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S59(A11.m(), A11.n());
#endif
    S_Add59(A31, A35, A52, A55, S59, x, sequential59);
#ifdef _PARALLEL_
    Matrix<Scalar> T59(mem_mngr.GetMem(start_index, 59, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T59(B11.m(), B11.n());
#endif
    T_Add59(B11, B51, T59, x, sequential59);
    FastMatmulRecursive(locker, mem_mngr, S59, T59, M59, total_steps, steps_left - 1, (start_index + 59 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S59.deallocate();
    T59.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 59, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M60 = (1.0 / (x) * A45 + -(x * x) * A51) * (1.0 / (x) * B13 + x * x * B51)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential60) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S60(mem_mngr.GetMem(start_index, 60, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S60(A11.m(), A11.n());
#endif
    S_Add60(A45, A51, S60, x, sequential60);
#ifdef _PARALLEL_
    Matrix<Scalar> T60(mem_mngr.GetMem(start_index, 60, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T60(B11.m(), B11.n());
#endif
    T_Add60(B13, B51, T60, x, sequential60);
    FastMatmulRecursive(locker, mem_mngr, S60, T60, M60, total_steps, steps_left - 1, (start_index + 60 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S60.deallocate();
    T60.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 60, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M61 = (1.0 / (x) * A13 + 1.0 / (x) * A52) * (1.0 / (x) * B22 + x * x * B33)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential61) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S61(mem_mngr.GetMem(start_index, 61, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S61(A11.m(), A11.n());
#endif
    S_Add61(A13, A52, S61, x, sequential61);
#ifdef _PARALLEL_
    Matrix<Scalar> T61(mem_mngr.GetMem(start_index, 61, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T61(B11.m(), B11.n());
#endif
    T_Add61(B22, B33, T61, x, sequential61);
    FastMatmulRecursive(locker, mem_mngr, S61, T61, M61, total_steps, steps_left - 1, (start_index + 61 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S61.deallocate();
    T61.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 61, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M62 = (1.0 * A13 + -1.0 * A32) * (1.0 * B22 + 1.0 * B31)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential62) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S62(mem_mngr.GetMem(start_index, 62, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S62(A11.m(), A11.n());
#endif
    S_Add62(A13, A32, S62, x, sequential62);
#ifdef _PARALLEL_
    Matrix<Scalar> T62(mem_mngr.GetMem(start_index, 62, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T62(B11.m(), B11.n());
#endif
    T_Add62(B22, B31, T62, x, sequential62);
    FastMatmulRecursive(locker, mem_mngr, S62, T62, M62, total_steps, steps_left - 1, (start_index + 62 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S62.deallocate();
    T62.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 62, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M63 = (0.2 * A23 + 1.0 / (x) * A24 + -(10.0 * (x)) * A53 + x * x * A54) * (1.0 / (x) * B42)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential63) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S63(mem_mngr.GetMem(start_index, 63, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S63(A11.m(), A11.n());
#endif
    S_Add63(A23, A24, A53, A54, S63, x, sequential63);
    M63.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, S63, B42, M63, total_steps, steps_left - 1, (start_index + 63 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S63.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 63, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M64 = (x * x * A13 + 1.0 / (x) * A35) * (1.0 / (x) * B34)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential64) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S64(mem_mngr.GetMem(start_index, 64, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S64(A11.m(), A11.n());
#endif
    S_Add64(A13, A35, S64, x, sequential64);
    M64.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, S64, B34, M64, total_steps, steps_left - 1, (start_index + 64 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S64.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 64, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M65 = (1.0 / (x) * A21 + x * x * A22 + 0.2 * (x) * A24 + -(x * x) * A25 + 1.0 / (x) * A41 + -(x * x) * A42 + -(1.0 / (x)) * A45 + x * x * A51) * (1.0 / (x) * B13)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential65) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S65(mem_mngr.GetMem(start_index, 65, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S65(A11.m(), A11.n());
#endif
    S_Add65(A21, A22, A24, A25, A41, A42, A45, A51, S65, x, sequential65);
    M65.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, S65, B13, M65, total_steps, steps_left - 1, (start_index + 65 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S65.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 65, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M66 = (1.0 / (x) * A21 + -(0.2 * (x)) * A54) * (1.0 / (x) * B14 + 5.0 * B41)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential66) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S66(mem_mngr.GetMem(start_index, 66, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S66(A11.m(), A11.n());
#endif
    S_Add66(A21, A54, S66, x, sequential66);
#ifdef _PARALLEL_
    Matrix<Scalar> T66(mem_mngr.GetMem(start_index, 66, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T66(B11.m(), B11.n());
#endif
    T_Add66(B14, B41, T66, x, sequential66);
    FastMatmulRecursive(locker, mem_mngr, S66, T66, M66, total_steps, steps_left - 1, (start_index + 66 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S66.deallocate();
    T66.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 66, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M67 = (10.0 * A32) * (-(x) * B11 + -(1.0 / (x)) * B21 + 0.1 * (x) * B22 + -(1.0 / (x)) * B24 + -(0.1 * (x)) * B25 + 0.1 * (x) * B31 + 0.25 * B41)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential67) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> T67(mem_mngr.GetMem(start_index, 67, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T67(B11.m(), B11.n());
#endif
    T_Add67(B11, B21, B22, B24, B25, B31, B41, T67, x, sequential67);
    M67.UpdateMultiplier(Scalar(10.0));
    FastMatmulRecursive(locker, mem_mngr, A32, T67, M67, total_steps, steps_left - 1, (start_index + 67 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    T67.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 67, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M68 = (-(1.0 / (x)) * A52) * (x * B13 + 1.0 / (x) * B21 + 1.0 / (x) * B22 + -(1.0 / (x)) * B23 + 1.0 / (x) * B24 + x * x * B33 + 1.0 / (x) * B51)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential68) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> T68(mem_mngr.GetMem(start_index, 68, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T68(B11.m(), B11.n());
#endif
    T_Add68(B13, B21, B22, B23, B24, B33, B51, T68, x, sequential68);
    M68.UpdateMultiplier(Scalar(-(1.0 / (x))));
    FastMatmulRecursive(locker, mem_mngr, A52, T68, M68, total_steps, steps_left - 1, (start_index + 68 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    T68.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 68, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M69 = (-1.0 * A11) * (1.0 / (x) * B24)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential69) shared(mem_mngr, locker) untied default(shared)
    {
#endif
    M69.UpdateMultiplier(Scalar(-1.0));
    M69.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, A11, B24, M69, total_steps, steps_left - 1, (start_index + 69 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 69, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M70 = (1.0 / (x) * A14 + -(1.0 / (x)) * A15 + -(1.0 / (x)) * A25 + 0.25 * (1.0 / (x)) * A45 + 0.25 * (x * x) * A54 + x * x * A55) * (1.0 / (x) * B52)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential70) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S70(mem_mngr.GetMem(start_index, 70, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S70(A11.m(), A11.n());
#endif
    S_Add70(A14, A15, A25, A45, A54, A55, S70, x, sequential70);
    M70.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, S70, B52, M70, total_steps, steps_left - 1, (start_index + 70 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S70.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 70, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M71 = (-(x * x) * A44 + 1.0 / (x) * A45 + x * x * A54) * (1.0 / (x) * B43)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential71) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S71(mem_mngr.GetMem(start_index, 71, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S71(A11.m(), A11.n());
#endif
    S_Add71(A44, A45, A54, S71, x, sequential71);
    M71.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, S71, B43, M71, total_steps, steps_left - 1, (start_index + 71 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S71.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 71, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M72 = (1.0 / (x) * A22) * (-(x * x) * B24 + 1.0 / (x) * B25)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential72) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> T72(mem_mngr.GetMem(start_index, 72, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T72(B11.m(), B11.n());
#endif
    T_Add72(B24, B25, T72, x, sequential72);
    M72.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, A22, T72, M72, total_steps, steps_left - 1, (start_index + 72 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    T72.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 72, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M73 = (-(1.0 / (x)) * A22 + x * x * A42) * (1.0 / (x) * B22 + -(x * x) * B24 + 1.0 / (x) * B25)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential73) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S73(mem_mngr.GetMem(start_index, 73, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S73(A11.m(), A11.n());
#endif
    S_Add73(A22, A42, S73, x, sequential73);
#ifdef _PARALLEL_
    Matrix<Scalar> T73(mem_mngr.GetMem(start_index, 73, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T73(B11.m(), B11.n());
#endif
    T_Add73(B22, B24, B25, T73, x, sequential73);
    FastMatmulRecursive(locker, mem_mngr, S73, T73, M73, total_steps, steps_left - 1, (start_index + 73 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S73.deallocate();
    T73.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 73, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M74 = (1.0 / (x) * A31) * (1.0 / (x) * B14 + -(1.0 / (x)) * B31 + 1.0 / (x) * B51)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential74) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> T74(mem_mngr.GetMem(start_index, 74, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T74(B11.m(), B11.n());
#endif
    T_Add74(B14, B31, B51, T74, x, sequential74);
    M74.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, A31, T74, M74, total_steps, steps_left - 1, (start_index + 74 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    T74.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 74, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M75 = (x * x * A13 + -(1.0 / (x)) * A43) * (1.0 / (x) * B23 + 1.0 / (x) * B32 + 5.0 * B35)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential75) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S75(mem_mngr.GetMem(start_index, 75, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S75(A11.m(), A11.n());
#endif
    S_Add75(A13, A43, S75, x, sequential75);
#ifdef _PARALLEL_
    Matrix<Scalar> T75(mem_mngr.GetMem(start_index, 75, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T75(B11.m(), B11.n());
#endif
    T_Add75(B23, B32, B35, T75, x, sequential75);
    FastMatmulRecursive(locker, mem_mngr, S75, T75, M75, total_steps, steps_left - 1, (start_index + 75 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S75.deallocate();
    T75.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 75, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M76 = (x * x * A12 + -(x * x) * A15 + 10.0 * A32) * (1.0 / (x) * B21)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential76) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S76(mem_mngr.GetMem(start_index, 76, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S76(A11.m(), A11.n());
#endif
    S_Add76(A12, A15, A32, S76, x, sequential76);
    M76.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, S76, B21, M76, total_steps, steps_left - 1, (start_index + 76 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S76.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 76, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M77 = (1.0 / (x) * A41) * (x * x * B11 + 0.5 * (x * x) * B12 + -(1.0 / (x)) * B13)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential77) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> T77(mem_mngr.GetMem(start_index, 77, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T77(B11.m(), B11.n());
#endif
    T_Add77(B11, B12, B13, T77, x, sequential77);
    M77.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, A41, T77, M77, total_steps, steps_left - 1, (start_index + 77 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    T77.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 77, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M78 = (1.0 / (x) * A35) * (-(1.0 / (x)) * B15 + -(1.0 / (x)) * B34 + 1.0 / (x) * B54)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential78) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> T78(mem_mngr.GetMem(start_index, 78, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T78(B11.m(), B11.n());
#endif
    T_Add78(B15, B34, B54, T78, x, sequential78);
    M78.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, A35, T78, M78, total_steps, steps_left - 1, (start_index + 78 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    T78.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 78, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M79 = (1.0 / (x) * A21 + 0.2 * (x) * A24) * (x * x * B11 + -(1.0 / (x)) * B21 + 5.0 * B41)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential79) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S79(mem_mngr.GetMem(start_index, 79, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S79(A11.m(), A11.n());
#endif
    S_Add79(A21, A24, S79, x, sequential79);
#ifdef _PARALLEL_
    Matrix<Scalar> T79(mem_mngr.GetMem(start_index, 79, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T79(B11.m(), B11.n());
#endif
    T_Add79(B11, B21, B41, T79, x, sequential79);
    FastMatmulRecursive(locker, mem_mngr, S79, T79, M79, total_steps, steps_left - 1, (start_index + 79 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S79.deallocate();
    T79.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 79, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M80 = (1.0 / (x) * A45) * (1.0 / (x) * B13 + -(1.0 / (x)) * B43 + x * x * B51 + 1.0 / (x) * B54 + -(1.0 / (x)) * B55)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential80) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> T80(mem_mngr.GetMem(start_index, 80, total_steps - steps_left, T), B11.m(), B11.m(), B11.n());
#else
    Matrix<Scalar> T80(B11.m(), B11.n());
#endif
    T_Add80(B13, B43, B51, B54, B55, T80, x, sequential80);
    M80.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, A45, T80, M80, total_steps, steps_left - 1, (start_index + 80 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    T80.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 80, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M81 = (1.0 / (x) * A22 + 1.0 / (x) * A44) * (1.0 / (x) * B22)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential81) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S81(mem_mngr.GetMem(start_index, 81, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S81(A11.m(), A11.n());
#endif
    S_Add81(A22, A44, S81, x, sequential81);
    M81.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, S81, B22, M81, total_steps, steps_left - 1, (start_index + 81 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S81.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 81, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M82 = (1.0 / (x) * A21) * (5.0 * B41)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential82) shared(mem_mngr, locker) untied default(shared)
    {
#endif
    M82.UpdateMultiplier(Scalar(1.0 / (x)));
    M82.UpdateMultiplier(Scalar(5.0));
    FastMatmulRecursive(locker, mem_mngr, A21, B41, M82, total_steps, steps_left - 1, (start_index + 82 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 82, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M83 = (1.0 / (x) * A13) * (1.0 / (x) * B22)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential83) shared(mem_mngr, locker) untied default(shared)
    {
#endif
    M83.UpdateMultiplier(Scalar(1.0 / (x)));
    M83.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, A13, B22, M83, total_steps, steps_left - 1, (start_index + 83 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 83, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M84 = (-(x * x) * A12 + 1.0 / (x) * A52) * (1.0 / (x) * B23)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential84) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S84(mem_mngr.GetMem(start_index, 84, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S84(A11.m(), A11.n());
#endif
    S_Add84(A12, A52, S84, x, sequential84);
    M84.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, S84, B23, M84, total_steps, steps_left - 1, (start_index + 84 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S84.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 84, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M85 = (1.0 / (x) * A21 + 1.0 / (x) * A31 + 1.0 / (x) * A51 + -(0.2 * (x)) * A54) * (1.0 / (x) * B14)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential85) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S85(mem_mngr.GetMem(start_index, 85, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S85(A11.m(), A11.n());
#endif
    S_Add85(A21, A31, A51, A54, S85, x, sequential85);
    M85.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, S85, B14, M85, total_steps, steps_left - 1, (start_index + 85 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S85.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 85, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M86 = (1.0 * A11 + -(x * x) * A42 + 1.0 / (x) * A44 + 1.0 / (x) * A51) * (-(1.0 / (x)) * B12)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential86) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S86(mem_mngr.GetMem(start_index, 86, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S86(A11.m(), A11.n());
#endif
    S_Add86(A11, A42, A44, A51, S86, x, sequential86);
    M86.UpdateMultiplier(Scalar(-(1.0 / (x))));
    FastMatmulRecursive(locker, mem_mngr, S86, B12, M86, total_steps, steps_left - 1, (start_index + 86 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S86.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 86, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M87 = (x * x * A12 + 1.0 / (x) * A22 + -(1.0 / (x)) * A25) * (1.0 / (x) * B25)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential87) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S87(mem_mngr.GetMem(start_index, 87, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S87(A11.m(), A11.n());
#endif
    S_Add87(A12, A22, A25, S87, x, sequential87);
    M87.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, S87, B25, M87, total_steps, steps_left - 1, (start_index + 87 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S87.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 87, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M88 = (x * x * A15 + -(1.0 / (x)) * A35) * (1.0 / (x) * B54)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential88) shared(mem_mngr, locker) untied default(shared)
    {
#endif
#ifdef _PARALLEL_
    Matrix<Scalar> S88(mem_mngr.GetMem(start_index, 88, total_steps - steps_left, S), A11.m(), A11.m(), A11.n());
#else
    Matrix<Scalar> S88(A11.m(), A11.n());
#endif
    S_Add88(A15, A35, S88, x, sequential88);
    M88.UpdateMultiplier(Scalar(1.0 / (x)));
    FastMatmulRecursive(locker, mem_mngr, S88, B54, M88, total_steps, steps_left - 1, (start_index + 88 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
    S88.deallocate();
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 88, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M89 = (1.0 * A41) * (1.0 * B15)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential89) shared(mem_mngr, locker) untied default(shared)
    {
#endif
    M89.UpdateMultiplier(Scalar(1.0));
    M89.UpdateMultiplier(Scalar(1.0));
    FastMatmulRecursive(locker, mem_mngr, A41, B15, M89, total_steps, steps_left - 1, (start_index + 89 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 89, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    SwitchToDFS(locker, num_threads);
# endif
    }
#endif

    // M90 = (1.0 * A51) * (1.0 * B15)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(sequential90) shared(mem_mngr, locker) untied default(shared)
    {
#endif
    M90.UpdateMultiplier(Scalar(1.0));
    M90.UpdateMultiplier(Scalar(1.0));
    FastMatmulRecursive(locker, mem_mngr, A51, B15, M90, total_steps, steps_left - 1, (start_index + 90 - 1) * 90, x, num_threads, Scalar(0.0));
#ifndef _PARALLEL_
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
locker.Decrement();
    }
    if (should_task_wait(90, total_steps, steps_left, start_index, 90, num_threads)) {
# pragma omp taskwait
    }
#endif

    M_Add1(M1, M9, M15, M21, M51, M57, M58, M62, M67, M69, M76, C11, x, false, beta);
    M_Add2(M1, M9, M15, M52, M54, M58, M69, M83, C12, x, false, beta);
    M_Add3(M1, M9, M10, M15, M21, M24, M29, M33, M56, M61, M68, M69, M83, M84, C13, x, false, beta);
    M_Add4(M22, M24, M43, M44, M64, M69, M78, M88, C14, x, false, beta);
    M_Add5(M9, M10, M20, M22, M23, M24, M29, M47, M69, M70, C15, x, false, beta);
    M_Add6(M11, M23, M40, M41, M42, M79, M82, C21, x, false, beta);
    M_Add7(M2, M20, M23, M26, M47, M49, M63, M81, C22, x, false, beta);
    M_Add8(M2, M3, M6, M19, M20, M26, M27, M41, M45, M46, M75, C23, x, false, beta);
    M_Add9(M6, M12, M18, M23, M39, M64, M66, M72, M82, M87, C24, x, false, beta);
    M_Add10(M6, M8, M16, M17, M44, M72, C25, x, false, beta);
    M_Add11(M5, M6, M7, M12, M14, M40, M50, M59, M76, C31, x, false, beta);
    M_Add12(M18, M26, M31, M39, M40, M42, M62, M64, M83, C32, x, false, beta);
    M_Add13(M4, M14, M25, M34, M38, M40, M48, M64, C33, x, false, beta);
    M_Add14(M5, M6, M7, M12, M22, M25, M28, M40, M44, M50, M57, M69, M74, M78, C34, x, false, beta);
    M_Add15(M5, M6, M7, M8, M12, M30, M44, M51, M82, C35, x, false, beta);
    M_Add16(M2, M3, M4, M27, M34, M41, M60, M65, M77, C41, x, false, beta);
    M_Add17(M1, M19, M26, M27, M49, M53, M71, C42, x, false, beta);
    M_Add18(M15, M34, M45, M71, M77, C43, x, false, beta);
    M_Add19(M18, M27, M32, M35, M37, M38, M54, M55, M69, M86, C44, x, false, beta);
    M_Add20(M6, M10, M20, M29, M46, M53, M72, M73, M81, M89, C45, x, false, beta);
    M_Add21(M1, M9, M15, M21, M33, M40, M50, M55, M66, M74, M82, M85, C51, x, false, beta);
    M_Add22(M12, M24, M35, M36, M43, M52, M55, M61, M63, M83, C52, x, false, beta);
    M_Add23(M10, M15, M29, M32, M37, M60, M80, M84, C53, x, false, beta);
    M_Add24(M7, M18, M28, M32, M55, M56, M69, C54, x, false, beta);
    M_Add25(M1, M9, M10, M11, M13, M15, M17, M20, M21, M33, M90, C55, x, false, beta);

    // Handle edge cases with dynamic peeling
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
    if (total_steps == steps_left) {
        mkl_set_num_threads_local(num_threads);
        mkl_set_dynamic(0);
    }
#endif
    DynamicPeeling(A, B, C, 5, 5, 5, beta);
}

// C := alpha * A * B + beta * C
template <typename Scalar>
double FastMatmul(Matrix<Scalar>& A, Matrix<Scalar>& B, Matrix<Scalar>& C,
    int num_steps, double x=1e-8, int num_threads=12, Scalar alpha=Scalar(1.0), Scalar beta=Scalar(0.0)) {
    MemoryManager<Scalar> mem_mngr;
#ifdef _PARALLEL_
    mem_mngr.Allocate(5, 5, 5, 90, num_steps, A.m(), A.n(), B.n());
#endif
    A.set_multiplier(alpha);
    int num_multiplies_per_step = 90;
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

}  // namespace smirnov555_90_710_approx

#endif  // _smirnov555_90_710_approx_HPP_
