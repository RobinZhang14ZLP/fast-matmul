#ifndef _grey424_26_206_HPP_
#define _grey424_26_206_HPP_

// This is an automatically generated file from gen.py.
#include "common.hpp"

namespace grey424_26_206 {

template <typename Scalar>
void S_Add1(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& C, double x=1e-8) {
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
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataS1[i + j * strideS1] + dataS2[i + j * strideS2] + dataS3[i + j * strideS3] + dataS4[i + j * strideS4];
        }
    }
}

template <typename Scalar>
void S_Add2(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& C, double x=1e-8) {
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
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(0.5) * dataS1[i + j * strideS1] -dataS2[i + j * strideS2] + Scalar(0.5) * dataS3[i + j * strideS3] + Scalar(0.5) * dataS4[i + j * strideS4];
        }
    }
}

template <typename Scalar>
void S_Add3(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& C, double x=1e-8) {
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
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataS1[i + j * strideS1] + dataS2[i + j * strideS2] + dataS3[i + j * strideS3] + dataS4[i + j * strideS4];
        }
    }
}

template <typename Scalar>
void S_Add4(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& C, double x=1e-8) {
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
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataS1[i + j * strideS1] + dataS2[i + j * strideS2] + dataS3[i + j * strideS3] + dataS4[i + j * strideS4];
        }
    }
}

template <typename Scalar>
void S_Add5(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x=1e-8) {
    const int strideS1 = S1.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataS1[i + j * strideS1];
        }
    }
}

template <typename Scalar>
void S_Add6(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x=1e-8) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataS1[i + j * strideS1] -dataS2[i + j * strideS2] -dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add7(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x=1e-8) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataS1[i + j * strideS1] -dataS2[i + j * strideS2] + dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add8(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x=1e-8) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataS1[i + j * strideS1] + dataS2[i + j * strideS2] + dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add9(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& C, double x=1e-8) {
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
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataS1[i + j * strideS1] -dataS2[i + j * strideS2] + dataS3[i + j * strideS3] + dataS4[i + j * strideS4];
        }
    }
}

template <typename Scalar>
void S_Add10(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x=1e-8) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataS1[i + j * strideS1] -dataS2[i + j * strideS2] + dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add11(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataS1[i + j * strideS1] -dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add12(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x=1e-8) {
    const int strideS1 = S1.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataS1[i + j * strideS1];
        }
    }
}

template <typename Scalar>
void S_Add13(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x=1e-8) {
    const int strideS1 = S1.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataS1[i + j * strideS1];
        }
    }
}

template <typename Scalar>
void S_Add14(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x=1e-8) {
    const int strideS1 = S1.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataS1[i + j * strideS1];
        }
    }
}

template <typename Scalar>
void S_Add15(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataS1[i + j * strideS1] + dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add16(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x=1e-8) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataS1[i + j * strideS1] -dataS2[i + j * strideS2] + dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add17(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x=1e-8) {
    const int strideS1 = S1.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataS1[i + j * strideS1];
        }
    }
}

template <typename Scalar>
void S_Add18(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x=1e-8) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataS1[i + j * strideS1] -dataS2[i + j * strideS2] + dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add19(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x=1e-8) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataS1[i + j * strideS1] -dataS2[i + j * strideS2] -dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add20(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x=1e-8) {
    const int strideS1 = S1.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataS1[i + j * strideS1];
        }
    }
}

template <typename Scalar>
void S_Add21(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& C, double x=1e-8) {
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
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataS1[i + j * strideS1] + dataS2[i + j * strideS2] + dataS3[i + j * strideS3] + dataS4[i + j * strideS4];
        }
    }
}

template <typename Scalar>
void S_Add22(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x=1e-8) {
    const int strideS1 = S1.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataS1[i + j * strideS1];
        }
    }
}

template <typename Scalar>
void S_Add23(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& C, double x=1e-8) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideS3 = S3.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    const Scalar *dataS3 = S3.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataS1[i + j * strideS1] -dataS2[i + j * strideS2] -dataS3[i + j * strideS3];
        }
    }
}

template <typename Scalar>
void S_Add24(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& C, double x=1e-8) {
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
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataS1[i + j * strideS1] + dataS2[i + j * strideS2] -dataS3[i + j * strideS3] -dataS4[i + j * strideS4];
        }
    }
}

template <typename Scalar>
void S_Add25(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& S3, Matrix<Scalar>& S4, Matrix<Scalar>& C, double x=1e-8) {
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
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(-0.5) * dataS1[i + j * strideS1] + Scalar(0.5) * dataS2[i + j * strideS2] + Scalar(-0.5) * dataS3[i + j * strideS3] -dataS4[i + j * strideS4];
        }
    }
}

template <typename Scalar>
void S_Add26(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideS1 = S1.stride();
    const int strideS2 = S2.stride();
    const int strideC = C.stride();
    const Scalar *dataS1 = S1.data();
    const Scalar *dataS2 = S2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataS1[i + j * strideS1] -dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void T_Add1(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x=1e-8) {
    const int strideT1 = T1.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataT1[i + j * strideT1];
        }
    }
}

template <typename Scalar>
void T_Add2(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataT1[i + j * strideT1] + dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add3(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataT1[i + j * strideT1] + dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add4(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataT1[i + j * strideT1] + dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add5(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x=1e-8) {
    const int strideT1 = T1.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataT1[i + j * strideT1];
        }
    }
}

template <typename Scalar>
void T_Add6(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataT1[i + j * strideT1] + dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add7(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& C, double x=1e-8) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideT3 = T3.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    const Scalar *dataT3 = T3.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataT1[i + j * strideT1] + dataT2[i + j * strideT2] + dataT3[i + j * strideT3];
        }
    }
}

template <typename Scalar>
void T_Add8(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataT1[i + j * strideT1] + dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add9(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x=1e-8) {
    const int strideT1 = T1.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataT1[i + j * strideT1];
        }
    }
}

template <typename Scalar>
void T_Add10(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& C, double x=1e-8) {
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
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataT1[i + j * strideT1] + dataT2[i + j * strideT2] + dataT3[i + j * strideT3] + dataT4[i + j * strideT4];
        }
    }
}

template <typename Scalar>
void T_Add11(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataT1[i + j * strideT1] -dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add12(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x=1e-8) {
    const int strideT1 = T1.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataT1[i + j * strideT1];
        }
    }
}

template <typename Scalar>
void T_Add13(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x=1e-8) {
    const int strideT1 = T1.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataT1[i + j * strideT1];
        }
    }
}

template <typename Scalar>
void T_Add14(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x=1e-8) {
    const int strideT1 = T1.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataT1[i + j * strideT1];
        }
    }
}

template <typename Scalar>
void T_Add15(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& C, double x=1e-8) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideT3 = T3.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    const Scalar *dataT3 = T3.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataT1[i + j * strideT1] -dataT2[i + j * strideT2] + dataT3[i + j * strideT3];
        }
    }
}

template <typename Scalar>
void T_Add16(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& C, double x=1e-8) {
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
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataT1[i + j * strideT1] + dataT2[i + j * strideT2] -dataT3[i + j * strideT3] + dataT4[i + j * strideT4];
        }
    }
}

template <typename Scalar>
void T_Add17(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& C, double x=1e-8) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideT3 = T3.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    const Scalar *dataT3 = T3.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataT1[i + j * strideT1] + dataT2[i + j * strideT2] + dataT3[i + j * strideT3];
        }
    }
}

template <typename Scalar>
void T_Add18(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& C, double x=1e-8) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideT3 = T3.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    const Scalar *dataT3 = T3.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataT1[i + j * strideT1] + dataT2[i + j * strideT2] + dataT3[i + j * strideT3];
        }
    }
}

template <typename Scalar>
void T_Add19(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x=1e-8) {
    const int strideT1 = T1.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataT1[i + j * strideT1];
        }
    }
}

template <typename Scalar>
void T_Add20(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& T4, Matrix<Scalar>& C, double x=1e-8) {
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
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataT1[i + j * strideT1] + dataT2[i + j * strideT2] -dataT3[i + j * strideT3] -dataT4[i + j * strideT4];
        }
    }
}

template <typename Scalar>
void T_Add21(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& T3, Matrix<Scalar>& C, double x=1e-8) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideT3 = T3.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    const Scalar *dataT3 = T3.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataT1[i + j * strideT1] + dataT2[i + j * strideT2] + dataT3[i + j * strideT3];
        }
    }
}

template <typename Scalar>
void T_Add22(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataT1[i + j * strideT1] -dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add23(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x=1e-8) {
    const int strideT1 = T1.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataT1[i + j * strideT1];
        }
    }
}

template <typename Scalar>
void T_Add24(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataT1[i + j * strideT1] -dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add25(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideT1 = T1.stride();
    const int strideT2 = T2.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    const Scalar *dataT2 = T2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataT1[i + j * strideT1] + dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add26(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x=1e-8) {
    const int strideT1 = T1.stride();
    const int strideC = C.stride();
    const Scalar *dataT1 = T1.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataT1[i + j * strideT1];
        }
    }
}

template <typename Scalar>
void M_Add1(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& C, double x=1e-8, Scalar beta=Scalar(0.0)) {
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
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = dataM1[i + j * strideM1] -dataM2[i + j * strideM2] -dataM3[i + j * strideM3] -dataM4[i + j * strideM4] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = dataM1[i + j * strideM1] -dataM2[i + j * strideM2] -dataM3[i + j * strideM3] -dataM4[i + j * strideM4];
            }
        }
    }
}

template <typename Scalar>
void M_Add2(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& C, double x=1e-8, Scalar beta=Scalar(0.0)) {
    const int strideM1 = M1.stride();
    const int strideM2 = M2.stride();
    const int strideC = C.stride();
    const Scalar *dataM1 = M1.data();
    const Scalar *dataM2 = M2.data();
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = dataM1[i + j * strideM1] + dataM2[i + j * strideM2] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = dataM1[i + j * strideM1] + dataM2[i + j * strideM2];
            }
        }
    }
}

template <typename Scalar>
void M_Add3(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& C, double x=1e-8, Scalar beta=Scalar(0.0)) {
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
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = dataM1[i + j * strideM1] -dataM2[i + j * strideM2] -dataM3[i + j * strideM3] + dataM4[i + j * strideM4] + dataM5[i + j * strideM5] + dataM6[i + j * strideM6] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = dataM1[i + j * strideM1] -dataM2[i + j * strideM2] -dataM3[i + j * strideM3] + dataM4[i + j * strideM4] + dataM5[i + j * strideM5] + dataM6[i + j * strideM6];
            }
        }
    }
}

template <typename Scalar>
void M_Add4(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& C, double x=1e-8, Scalar beta=Scalar(0.0)) {
    const int strideM1 = M1.stride();
    const int strideM2 = M2.stride();
    const int strideC = C.stride();
    const Scalar *dataM1 = M1.data();
    const Scalar *dataM2 = M2.data();
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = -dataM1[i + j * strideM1] + dataM2[i + j * strideM2] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = -dataM1[i + j * strideM1] + dataM2[i + j * strideM2];
            }
        }
    }
}

template <typename Scalar>
void M_Add5(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& C, double x=1e-8, Scalar beta=Scalar(0.0)) {
    const int strideM1 = M1.stride();
    const int strideM2 = M2.stride();
    const int strideC = C.stride();
    const Scalar *dataM1 = M1.data();
    const Scalar *dataM2 = M2.data();
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = dataM1[i + j * strideM1] + dataM2[i + j * strideM2] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = dataM1[i + j * strideM1] + dataM2[i + j * strideM2];
            }
        }
    }
}

template <typename Scalar>
void M_Add6(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& C, double x=1e-8, Scalar beta=Scalar(0.0)) {
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
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = -dataM1[i + j * strideM1] + dataM2[i + j * strideM2] -dataM3[i + j * strideM3] + dataM4[i + j * strideM4] + dataM5[i + j * strideM5] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = -dataM1[i + j * strideM1] + dataM2[i + j * strideM2] -dataM3[i + j * strideM3] + dataM4[i + j * strideM4] + dataM5[i + j * strideM5];
            }
        }
    }
}

template <typename Scalar>
void M_Add7(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& C, double x=1e-8, Scalar beta=Scalar(0.0)) {
    const int strideM1 = M1.stride();
    const int strideM2 = M2.stride();
    const int strideM3 = M3.stride();
    const int strideC = C.stride();
    const Scalar *dataM1 = M1.data();
    const Scalar *dataM2 = M2.data();
    const Scalar *dataM3 = M3.data();
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(-0.5) * dataM1[i + j * strideM1] -dataM2[i + j * strideM2] + dataM3[i + j * strideM3] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(-0.5) * dataM1[i + j * strideM1] -dataM2[i + j * strideM2] + dataM3[i + j * strideM3];
            }
        }
    }
}

template <typename Scalar>
void M_Add8(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& M6, Matrix<Scalar>& M7, Matrix<Scalar>& C, double x=1e-8, Scalar beta=Scalar(0.0)) {
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
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(-0.5) * dataM1[i + j * strideM1] + dataM2[i + j * strideM2] + Scalar(0.5) * dataM3[i + j * strideM3] + dataM4[i + j * strideM4] + dataM5[i + j * strideM5] + dataM6[i + j * strideM6] + dataM7[i + j * strideM7] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(-0.5) * dataM1[i + j * strideM1] + dataM2[i + j * strideM2] + Scalar(0.5) * dataM3[i + j * strideM3] + dataM4[i + j * strideM4] + dataM5[i + j * strideM5] + dataM6[i + j * strideM6] + dataM7[i + j * strideM7];
            }
        }
    }
}

template <typename Scalar>
void M_Add9(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& C, double x=1e-8, Scalar beta=Scalar(0.0)) {
    const int strideM1 = M1.stride();
    const int strideM2 = M2.stride();
    const int strideM3 = M3.stride();
    const int strideC = C.stride();
    const Scalar *dataM1 = M1.data();
    const Scalar *dataM2 = M2.data();
    const Scalar *dataM3 = M3.data();
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(0.5) * dataM1[i + j * strideM1] + Scalar(0.5) * dataM2[i + j * strideM2] + dataM3[i + j * strideM3] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(0.5) * dataM1[i + j * strideM1] + Scalar(0.5) * dataM2[i + j * strideM2] + dataM3[i + j * strideM3];
            }
        }
    }
}

template <typename Scalar>
void M_Add10(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& C, double x=1e-8, Scalar beta=Scalar(0.0)) {
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
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = dataM1[i + j * strideM1] + dataM2[i + j * strideM2] + dataM3[i + j * strideM3] + dataM4[i + j * strideM4] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = dataM1[i + j * strideM1] + dataM2[i + j * strideM2] + dataM3[i + j * strideM3] + dataM4[i + j * strideM4];
            }
        }
    }
}

template <typename Scalar>
void M_Add11(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& C, double x=1e-8, Scalar beta=Scalar(0.0)) {
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
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = -dataM1[i + j * strideM1] + dataM2[i + j * strideM2] + dataM3[i + j * strideM3] + dataM4[i + j * strideM4] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = -dataM1[i + j * strideM1] + dataM2[i + j * strideM2] + dataM3[i + j * strideM3] + dataM4[i + j * strideM4];
            }
        }
    }
}

template <typename Scalar>
void M_Add12(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& C, double x=1e-8, Scalar beta=Scalar(0.0)) {
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
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(-0.5) * dataM1[i + j * strideM1] + Scalar(-0.5) * dataM2[i + j * strideM2] + dataM3[i + j * strideM3] -dataM4[i + j * strideM4] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(-0.5) * dataM1[i + j * strideM1] + Scalar(-0.5) * dataM2[i + j * strideM2] + dataM3[i + j * strideM3] -dataM4[i + j * strideM4];
            }
        }
    }
}

template <typename Scalar>
void M_Add13(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& M5, Matrix<Scalar>& C, double x=1e-8, Scalar beta=Scalar(0.0)) {
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
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = -dataM1[i + j * strideM1] + dataM2[i + j * strideM2] + Scalar(-0.5) * dataM3[i + j * strideM3] -dataM4[i + j * strideM4] -dataM5[i + j * strideM5] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = -dataM1[i + j * strideM1] + dataM2[i + j * strideM2] + Scalar(-0.5) * dataM3[i + j * strideM3] -dataM4[i + j * strideM4] -dataM5[i + j * strideM5];
            }
        }
    }
}

template <typename Scalar>
void M_Add14(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& C, double x=1e-8, Scalar beta=Scalar(0.0)) {
    const int strideM1 = M1.stride();
    const int strideM2 = M2.stride();
    const int strideC = C.stride();
    const Scalar *dataM1 = M1.data();
    const Scalar *dataM2 = M2.data();
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = dataM1[i + j * strideM1] + dataM2[i + j * strideM2] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = dataM1[i + j * strideM1] + dataM2[i + j * strideM2];
            }
        }
    }
}

template <typename Scalar>
void M_Add15(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& C, double x=1e-8, Scalar beta=Scalar(0.0)) {
    const int strideM1 = M1.stride();
    const int strideM2 = M2.stride();
    const int strideM3 = M3.stride();
    const int strideC = C.stride();
    const Scalar *dataM1 = M1.data();
    const Scalar *dataM2 = M2.data();
    const Scalar *dataM3 = M3.data();
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = dataM1[i + j * strideM1] + Scalar(0.5) * dataM2[i + j * strideM2] -dataM3[i + j * strideM3] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = dataM1[i + j * strideM1] + Scalar(0.5) * dataM2[i + j * strideM2] -dataM3[i + j * strideM3];
            }
        }
    }
}

template <typename Scalar>
void M_Add16(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& C, double x=1e-8, Scalar beta=Scalar(0.0)) {
    const int strideM1 = M1.stride();
    const int strideM2 = M2.stride();
    const int strideC = C.stride();
    const Scalar *dataM1 = M1.data();
    const Scalar *dataM2 = M2.data();
    Scalar *dataC = C.data();
    if (beta != Scalar(0.0)) {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = -dataM1[i + j * strideM1] -dataM2[i + j * strideM2] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = -dataM1[i + j * strideM1] -dataM2[i + j * strideM2];
            }
        }
    }
}

template <typename Scalar>
void SX_Add1(Matrix<Scalar>& SX1, Matrix<Scalar>& SX2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideSX1 = SX1.stride();
    const int strideSX2 = SX2.stride();
    const int strideC = C.stride();
    const Scalar *dataSX1 = SX1.data();
    const Scalar *dataSX2 = SX2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataSX1[i + j * strideSX1] + dataSX2[i + j * strideSX2];
        }
    }
}

template <typename Scalar>
void SX_Add2(Matrix<Scalar>& SX1, Matrix<Scalar>& SX2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideSX1 = SX1.stride();
    const int strideSX2 = SX2.stride();
    const int strideC = C.stride();
    const Scalar *dataSX1 = SX1.data();
    const Scalar *dataSX2 = SX2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataSX1[i + j * strideSX1] -dataSX2[i + j * strideSX2];
        }
    }
}

template <typename Scalar>
void SX_Add3(Matrix<Scalar>& SX1, Matrix<Scalar>& SX2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideSX1 = SX1.stride();
    const int strideSX2 = SX2.stride();
    const int strideC = C.stride();
    const Scalar *dataSX1 = SX1.data();
    const Scalar *dataSX2 = SX2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataSX1[i + j * strideSX1] -dataSX2[i + j * strideSX2];
        }
    }
}

template <typename Scalar>
void SX_Add4(Matrix<Scalar>& SX1, Matrix<Scalar>& SX2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideSX1 = SX1.stride();
    const int strideSX2 = SX2.stride();
    const int strideC = C.stride();
    const Scalar *dataSX1 = SX1.data();
    const Scalar *dataSX2 = SX2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataSX1[i + j * strideSX1] -dataSX2[i + j * strideSX2];
        }
    }
}

template <typename Scalar>
void SX_Add5(Matrix<Scalar>& SX1, Matrix<Scalar>& SX2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideSX1 = SX1.stride();
    const int strideSX2 = SX2.stride();
    const int strideC = C.stride();
    const Scalar *dataSX1 = SX1.data();
    const Scalar *dataSX2 = SX2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataSX1[i + j * strideSX1] + dataSX2[i + j * strideSX2];
        }
    }
}

template <typename Scalar>
void SX_Add6(Matrix<Scalar>& SX1, Matrix<Scalar>& SX2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideSX1 = SX1.stride();
    const int strideSX2 = SX2.stride();
    const int strideC = C.stride();
    const Scalar *dataSX1 = SX1.data();
    const Scalar *dataSX2 = SX2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataSX1[i + j * strideSX1] -dataSX2[i + j * strideSX2];
        }
    }
}

template <typename Scalar>
void TX_Add1(Matrix<Scalar>& TX1, Matrix<Scalar>& TX2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideTX1 = TX1.stride();
    const int strideTX2 = TX2.stride();
    const int strideC = C.stride();
    const Scalar *dataTX1 = TX1.data();
    const Scalar *dataTX2 = TX2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataTX1[i + j * strideTX1] + dataTX2[i + j * strideTX2];
        }
    }
}

template <typename Scalar>
void TX_Add2(Matrix<Scalar>& TX1, Matrix<Scalar>& TX2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideTX1 = TX1.stride();
    const int strideTX2 = TX2.stride();
    const int strideC = C.stride();
    const Scalar *dataTX1 = TX1.data();
    const Scalar *dataTX2 = TX2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataTX1[i + j * strideTX1] + dataTX2[i + j * strideTX2];
        }
    }
}

template <typename Scalar>
void TX_Add3(Matrix<Scalar>& TX1, Matrix<Scalar>& TX2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideTX1 = TX1.stride();
    const int strideTX2 = TX2.stride();
    const int strideC = C.stride();
    const Scalar *dataTX1 = TX1.data();
    const Scalar *dataTX2 = TX2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataTX1[i + j * strideTX1] -dataTX2[i + j * strideTX2];
        }
    }
}

template <typename Scalar>
void TX_Add4(Matrix<Scalar>& TX1, Matrix<Scalar>& TX2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideTX1 = TX1.stride();
    const int strideTX2 = TX2.stride();
    const int strideC = C.stride();
    const Scalar *dataTX1 = TX1.data();
    const Scalar *dataTX2 = TX2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataTX1[i + j * strideTX1] -dataTX2[i + j * strideTX2];
        }
    }
}

template <typename Scalar>
void TX_Add5(Matrix<Scalar>& TX1, Matrix<Scalar>& TX2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideTX1 = TX1.stride();
    const int strideTX2 = TX2.stride();
    const int strideC = C.stride();
    const Scalar *dataTX1 = TX1.data();
    const Scalar *dataTX2 = TX2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataTX1[i + j * strideTX1] -dataTX2[i + j * strideTX2];
        }
    }
}

template <typename Scalar>
void TX_Add6(Matrix<Scalar>& TX1, Matrix<Scalar>& TX2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideTX1 = TX1.stride();
    const int strideTX2 = TX2.stride();
    const int strideC = C.stride();
    const Scalar *dataTX1 = TX1.data();
    const Scalar *dataTX2 = TX2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataTX1[i + j * strideTX1] + dataTX2[i + j * strideTX2];
        }
    }
}

template <typename Scalar>
void MX_Add1(Matrix<Scalar>& MX1, Matrix<Scalar>& MX2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideMX1 = MX1.stride();
    const int strideMX2 = MX2.stride();
    const int strideC = C.stride();
    const Scalar *dataMX1 = MX1.data();
    const Scalar *dataMX2 = MX2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataMX1[i + j * strideMX1] + dataMX2[i + j * strideMX2];
        }
    }
}

template <typename Scalar>
void MX_Add2(Matrix<Scalar>& MX1, Matrix<Scalar>& MX2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideMX1 = MX1.stride();
    const int strideMX2 = MX2.stride();
    const int strideC = C.stride();
    const Scalar *dataMX1 = MX1.data();
    const Scalar *dataMX2 = MX2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataMX1[i + j * strideMX1] + dataMX2[i + j * strideMX2];
        }
    }
}

template <typename Scalar>
void MX_Add3(Matrix<Scalar>& MX1, Matrix<Scalar>& MX2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideMX1 = MX1.stride();
    const int strideMX2 = MX2.stride();
    const int strideC = C.stride();
    const Scalar *dataMX1 = MX1.data();
    const Scalar *dataMX2 = MX2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataMX1[i + j * strideMX1] + Scalar(0.5) * dataMX2[i + j * strideMX2];
        }
    }
}

template <typename Scalar>
void MX_Add4(Matrix<Scalar>& MX1, Matrix<Scalar>& MX2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideMX1 = MX1.stride();
    const int strideMX2 = MX2.stride();
    const int strideC = C.stride();
    const Scalar *dataMX1 = MX1.data();
    const Scalar *dataMX2 = MX2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(-0.5) * dataMX1[i + j * strideMX1] + dataMX2[i + j * strideMX2];
        }
    }
}

template <typename Scalar>
void MX_Add5(Matrix<Scalar>& MX1, Matrix<Scalar>& MX2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideMX1 = MX1.stride();
    const int strideMX2 = MX2.stride();
    const int strideC = C.stride();
    const Scalar *dataMX1 = MX1.data();
    const Scalar *dataMX2 = MX2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataMX1[i + j * strideMX1] + Scalar(0.5) * dataMX2[i + j * strideMX2];
        }
    }
}

template <typename Scalar>
void MX_Add6(Matrix<Scalar>& MX1, Matrix<Scalar>& MX2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideMX1 = MX1.stride();
    const int strideMX2 = MX2.stride();
    const int strideC = C.stride();
    const Scalar *dataMX1 = MX1.data();
    const Scalar *dataMX2 = MX2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(0.5) * dataMX1[i + j * strideMX1] -dataMX2[i + j * strideMX2];
        }
    }
}

template <typename Scalar>
void MX_Add7(Matrix<Scalar>& MX1, Matrix<Scalar>& MX2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideMX1 = MX1.stride();
    const int strideMX2 = MX2.stride();
    const int strideC = C.stride();
    const Scalar *dataMX1 = MX1.data();
    const Scalar *dataMX2 = MX2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(0.5) * dataMX1[i + j * strideMX1] + Scalar(-0.5) * dataMX2[i + j * strideMX2];
        }
    }
}

template <typename Scalar>
void MX_Add8(Matrix<Scalar>& MX1, Matrix<Scalar>& MX2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideMX1 = MX1.stride();
    const int strideMX2 = MX2.stride();
    const int strideC = C.stride();
    const Scalar *dataMX1 = MX1.data();
    const Scalar *dataMX2 = MX2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = dataMX1[i + j * strideMX1] -dataMX2[i + j * strideMX2];
        }
    }
}

template <typename Scalar>
void MX_Add9(Matrix<Scalar>& MX1, Matrix<Scalar>& MX2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideMX1 = MX1.stride();
    const int strideMX2 = MX2.stride();
    const int strideC = C.stride();
    const Scalar *dataMX1 = MX1.data();
    const Scalar *dataMX2 = MX2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(-0.5) * dataMX1[i + j * strideMX1] + Scalar(0.5) * dataMX2[i + j * strideMX2];
        }
    }
}

template <typename Scalar>
void MX_Add10(Matrix<Scalar>& MX1, Matrix<Scalar>& MX2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideMX1 = MX1.stride();
    const int strideMX2 = MX2.stride();
    const int strideC = C.stride();
    const Scalar *dataMX1 = MX1.data();
    const Scalar *dataMX2 = MX2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataMX1[i + j * strideMX1] -dataMX2[i + j * strideMX2];
        }
    }
}

template <typename Scalar>
void MX_Add11(Matrix<Scalar>& MX1, Matrix<Scalar>& MX2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideMX1 = MX1.stride();
    const int strideMX2 = MX2.stride();
    const int strideC = C.stride();
    const Scalar *dataMX1 = MX1.data();
    const Scalar *dataMX2 = MX2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(-0.5) * dataMX1[i + j * strideMX1] + Scalar(0.5) * dataMX2[i + j * strideMX2];
        }
    }
}

template <typename Scalar>
void MX_Add12(Matrix<Scalar>& MX1, Matrix<Scalar>& MX2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideMX1 = MX1.stride();
    const int strideMX2 = MX2.stride();
    const int strideC = C.stride();
    const Scalar *dataMX1 = MX1.data();
    const Scalar *dataMX2 = MX2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = -dataMX1[i + j * strideMX1] + dataMX2[i + j * strideMX2];
        }
    }
}

template <typename Scalar>
void MX_Add13(Matrix<Scalar>& MX1, Matrix<Scalar>& MX2, Matrix<Scalar>& C, double x=1e-8) {
    const int strideMX1 = MX1.stride();
    const int strideMX2 = MX2.stride();
    const int strideC = C.stride();
    const Scalar *dataMX1 = MX1.data();
    const Scalar *dataMX2 = MX2.data();
    Scalar *dataC = C.data();
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
    for (int j = 0; j < C.n(); ++j) {
        for (int i = 0; i < C.m(); ++i) {
            dataC[i + j * strideC] = Scalar(0.5) * dataMX1[i + j * strideMX1] + dataMX2[i + j * strideMX2];
        }
    }
}

template <typename Scalar>
void FastMatmulRecursive(Matrix<Scalar>& A, Matrix<Scalar>& B, Matrix<Scalar>& C, int total_steps, int steps_left, int start_index, double x, int num_threads, Scalar beta) {
    // Update multipliers
    C.UpdateMultiplier(A.multiplier());
    C.UpdateMultiplier(B.multiplier());
    A.set_multiplier(Scalar(1.0));
    B.set_multiplier(Scalar(1.0));
    // Base case for recursion
    if (steps_left == 0) {
        Gemm(A, B, C);
        return;
    }

    Matrix<Scalar> A11 = A.Subblock(4, 2, 1, 1);
    Matrix<Scalar> A12 = A.Subblock(4, 2, 1, 2);
    Matrix<Scalar> A21 = A.Subblock(4, 2, 2, 1);
    Matrix<Scalar> A22 = A.Subblock(4, 2, 2, 2);
    Matrix<Scalar> A31 = A.Subblock(4, 2, 3, 1);
    Matrix<Scalar> A32 = A.Subblock(4, 2, 3, 2);
    Matrix<Scalar> A41 = A.Subblock(4, 2, 4, 1);
    Matrix<Scalar> A42 = A.Subblock(4, 2, 4, 2);
    Matrix<Scalar> B11 = B.Subblock(2, 4, 1, 1);
    Matrix<Scalar> B12 = B.Subblock(2, 4, 1, 2);
    Matrix<Scalar> B13 = B.Subblock(2, 4, 1, 3);
    Matrix<Scalar> B14 = B.Subblock(2, 4, 1, 4);
    Matrix<Scalar> B21 = B.Subblock(2, 4, 2, 1);
    Matrix<Scalar> B22 = B.Subblock(2, 4, 2, 2);
    Matrix<Scalar> B23 = B.Subblock(2, 4, 2, 3);
    Matrix<Scalar> B24 = B.Subblock(2, 4, 2, 4);
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
    Matrix<Scalar> A_X1(A11.m(), A11.n());
    SX_Add1(A11, A21, A_X1);
    Matrix<Scalar> A_X2(A11.m(), A11.n());
    SX_Add2(A12, A22, A_X2);
    Matrix<Scalar> A_X3(A11.m(), A11.n());
    SX_Add3(A31, A32, A_X3);
    Matrix<Scalar> A_X4(A11.m(), A11.n());
    SX_Add4(A41, A42, A_X4);
    Matrix<Scalar> A_X5(A11.m(), A11.n());
    SX_Add5(A31, A32, A_X5);
    Matrix<Scalar> A_X6(A11.m(), A11.n());
    SX_Add6(A41, A42, A_X6);

    Matrix<Scalar> B_X1(B11.m(), B11.n());
    TX_Add1(B11, B14, B_X1);
    Matrix<Scalar> B_X2(B11.m(), B11.n());
    TX_Add2(B11, B13, B_X2);
    Matrix<Scalar> B_X3(B11.m(), B11.n());
    TX_Add3(B12, B13, B_X3);
    Matrix<Scalar> B_X4(B11.m(), B11.n());
    TX_Add4(B12, B22, B_X4);
    Matrix<Scalar> B_X5(B11.m(), B11.n());
    TX_Add5(B14, B22, B_X5);
    Matrix<Scalar> B_X6(B11.m(), B11.n());
    TX_Add6(B21, B23, B_X6);


    // M1 = (1.0 * A_X1 + 1.0 * A_X2 + 1.0 * A_X3 + 1.0 * A_X4) * (1.0 * B_X1)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(26, total_steps, steps_left, start_index, 1, num_threads)) untied
    {
#endif
    Matrix<Scalar> S1(A11.m(), A11.n());
    S_Add1(A_X1, A_X2, A_X3, A_X4, S1, x);
    FastMatmulRecursive(S1, B_X1, M1, total_steps, steps_left - 1, (start_index + 1 - 1) * 26, x, num_threads, Scalar(0.0));
    S1.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(26, total_steps, steps_left, start_index, 1, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M2 = (0.5 * A21 + -1.0 * A22 + 0.5 * A31 + 0.5 * A41) * (-1.0 * B21 + 1.0 * B_X2)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(26, total_steps, steps_left, start_index, 2, num_threads)) untied
    {
#endif
    Matrix<Scalar> S2(A11.m(), A11.n());
    S_Add2(A21, A22, A31, A41, S2, x);
    Matrix<Scalar> T2(B11.m(), B11.n());
    T_Add2(B21, B_X2, T2, x);
    FastMatmulRecursive(S2, T2, M2, total_steps, steps_left - 1, (start_index + 2 - 1) * 26, x, num_threads, Scalar(0.0));
    S2.deallocate();
    T2.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(26, total_steps, steps_left, start_index, 2, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M3 = (-1.0 * A_X1 + 1.0 * A_X2 + 1.0 * A_X5 + 1.0 * A_X6) * (-1.0 * B14 + 1.0 * B_X3)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(26, total_steps, steps_left, start_index, 3, num_threads)) untied
    {
#endif
    Matrix<Scalar> S3(A11.m(), A11.n());
    S_Add3(A_X1, A_X2, A_X5, A_X6, S3, x);
    Matrix<Scalar> T3(B11.m(), B11.n());
    T_Add3(B14, B_X3, T3, x);
    FastMatmulRecursive(S3, T3, M3, total_steps, steps_left - 1, (start_index + 3 - 1) * 26, x, num_threads, Scalar(0.0));
    S3.deallocate();
    T3.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(26, total_steps, steps_left, start_index, 3, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M4 = (1.0 * A21 + 1.0 * A_X2 + 1.0 * A_X3 + 1.0 * A_X4) * (1.0 * B21 + 1.0 * B_X1)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(26, total_steps, steps_left, start_index, 4, num_threads)) untied
    {
#endif
    Matrix<Scalar> S4(A11.m(), A11.n());
    S_Add4(A21, A_X2, A_X3, A_X4, S4, x);
    Matrix<Scalar> T4(B11.m(), B11.n());
    T_Add4(B21, B_X1, T4, x);
    FastMatmulRecursive(S4, T4, M4, total_steps, steps_left - 1, (start_index + 4 - 1) * 26, x, num_threads, Scalar(0.0));
    S4.deallocate();
    T4.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(26, total_steps, steps_left, start_index, 4, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M5 = (-1.0 * A42) * (1.0 * B_X4)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(26, total_steps, steps_left, start_index, 5, num_threads)) untied
    {
#endif
    M5.UpdateMultiplier(Scalar(-1));
    FastMatmulRecursive(A42, B_X4, M5, total_steps, steps_left - 1, (start_index + 5 - 1) * 26, x, num_threads, Scalar(0.0));
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(26, total_steps, steps_left, start_index, 5, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M6 = (-1.0 * A11 + -1.0 * A12 + -1.0 * A_X6) * (1.0 * B22 + 1.0 * B24)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(26, total_steps, steps_left, start_index, 6, num_threads)) untied
    {
#endif
    Matrix<Scalar> S6(A11.m(), A11.n());
    S_Add6(A11, A12, A_X6, S6, x);
    Matrix<Scalar> T6(B11.m(), B11.n());
    T_Add6(B22, B24, T6, x);
    FastMatmulRecursive(S6, T6, M6, total_steps, steps_left - 1, (start_index + 6 - 1) * 26, x, num_threads, Scalar(0.0));
    S6.deallocate();
    T6.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(26, total_steps, steps_left, start_index, 6, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M7 = (-1.0 * A32 + -1.0 * A42 + 1.0 * A_X2) * (1.0 * B21 + 1.0 * B24 + 1.0 * B_X1)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(26, total_steps, steps_left, start_index, 7, num_threads)) untied
    {
#endif
    Matrix<Scalar> S7(A11.m(), A11.n());
    S_Add7(A32, A42, A_X2, S7, x);
    Matrix<Scalar> T7(B11.m(), B11.n());
    T_Add7(B21, B24, B_X1, T7, x);
    FastMatmulRecursive(S7, T7, M7, total_steps, steps_left - 1, (start_index + 7 - 1) * 26, x, num_threads, Scalar(0.0));
    S7.deallocate();
    T7.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(26, total_steps, steps_left, start_index, 7, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M8 = (1.0 * A21 + 1.0 * A31 + 1.0 * A41) * (1.0 * B11 + 1.0 * B21)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(26, total_steps, steps_left, start_index, 8, num_threads)) untied
    {
#endif
    Matrix<Scalar> S8(A11.m(), A11.n());
    S_Add8(A21, A31, A41, S8, x);
    Matrix<Scalar> T8(B11.m(), B11.n());
    T_Add8(B11, B21, T8, x);
    FastMatmulRecursive(S8, T8, M8, total_steps, steps_left - 1, (start_index + 8 - 1) * 26, x, num_threads, Scalar(0.0));
    S8.deallocate();
    T8.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(26, total_steps, steps_left, start_index, 8, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M9 = (1.0 * A21 + -1.0 * A22 + 1.0 * A_X3 + 1.0 * A_X4) * (-1.0 * B21)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(26, total_steps, steps_left, start_index, 9, num_threads)) untied
    {
#endif
    Matrix<Scalar> S9(A11.m(), A11.n());
    S_Add9(A21, A22, A_X3, A_X4, S9, x);
    M9.UpdateMultiplier(Scalar(-1));
    FastMatmulRecursive(S9, B21, M9, total_steps, steps_left - 1, (start_index + 9 - 1) * 26, x, num_threads, Scalar(0.0));
    S9.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(26, total_steps, steps_left, start_index, 9, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M10 = (1.0 * A32 + -1.0 * A42 + 1.0 * A_X2) * (1.0 * B23 + 1.0 * B24 + 1.0 * B_X3 + 1.0 * B_X5)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(26, total_steps, steps_left, start_index, 10, num_threads)) untied
    {
#endif
    Matrix<Scalar> S10(A11.m(), A11.n());
    S_Add10(A32, A42, A_X2, S10, x);
    Matrix<Scalar> T10(B11.m(), B11.n());
    T_Add10(B23, B24, B_X3, B_X5, T10, x);
    FastMatmulRecursive(S10, T10, M10, total_steps, steps_left - 1, (start_index + 10 - 1) * 26, x, num_threads, Scalar(0.0));
    S10.deallocate();
    T10.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(26, total_steps, steps_left, start_index, 10, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M11 = (-1.0 * A11 + -1.0 * A_X6) * (1.0 * B24 + -1.0 * B_X4)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(26, total_steps, steps_left, start_index, 11, num_threads)) untied
    {
#endif
    Matrix<Scalar> S11(A11.m(), A11.n());
    S_Add11(A11, A_X6, S11, x);
    Matrix<Scalar> T11(B11.m(), B11.n());
    T_Add11(B24, B_X4, T11, x);
    FastMatmulRecursive(S11, T11, M11, total_steps, steps_left - 1, (start_index + 11 - 1) * 26, x, num_threads, Scalar(0.0));
    S11.deallocate();
    T11.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(26, total_steps, steps_left, start_index, 11, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M12 = (-1.0 * A12) * (1.0 * B24)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(26, total_steps, steps_left, start_index, 12, num_threads)) untied
    {
#endif
    M12.UpdateMultiplier(Scalar(-1));
    FastMatmulRecursive(A12, B24, M12, total_steps, steps_left - 1, (start_index + 12 - 1) * 26, x, num_threads, Scalar(0.0));
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(26, total_steps, steps_left, start_index, 12, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M13 = (1.0 * A_X6) * (-1.0 * B12)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(26, total_steps, steps_left, start_index, 13, num_threads)) untied
    {
#endif
    M13.UpdateMultiplier(Scalar(-1));
    FastMatmulRecursive(A_X6, B12, M13, total_steps, steps_left - 1, (start_index + 13 - 1) * 26, x, num_threads, Scalar(0.0));
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(26, total_steps, steps_left, start_index, 13, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M14 = (-1.0 * A11) * (-1.0 * B14)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(26, total_steps, steps_left, start_index, 14, num_threads)) untied
    {
#endif
    M14.UpdateMultiplier(Scalar(-1));
    M14.UpdateMultiplier(Scalar(-1));
    FastMatmulRecursive(A11, B14, M14, total_steps, steps_left - 1, (start_index + 14 - 1) * 26, x, num_threads, Scalar(0.0));
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(26, total_steps, steps_left, start_index, 14, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M15 = (-1.0 * A11 + 1.0 * A41) * (1.0 * B14 + -1.0 * B24 + 1.0 * B_X4)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(26, total_steps, steps_left, start_index, 15, num_threads)) untied
    {
#endif
    Matrix<Scalar> S15(A11.m(), A11.n());
    S_Add15(A11, A41, S15, x);
    Matrix<Scalar> T15(B11.m(), B11.n());
    T_Add15(B14, B24, B_X4, T15, x);
    FastMatmulRecursive(S15, T15, M15, total_steps, steps_left - 1, (start_index + 15 - 1) * 26, x, num_threads, Scalar(0.0));
    S15.deallocate();
    T15.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(26, total_steps, steps_left, start_index, 15, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M16 = (1.0 * A21 + -1.0 * A32 + 1.0 * A41) * (1.0 * B21 + 1.0 * B23 + -1.0 * B_X2 + 1.0 * B_X4)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(26, total_steps, steps_left, start_index, 16, num_threads)) untied
    {
#endif
    Matrix<Scalar> S16(A11.m(), A11.n());
    S_Add16(A21, A32, A41, S16, x);
    Matrix<Scalar> T16(B11.m(), B11.n());
    T_Add16(B21, B23, B_X2, B_X4, T16, x);
    FastMatmulRecursive(S16, T16, M16, total_steps, steps_left - 1, (start_index + 16 - 1) * 26, x, num_threads, Scalar(0.0));
    S16.deallocate();
    T16.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(26, total_steps, steps_left, start_index, 16, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M17 = (1.0 * A_X3) * (-1.0 * B_X2 + 1.0 * B_X4 + 1.0 * B_X6)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(26, total_steps, steps_left, start_index, 17, num_threads)) untied
    {
#endif
    Matrix<Scalar> T17(B11.m(), B11.n());
    T_Add17(B_X2, B_X4, B_X6, T17, x);
    FastMatmulRecursive(A_X3, T17, M17, total_steps, steps_left - 1, (start_index + 17 - 1) * 26, x, num_threads, Scalar(0.0));
    T17.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(26, total_steps, steps_left, start_index, 17, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M18 = (1.0 * A21 + -1.0 * A31 + 1.0 * A41) * (-1.0 * B22 + 1.0 * B23 + 1.0 * B_X3)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(26, total_steps, steps_left, start_index, 18, num_threads)) untied
    {
#endif
    Matrix<Scalar> S18(A11.m(), A11.n());
    S_Add18(A21, A31, A41, S18, x);
    Matrix<Scalar> T18(B11.m(), B11.n());
    T_Add18(B22, B23, B_X3, T18, x);
    FastMatmulRecursive(S18, T18, M18, total_steps, steps_left - 1, (start_index + 18 - 1) * 26, x, num_threads, Scalar(0.0));
    S18.deallocate();
    T18.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(26, total_steps, steps_left, start_index, 18, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M19 = (1.0 * A21 + -1.0 * A31 + -1.0 * A41) * (-1.0 * B13)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(26, total_steps, steps_left, start_index, 19, num_threads)) untied
    {
#endif
    Matrix<Scalar> S19(A11.m(), A11.n());
    S_Add19(A21, A31, A41, S19, x);
    M19.UpdateMultiplier(Scalar(-1));
    FastMatmulRecursive(S19, B13, M19, total_steps, steps_left - 1, (start_index + 19 - 1) * 26, x, num_threads, Scalar(0.0));
    S19.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(26, total_steps, steps_left, start_index, 19, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M20 = (1.0 * A_X5) * (1.0 * B12 + 1.0 * B22 + -1.0 * B_X2 + -1.0 * B_X6)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(26, total_steps, steps_left, start_index, 20, num_threads)) untied
    {
#endif
    Matrix<Scalar> T20(B11.m(), B11.n());
    T_Add20(B12, B22, B_X2, B_X6, T20, x);
    FastMatmulRecursive(A_X5, T20, M20, total_steps, steps_left - 1, (start_index + 20 - 1) * 26, x, num_threads, Scalar(0.0));
    T20.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(26, total_steps, steps_left, start_index, 20, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M21 = (-1.0 * A21 + 1.0 * A_X2 + 1.0 * A_X5 + 1.0 * A_X6) * (1.0 * B23 + 1.0 * B_X3 + 1.0 * B_X5)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(26, total_steps, steps_left, start_index, 21, num_threads)) untied
    {
#endif
    Matrix<Scalar> S21(A11.m(), A11.n());
    S_Add21(A21, A_X2, A_X5, A_X6, S21, x);
    Matrix<Scalar> T21(B11.m(), B11.n());
    T_Add21(B23, B_X3, B_X5, T21, x);
    FastMatmulRecursive(S21, T21, M21, total_steps, steps_left - 1, (start_index + 21 - 1) * 26, x, num_threads, Scalar(0.0));
    S21.deallocate();
    T21.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(26, total_steps, steps_left, start_index, 21, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M22 = (1.0 * A22) * (-1.0 * B_X2 + -1.0 * B_X6)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(26, total_steps, steps_left, start_index, 22, num_threads)) untied
    {
#endif
    Matrix<Scalar> T22(B11.m(), B11.n());
    T_Add22(B_X2, B_X6, T22, x);
    FastMatmulRecursive(A22, T22, M22, total_steps, steps_left - 1, (start_index + 22 - 1) * 26, x, num_threads, Scalar(0.0));
    T22.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(26, total_steps, steps_left, start_index, 22, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M23 = (1.0 * A22 + -1.0 * A32 + -1.0 * A42) * (-1.0 * B23)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(26, total_steps, steps_left, start_index, 23, num_threads)) untied
    {
#endif
    Matrix<Scalar> S23(A11.m(), A11.n());
    S_Add23(A22, A32, A42, S23, x);
    M23.UpdateMultiplier(Scalar(-1));
    FastMatmulRecursive(S23, B23, M23, total_steps, steps_left - 1, (start_index + 23 - 1) * 26, x, num_threads, Scalar(0.0));
    S23.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(26, total_steps, steps_left, start_index, 23, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M24 = (1.0 * A21 + 1.0 * A22 + -1.0 * A_X5 + -1.0 * A_X6) * (1.0 * B22 + -1.0 * B23)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(26, total_steps, steps_left, start_index, 24, num_threads)) untied
    {
#endif
    Matrix<Scalar> S24(A11.m(), A11.n());
    S_Add24(A21, A22, A_X5, A_X6, S24, x);
    Matrix<Scalar> T24(B11.m(), B11.n());
    T_Add24(B22, B23, T24, x);
    FastMatmulRecursive(S24, T24, M24, total_steps, steps_left - 1, (start_index + 24 - 1) * 26, x, num_threads, Scalar(0.0));
    S24.deallocate();
    T24.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(26, total_steps, steps_left, start_index, 24, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M25 = (-0.5 * A21 + 0.5 * A31 + -0.5 * A41 + -1.0 * A42) * (1.0 * B23 + 1.0 * B_X4)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(26, total_steps, steps_left, start_index, 25, num_threads)) untied
    {
#endif
    Matrix<Scalar> S25(A11.m(), A11.n());
    S_Add25(A21, A31, A41, A42, S25, x);
    Matrix<Scalar> T25(B11.m(), B11.n());
    T_Add25(B23, B_X4, T25, x);
    FastMatmulRecursive(S25, T25, M25, total_steps, steps_left - 1, (start_index + 25 - 1) * 26, x, num_threads, Scalar(0.0));
    S25.deallocate();
    T25.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(26, total_steps, steps_left, start_index, 25, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M26 = (1.0 * A21 + -1.0 * A22) * (-1.0 * B_X2)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(26, total_steps, steps_left, start_index, 26, num_threads)) untied
    {
#endif
    Matrix<Scalar> S26(A11.m(), A11.n());
    S_Add26(A21, A22, S26, x);
    M26.UpdateMultiplier(Scalar(-1));
    FastMatmulRecursive(S26, B_X2, M26, total_steps, steps_left - 1, (start_index + 26 - 1) * 26, x, num_threads, Scalar(0.0));
    S26.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(26, total_steps, steps_left, start_index, 26, num_threads)) {
# pragma omp taskwait
    }
#endif

    Matrix<Scalar> M_X1(M11.m(), M11.n());
    MX_Add1(M6, M11, M_X1);
    Matrix<Scalar> M_X2(M11.m(), M11.n());
    MX_Add2(M12, M13, M_X2);
    Matrix<Scalar> M_X3(M11.m(), M11.n());
    MX_Add3(M2, M8, M_X3);
    Matrix<Scalar> M_X4(M11.m(), M11.n());
    MX_Add4(M19, M26, M_X4);
    Matrix<Scalar> M_X5(M11.m(), M11.n());
    MX_Add5(M2, M16, M_X5);
    Matrix<Scalar> M_X6(M11.m(), M11.n());
    MX_Add6(M20, M22, M_X6);
    Matrix<Scalar> M_X7(M11.m(), M11.n());
    MX_Add7(M4, M7, M_X7);
    Matrix<Scalar> M_X8(M11.m(), M11.n());
    MX_Add8(M5, M11, M_X8);
    Matrix<Scalar> M_X9(M11.m(), M11.n());
    MX_Add9(M10, M18, M_X9);
    Matrix<Scalar> M_X10(M11.m(), M11.n());
    MX_Add10(M14, M15, M_X10);
    Matrix<Scalar> M_X11(M11.m(), M11.n());
    MX_Add11(M16, M17, M_X11);
    Matrix<Scalar> M_X12(M11.m(), M11.n());
    MX_Add12(M5, M23, M_X12);
    Matrix<Scalar> M_X13(M11.m(), M11.n());
    MX_Add13(M18, M25, M_X13);
    M_Add1(M1, M4, M9, M14, C11, x, beta);
    M_Add2(M_X1, M_X2, C12, x, beta);
    M_Add3(M3, M14, M21, M24, M_X1, M_X2, C13, x, beta);
    M_Add4(M12, M14, C14, x, beta);
    M_Add5(M_X3, M_X4, C21, x, beta);
    M_Add6(M13, M24, M25, M_X5, M_X6, C22, x, beta);
    M_Add7(M19, M22, M_X3, C23, x, beta);
    M_Add8(M8, M12, M21, M_X7, M_X8, M_X9, M_X10, C24, x, beta);
    M_Add9(M8, M18, M_X11, C31, x, beta);
    M_Add10(M25, M_X5, M_X6, M_X12, C32, x, beta);
    M_Add11(M22, M_X3, M_X12, M_X13, C33, x, beta);
    M_Add12(M8, M21, M_X7, M_X9, C34, x, beta);
    M_Add13(M2, M9, M18, M_X4, M_X11, C41, x, beta);
    M_Add14(M5, M13, C42, x, beta);
    M_Add15(M5, M19, M_X13, C43, x, beta);
    M_Add16(M_X8, M_X10, C44, x, beta);

    // Handle edge cases with dynamic peeling
    DynamicPeeling(A, B, C, 4, 2, 4, beta);
}

// C := alpha * A * B + beta * C
template <typename Scalar>
void FastMatmul(Matrix<Scalar>& A, Matrix<Scalar>& B, Matrix<Scalar>& C,
    int num_steps, double x=1e-8, Scalar alpha=Scalar(1.0), Scalar beta=Scalar(0.0)) {
    A.set_multiplier(alpha);
    int num_multiplies_per_step = 26;
    int total_multiplies = pow(num_multiplies_per_step, num_steps);
#ifdef _PARALLEL_
    int num_threads = -1;
# pragma omp parallel
    {
        if (omp_get_thread_num() == 0) { num_threads = omp_get_num_threads(); }
    }
#else
    int num_threads = 0;
#endif
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
    omp_set_nested(1);
    if (num_threads > total_multiplies) {
        mkl_set_dynamic(0);
    }
# pragma omp parallel
    {
# pragma omp single
#endif
        FastMatmulRecursive(A, B, C, num_steps, num_steps, 0, x, num_threads, beta);
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ || _PARALLEL_ == _HYBRID_PAR_)
    }
#endif
}

}  // namespace grey424_26_206

#endif  // _grey424_26_206_HPP_
