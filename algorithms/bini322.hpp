#ifndef _bini322_HPP_
#define _bini322_HPP_

// This is an automatically generated file from gen.py.
#include "common.hpp"

namespace bini322 {

template <typename Scalar>
void S_Add1(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x=1e-8) {
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
            dataC[i + j * strideC] = dataS1[i + j * strideS1] + dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add2(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x=1e-8) {
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
void S_Add3(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x=1e-8) {
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
void S_Add4(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x=1e-8) {
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
            dataC[i + j * strideC] = Scalar(x) * dataS1[i + j * strideS1] + dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add5(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x=1e-8) {
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
            dataC[i + j * strideC] = dataS1[i + j * strideS1] + Scalar(x) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add6(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x=1e-8) {
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
            dataC[i + j * strideC] = dataS1[i + j * strideS1] + dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add7(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x=1e-8) {
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
void S_Add8(Matrix<Scalar>& S1, Matrix<Scalar>& C, double x=1e-8) {
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
void S_Add9(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x=1e-8) {
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
            dataC[i + j * strideC] = Scalar(x) * dataS1[i + j * strideS1] + dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void S_Add10(Matrix<Scalar>& S1, Matrix<Scalar>& S2, Matrix<Scalar>& C, double x=1e-8) {
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
            dataC[i + j * strideC] = dataS1[i + j * strideS1] + Scalar(x) * dataS2[i + j * strideS2];
        }
    }
}

template <typename Scalar>
void T_Add1(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x=1e-8) {
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
            dataC[i + j * strideC] = Scalar(x) * dataT1[i + j * strideT1] + dataT2[i + j * strideT2];
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
            dataC[i + j * strideC] = -dataT1[i + j * strideT1] -dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add3(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x=1e-8) {
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
            dataC[i + j * strideC] = Scalar(-(x)) * dataT1[i + j * strideT1] + dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add5(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x=1e-8) {
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
            dataC[i + j * strideC] = Scalar(x) * dataT1[i + j * strideT1] + dataT2[i + j * strideT2];
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
            dataC[i + j * strideC] = dataT1[i + j * strideT1] + Scalar(x) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add7(Matrix<Scalar>& T1, Matrix<Scalar>& C, double x=1e-8) {
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
            dataC[i + j * strideC] = -dataT1[i + j * strideT1] -dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add9(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x=1e-8) {
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
            dataC[i + j * strideC] = dataT1[i + j * strideT1] + Scalar(x) * dataT2[i + j * strideT2];
        }
    }
}

template <typename Scalar>
void T_Add10(Matrix<Scalar>& T1, Matrix<Scalar>& T2, Matrix<Scalar>& C, double x=1e-8) {
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
            dataC[i + j * strideC] = dataT1[i + j * strideT1] + Scalar(-(x)) * dataT2[i + j * strideT2];
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
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(-(1.0 / (x))) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + Scalar(-(1.0 / (x))) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4];
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
                dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2];
            }
        }
    }
}

template <typename Scalar>
void M_Add3(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& C, double x=1e-8, Scalar beta=Scalar(0.0)) {
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
                dataC[i + j * strideC] = dataM1[i + j * strideM1] + dataM2[i + j * strideM2] -dataM3[i + j * strideM3] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = dataM1[i + j * strideM1] + dataM2[i + j * strideM2] -dataM3[i + j * strideM3];
            }
        }
    }
}

template <typename Scalar>
void M_Add4(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& C, double x=1e-8, Scalar beta=Scalar(0.0)) {
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
                dataC[i + j * strideC] = dataM1[i + j * strideM1] -dataM2[i + j * strideM2] + dataM3[i + j * strideM3] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = dataM1[i + j * strideM1] -dataM2[i + j * strideM2] + dataM3[i + j * strideM3];
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
                dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(-(1.0 / (x))) * dataM1[i + j * strideM1] + Scalar(1.0 / (x)) * dataM2[i + j * strideM2];
            }
        }
    }
}

template <typename Scalar>
void M_Add6(Matrix<Scalar>& M1, Matrix<Scalar>& M2, Matrix<Scalar>& M3, Matrix<Scalar>& M4, Matrix<Scalar>& C, double x=1e-8, Scalar beta=Scalar(0.0)) {
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
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(-(1.0 / (x))) * dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4] + beta * dataC[i + j * strideC];
            }
        }
    } else {
#ifdef _PARALLEL_
# pragma omp parallel for collapse(2)
#endif
        for (int j = 0; j < C.n(); ++j) {
            for (int i = 0; i < C.m(); ++i) {
                dataC[i + j * strideC] = Scalar(1.0 / (x)) * dataM1[i + j * strideM1] + Scalar(-(1.0 / (x))) * dataM2[i + j * strideM2] + Scalar(1.0 / (x)) * dataM3[i + j * strideM3] + Scalar(1.0 / (x)) * dataM4[i + j * strideM4];
            }
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

    Matrix<Scalar> A11 = A.Subblock(3, 2, 1, 1);
    Matrix<Scalar> A12 = A.Subblock(3, 2, 1, 2);
    Matrix<Scalar> A21 = A.Subblock(3, 2, 2, 1);
    Matrix<Scalar> A22 = A.Subblock(3, 2, 2, 2);
    Matrix<Scalar> A31 = A.Subblock(3, 2, 3, 1);
    Matrix<Scalar> A32 = A.Subblock(3, 2, 3, 2);
    Matrix<Scalar> B11 = B.Subblock(2, 2, 1, 1);
    Matrix<Scalar> B12 = B.Subblock(2, 2, 1, 2);
    Matrix<Scalar> B21 = B.Subblock(2, 2, 2, 1);
    Matrix<Scalar> B22 = B.Subblock(2, 2, 2, 2);
    Matrix<Scalar> C11 = C.Subblock(3, 2, 1, 1);
    Matrix<Scalar> C12 = C.Subblock(3, 2, 1, 2);
    Matrix<Scalar> C21 = C.Subblock(3, 2, 2, 1);
    Matrix<Scalar> C22 = C.Subblock(3, 2, 2, 2);
    Matrix<Scalar> C31 = C.Subblock(3, 2, 3, 1);
    Matrix<Scalar> C32 = C.Subblock(3, 2, 3, 2);


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



    // M1 = (1.0 * A11 + 1.0 * A22) * (x * B11 + 1.0 * B22)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(10, total_steps, steps_left, start_index, 1, num_threads)) untied
    {
#endif
    Matrix<Scalar> S1(A11.m(), A11.n());
    S_Add1(A11, A22, S1, x);
    Matrix<Scalar> T1(B11.m(), B11.n());
    T_Add1(B11, B22, T1, x);
    FastMatmulRecursive(S1, T1, M1, total_steps, steps_left - 1, (start_index + 1 - 1) * 10, x, num_threads, Scalar(0.0));
    S1.deallocate();
    T1.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(10, total_steps, steps_left, start_index, 1, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M2 = (1.0 * A22) * (-1.0 * B21 + -1.0 * B22)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(10, total_steps, steps_left, start_index, 2, num_threads)) untied
    {
#endif
    Matrix<Scalar> T2(B11.m(), B11.n());
    T_Add2(B21, B22, T2, x);
    FastMatmulRecursive(A22, T2, M2, total_steps, steps_left - 1, (start_index + 2 - 1) * 10, x, num_threads, Scalar(0.0));
    T2.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(10, total_steps, steps_left, start_index, 2, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M3 = (1.0 * A11) * (1.0 * B22)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(10, total_steps, steps_left, start_index, 3, num_threads)) untied
    {
#endif
    FastMatmulRecursive(A11, B22, M3, total_steps, steps_left - 1, (start_index + 3 - 1) * 10, x, num_threads, Scalar(0.0));
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(10, total_steps, steps_left, start_index, 3, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M4 = (x * A12 + 1.0 * A22) * (-(x) * B11 + 1.0 * B21)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(10, total_steps, steps_left, start_index, 4, num_threads)) untied
    {
#endif
    Matrix<Scalar> S4(A11.m(), A11.n());
    S_Add4(A12, A22, S4, x);
    Matrix<Scalar> T4(B11.m(), B11.n());
    T_Add4(B11, B21, T4, x);
    FastMatmulRecursive(S4, T4, M4, total_steps, steps_left - 1, (start_index + 4 - 1) * 10, x, num_threads, Scalar(0.0));
    S4.deallocate();
    T4.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(10, total_steps, steps_left, start_index, 4, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M5 = (1.0 * A11 + x * A12) * (x * B12 + 1.0 * B22)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(10, total_steps, steps_left, start_index, 5, num_threads)) untied
    {
#endif
    Matrix<Scalar> S5(A11.m(), A11.n());
    S_Add5(A11, A12, S5, x);
    Matrix<Scalar> T5(B11.m(), B11.n());
    T_Add5(B12, B22, T5, x);
    FastMatmulRecursive(S5, T5, M5, total_steps, steps_left - 1, (start_index + 5 - 1) * 10, x, num_threads, Scalar(0.0));
    S5.deallocate();
    T5.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(10, total_steps, steps_left, start_index, 5, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M6 = (1.0 * A21 + 1.0 * A32) * (1.0 * B11 + x * B22)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(10, total_steps, steps_left, start_index, 6, num_threads)) untied
    {
#endif
    Matrix<Scalar> S6(A11.m(), A11.n());
    S_Add6(A21, A32, S6, x);
    Matrix<Scalar> T6(B11.m(), B11.n());
    T_Add6(B11, B22, T6, x);
    FastMatmulRecursive(S6, T6, M6, total_steps, steps_left - 1, (start_index + 6 - 1) * 10, x, num_threads, Scalar(0.0));
    S6.deallocate();
    T6.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(10, total_steps, steps_left, start_index, 6, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M7 = (1.0 * A32) * (1.0 * B11)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(10, total_steps, steps_left, start_index, 7, num_threads)) untied
    {
#endif
    FastMatmulRecursive(A32, B11, M7, total_steps, steps_left - 1, (start_index + 7 - 1) * 10, x, num_threads, Scalar(0.0));
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(10, total_steps, steps_left, start_index, 7, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M8 = (1.0 * A21) * (-1.0 * B11 + -1.0 * B12)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(10, total_steps, steps_left, start_index, 8, num_threads)) untied
    {
#endif
    Matrix<Scalar> T8(B11.m(), B11.n());
    T_Add8(B11, B12, T8, x);
    FastMatmulRecursive(A21, T8, M8, total_steps, steps_left - 1, (start_index + 8 - 1) * 10, x, num_threads, Scalar(0.0));
    T8.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(10, total_steps, steps_left, start_index, 8, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M9 = (x * A31 + 1.0 * A32) * (1.0 * B11 + x * B21)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(10, total_steps, steps_left, start_index, 9, num_threads)) untied
    {
#endif
    Matrix<Scalar> S9(A11.m(), A11.n());
    S_Add9(A31, A32, S9, x);
    Matrix<Scalar> T9(B11.m(), B11.n());
    T_Add9(B11, B21, T9, x);
    FastMatmulRecursive(S9, T9, M9, total_steps, steps_left - 1, (start_index + 9 - 1) * 10, x, num_threads, Scalar(0.0));
    S9.deallocate();
    T9.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(10, total_steps, steps_left, start_index, 9, num_threads)) {
# pragma omp taskwait
# if defined(_PARALLEL_) && (_PARALLEL_ == _HYBRID_PAR_)
    mkl_set_num_threads_local(num_threads);
# endif
    }
#endif

    // M10 = (1.0 * A21 + x * A31) * (1.0 * B12 + -(x) * B22)
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
# pragma omp task if(should_launch_task(10, total_steps, steps_left, start_index, 10, num_threads)) untied
    {
#endif
    Matrix<Scalar> S10(A11.m(), A11.n());
    S_Add10(A21, A31, S10, x);
    Matrix<Scalar> T10(B11.m(), B11.n());
    T_Add10(B12, B22, T10, x);
    FastMatmulRecursive(S10, T10, M10, total_steps, steps_left - 1, (start_index + 10 - 1) * 10, x, num_threads, Scalar(0.0));
    S10.deallocate();
    T10.deallocate();
#if defined(_PARALLEL_) && (_PARALLEL_ == _BFS_PAR_ ||  _PARALLEL_ == _HYBRID_PAR_)
    }
    if (should_task_wait(10, total_steps, steps_left, start_index, 10, num_threads)) {
# pragma omp taskwait
    }
#endif

    M_Add1(M1, M2, M3, M4, C11, x, beta);
    M_Add2(M3, M5, C12, x, beta);
    M_Add3(M4, M6, M9, C21, x, beta);
    M_Add4(M1, M5, M10, C22, x, beta);
    M_Add5(M7, M9, C31, x, beta);
    M_Add6(M6, M7, M8, M10, C32, x, beta);

    // Handle edge cases with dynamic peeling
    DynamicPeeling(A, B, C, 3, 2, 2, beta);
}

// C := alpha * A * B + beta * C
template <typename Scalar>
void FastMatmul(Matrix<Scalar>& A, Matrix<Scalar>& B, Matrix<Scalar>& C,
    int num_steps, double x=1e-8, Scalar alpha=Scalar(1.0), Scalar beta=Scalar(0.0)) {
    A.set_multiplier(alpha);
    int num_multiplies_per_step = 10;
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

}  // namespace bini322

#endif  // _bini322_HPP_
