/*!!--------------------------------------------------------------------------!
 *!! libNEGF: a general library for Non-Equilibrium Greens functions.         !
 *!! Copyright (C) 2024                                                       !
 *!!                                                                          !
 *!! This file is part of libNEGF: a library for                              !
 *!! Non Equilibrium Green's Function calculation                             !
 *!!                                                                          !
 *!! Developers: Alessandro Pecchia, Daniele Soccodato                        !
 *!! Former Contributors: Gabriele Penazzi, Luca Latessa, Aldo Di Carlo       !
 *!!                                                                          !
 *!! libNEGF is free software: you can redistribute it and/or modify          !
 *!! it under the terms of the GNU Lesse General Public License as published  !
 *!! by the Free Software Foundation, either version 3 of the License, or     !
 *!! (at your option) any later version.                                      !
 *!!                                                                          !
 *!!  You should have received a copy of the GNU Lesser General Public        !
 *!!  License along with libNEGF.  If not, see                                !
 *!!  <http://www.gnu.org/licenses/>.                                         !
 *!!--------------------------------------------------------------------------!
 */

#undef NDEBUG
#include <cassert>
#include <climits>

#include <cublas_v2.h>


/// This files contains C++ overloads of cuBLAS and cuSolver functions.


namespace libnegf {

cublasStatus_t cublasAsum(
    cublasHandle_t handle, size_t n, const cuComplex* x, size_t incx,
    float* result
) {
    assert(n <= INT_MAX);

    return cublasScasum(handle, n, x, incx, result);
}

cublasStatus_t cublasAsum(
    cublasHandle_t handle, size_t n, const cuDoubleComplex* x, size_t incx,
    double* result
) {
    assert(n <= INT_MAX);

    return cublasDzasum(handle, n, x, incx, result);
}


cublasStatus_t cublasAxpy(
    cublasHandle_t handle, size_t n, const cuComplex* alpha, const cuComplex* x,
    size_t incx, cuComplex* y, size_t incy
) {
    assert(n <= INT_MAX);
    assert(alpha);
    assert(x);
    assert(y);

    return cublasCaxpy(handle, n, alpha, x, incx, y, incy);
}

cublasStatus_t cublasAxpy(
    cublasHandle_t handle, size_t n, const cuDoubleComplex* alpha,
    const cuDoubleComplex* x, size_t incx, cuDoubleComplex* y, size_t incy
) {
    assert(n <= INT_MAX);
    assert(alpha);
    assert(x);
    assert(y);

    return cublasZaxpy(handle, n, alpha, x, incx, y, incy);
}


cublasStatus_t cublasCopy(
    cublasHandle_t handle, size_t n, const cuComplex* d_x, size_t incx,
    cuComplex* d_y, size_t incy
) {
    assert(n <= INT_MAX);
    assert(d_x);
    assert(d_y);

    return cublasCcopy(handle, n, d_x, incx, d_y, incy);
}

cublasStatus_t cublasCopy(
    cublasHandle_t handle, size_t n, const cuDoubleComplex* d_x, size_t incx,
    cuDoubleComplex* d_y, size_t incy
) {
    assert(n <= INT_MAX);
    assert(d_x);
    assert(d_y);

    return cublasZcopy(handle, n, d_x, incx, d_y, incy);
}


// cublasXdot

cublasStatus_t cublasDot(
    cublasHandle_t handle, size_t n, const float* x, size_t incx,
    const float* y, size_t incy, float* result
) {
    assert(n);
    assert(x);
    assert(y);
    assert(result);

    return cublasSdot(handle, n, x, incx, y, incy, result);
}

cublasStatus_t cublasDot(
    cublasHandle_t handle, size_t n, const double* x, size_t incx,
    const double* y, size_t incy, double* result
) {
    assert(n);
    assert(x);
    assert(y);
    assert(result);

    return cublasDdot(handle, n, x, incx, y, incy, result);
}


// cublasXgeam

cublasStatus_t cublasGeam(
    cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb,
    size_t m, size_t n, const float* alpha, const float* A, size_t lda,
    const float* beta, const float* B, size_t ldb, float* C, size_t ldc
) {
    assert(m <= INT_MAX);
    assert(n <= INT_MAX);
    assert(lda <= INT_MAX);
    assert(ldb <= INT_MAX);
    assert(ldc <= INT_MAX);
    return cublasSgeam(
        handle, transa, transb, m, n, alpha, A, lda, beta, B, ldb, C, ldc
    );
}

cublasStatus_t cublasGeam(
    cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb,
    size_t m, size_t n, const double* alpha, const double* A, size_t lda,
    const double* beta, const double* B, size_t ldb, double* C, size_t ldc
) {
    assert(m <= INT_MAX);
    assert(n <= INT_MAX);
    assert(lda <= INT_MAX);
    assert(ldb <= INT_MAX);
    assert(ldc <= INT_MAX);
    return cublasDgeam(
        handle, transa, transb, m, n, alpha, A, lda, beta, B, ldb, C, ldc
    );
}

cublasStatus_t cublasGeam(
    cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb,
    size_t m, size_t n, const cuComplex* alpha, const cuComplex* A, size_t lda,
    const cuComplex* beta, const cuComplex* B, size_t ldb, cuComplex* C,
    size_t ldc
) {
    assert(m <= INT_MAX);
    assert(n <= INT_MAX);
    assert(lda <= INT_MAX);
    assert(ldb <= INT_MAX);
    assert(ldc <= INT_MAX);
    return cublasCgeam(
        handle, transa, transb, m, n, alpha, A, lda, beta, B, ldb, C, ldc
    );
}

cublasStatus_t cublasGeam(
    cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb,
    size_t m, size_t n, const cuDoubleComplex* alpha, const cuDoubleComplex* A,
    size_t lda, const cuDoubleComplex* beta, const cuDoubleComplex* B,
    size_t ldb, cuDoubleComplex* C, size_t ldc
) {
    assert(m <= INT_MAX);
    assert(n <= INT_MAX);
    assert(lda <= INT_MAX);
    assert(ldb <= INT_MAX);
    assert(ldc <= INT_MAX);
    return cublasZgeam(
        handle, transa, transb, m, n, alpha, A, lda, beta, B, ldb, C, ldc
    );
}


// cublasXgemm

cublasStatus_t cublasGemm(
    cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb,
    size_t m, size_t n, size_t k, const float* alpha, const float* A,
    size_t lda, const float* B, size_t ldb, const float* beta, float* C,
    size_t ldc
) {
    assert(m <= INT_MAX);
    assert(n <= INT_MAX);
    assert(k <= INT_MAX);
    assert(lda <= INT_MAX);
    assert(ldb <= INT_MAX);
    assert(ldc <= INT_MAX);
    return cublasSgemm(
        handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc
    );
}

cublasStatus_t cublasGemm(
    cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb,
    size_t m, size_t n, size_t k, const double* alpha, const double* A,
    size_t lda, const double* B, size_t ldb, const double* beta, double* C,
    size_t ldc
) {
    assert(m <= INT_MAX);
    assert(n <= INT_MAX);
    assert(k <= INT_MAX);
    assert(lda <= INT_MAX);
    assert(ldb <= INT_MAX);
    assert(ldc <= INT_MAX);
    return cublasDgemm(
        handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc
    );
}

cublasStatus_t cublasGemm(
    cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb,
    size_t m, size_t n, size_t k, const cuComplex* alpha, const cuComplex* A,
    size_t lda, const cuComplex* B, size_t ldb, const cuComplex* beta,
    cuComplex* C, size_t ldc
) {
    assert(m <= INT_MAX);
    assert(n <= INT_MAX);
    assert(k <= INT_MAX);
    assert(lda <= INT_MAX);
    assert(ldb <= INT_MAX);
    assert(ldc <= INT_MAX);
    return cublasCgemm(
        handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc
    );
}

cublasStatus_t cublasGemm(
    cublasHandle_t handle, cublasOperation_t transa, cublasOperation_t transb,
    size_t m, size_t n, size_t k, const cuDoubleComplex* alpha,
    const cuDoubleComplex* A, size_t lda, const cuDoubleComplex* B, size_t ldb,
    const cuDoubleComplex* beta, cuDoubleComplex* C, size_t ldc
) {
    assert(m <= INT_MAX);
    assert(n <= INT_MAX);
    assert(k <= INT_MAX);
    assert(lda <= INT_MAX);
    assert(ldb <= INT_MAX);
    assert(ldc <= INT_MAX);
    return cublasZgemm(
        handle, transa, transb, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc
    );
}

} // namespace libnegf
