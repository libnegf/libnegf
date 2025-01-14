/*!!--------------------------------------------------------------------------!
 *!! libNEGF: a general library for Non-Equilibrium Greens functions.         !
 *!! Copyright (C) 2025                                                       !
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

#include <cassert>
#include <climits>

#include <cusolverDn.h>


namespace libnegf {

// Ensure cusolverDngetrf* is in a namespace to avoid clashing with
// the deprecated cuSolver function of the same name.

cusolverStatus_t cusolverDngetrf(
    cusolverDnHandle_t handle, size_t m, size_t n, cuComplex* A, size_t lda, cuComplex* workspace, int *d_ipiv, int* d_info)
{
    assert(m <= INT_MAX);
    assert(n <= INT_MAX);
    assert(A);
    assert(lda <= INT_MAX);
	assert(workspace);
	assert(d_ipiv);
	assert(d_info);

    return cusolverDnCgetrf(handle, m, n, A, lda, workspace, d_ipiv, d_info);
}

cusolverStatus_t cusolverDngetrf(
    cusolverDnHandle_t handle, size_t m, size_t n, cuDoubleComplex* A, size_t lda, cuDoubleComplex* workspace, int *d_ipiv, int* d_info)
{
    assert(m <= INT_MAX);
    assert(n <= INT_MAX);
    assert(A);
    assert(lda <= INT_MAX);
	assert(workspace);
	assert(d_ipiv);
	assert(d_info);

    return cusolverDnZgetrf(handle, m, n, A, lda, workspace, d_ipiv, d_info);
}


cusolverStatus_t cusolverDngetrs(
    cusolverDnHandle_t handle, cublasOperation_t trans, size_t n, size_t nrhs, const cuComplex* A, size_t lda, const int* d_ipiv, cuComplex* B, size_t ldb, int* d_info)
{
	assert(n <= INT_MAX);
	assert(nrhs <= INT_MAX);
	assert(A);
	assert(d_ipiv);
	assert(B);
	assert(d_info);

	return cusolverDnCgetrs(handle, trans, n, nrhs, A, lda, d_ipiv, B, ldb, d_info);
}

cusolverStatus_t cusolverDngetrs(
    cusolverDnHandle_t handle, cublasOperation_t trans, size_t n, size_t nrhs, const cuDoubleComplex* A, size_t lda, const int* d_ipiv, cuDoubleComplex* B, size_t ldb, int* d_info)
{
	assert(n <= INT_MAX);
	assert(nrhs <= INT_MAX);
	assert(A);
	assert(d_ipiv);
	assert(B);
	assert(d_info);

	return cusolverDnZgetrs(handle, trans, n, nrhs, A, lda, d_ipiv, B, ldb, d_info);
}


cusolverStatus_t cusolverDngetrf_bufferSize(
    cusolverDnHandle_t handle, size_t m, size_t n, cuComplex* A, size_t lda,
    int* lwork
) {
    assert(m <= INT_MAX);
    assert(n <= INT_MAX);
    assert(A);
    assert(lda <= INT_MAX);
    assert(lwork);

    return cusolverDnCgetrf_bufferSize(handle, m, n, A, lda, lwork);
}

cusolverStatus_t cusolverDngetrf_bufferSize(
    cusolverDnHandle_t handle, size_t m, size_t n, cuDoubleComplex* A,
    size_t lda, int* lwork
) {
    assert(m <= INT_MAX);
    assert(n <= INT_MAX);
    assert(A);
    assert(lda <= INT_MAX);
    assert(lwork);

    return cusolverDnZgetrf_bufferSize(handle, m, n, A, lda, lwork);
}

} // namespace libnegf
