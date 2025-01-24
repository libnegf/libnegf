/*!!--------------------------------------------------------------------------!
 *!! libNEGF: a general library for Non-Equilibrium Greens functions.         !
 *!! Copyright (C) 2012 - 2026                                                !
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

#include <cuda.h>
#include <cuda_runtime.h>

#include <cassert>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include <type_traits>

#include <libnegf/cublas_v2.h>
#include <libnegf/cusolverDn.h>


#if __CUDA_ARCH__
#define ENFORCE(expr) \
    do { \
        if(!(expr)) { \
            printf( \
                "%s:%u: ENFORCEMENT FAILURE: %s is FALSE\n", __FILE__, \
                __LINE__, #expr \
            ); \
            __trap(); \
        } \
    } while(false)
#else
#define ENFORCE(expr) \
    do { \
        if(!(expr)) { \
            std::fprintf( \
                stderr, "%s:%u: ENFORCEMENT FAILURE: %s is FALSE\n", __FILE__, \
                __LINE__, #expr \
            ); \
            std::exit(EXIT_FAILURE); \
        } \
    } while(false)
#endif


constexpr auto BLOCK_SIZE = std::size_t{1024};


/**
 * Given a floating-point type, returns the associated real-valued type; for
 * real-valued types, this is the type itself.
 */
template<typename>
struct get_real {};

template<>
struct get_real<float> {
    using type = float;
};

template<>
struct get_real<double> {
    using type = double;
};

template<>
struct get_real<cuComplex> {
    using type = float;
};

template<>
struct get_real<cuDoubleComplex> {
    using type = double;
};


/**
 * Computes c = α a · β b.
 */
template<typename Number>
__global__ void addKernel(
    Number* c, Number alpha, const Number* a, Number beta, const Number* b,
    size_t size
) {
    ENFORCE(c);
    ENFORCE(a);
    ENFORCE(b);

    auto i = blockIdx.x * blockDim.x + threadIdx.x;

    if(i < size) {
        c[i].x = (alpha.x * a[i].x - alpha.y * a[i].y) +
                 (beta.x * b[i].x - beta.y * b[i].y);
        c[i].y = (alpha.x * a[i].y + alpha.y * a[i].x) +
                 (beta.x * b[i].y + beta.y * b[i].x);
    }
}


/*
__global__ void hermitian(cuComplex *odata, const cuComplex *idata)
{
  __shared__ cuComplex tile[TILE_DIM][TILE_DIM];

  int x = blockIdx.x * TILE_DIM + threadIdx.x;
  int y = blockIdx.y * TILE_DIM + threadIdx.y;
  int width = gridDim.x * TILE_DIM;

  for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
     tile[threadIdx.y+j][threadIdx.x] = idata[(y+j)*width + x];

  __syncthreads();

  x = blockIdx.y * TILE_DIM + threadIdx.x;  // transpose block offset
  y = blockIdx.x * TILE_DIM + threadIdx.y;

  for (int j = 0; j < TILE_DIM; j += BLOCK_ROWS)
  {
     odata[(y+j)*width + x].x = tile[threadIdx.x][threadIdx.y + j].x;
     odata[(y+j)*width + x].y = -tile[threadIdx.x][threadIdx.y + j].y;
  }
}
*/


/**
 * Initializes a square complex matrix as the identity matrix.
 */
template<typename Number>
__global__ void initKernel(Number* a, size_t nrow) {
    using Real = typename get_real<Number>::type;

    ENFORCE(a);

    auto size = nrow * nrow;
    auto i = blockDim.x * blockIdx.x + threadIdx.x;
    auto one = Number{Real{1}};
    auto zero = Number{Real{0}};

    if(i < size) {
        if(i % (nrow + 1) == 0) {
            a[i] = one;
        } else {
            a[i] = zero;
        }
    }
}


template<typename Number>
__global__ void initArrayWithOnes(Number* a, size_t nrow) {
    static_assert(std::is_same_v<Number, typename get_real<Number>::type>);

    assert(a);

    auto i = blockDim.x * blockIdx.x + threadIdx.x;

    if(i < nrow) {
        a[i] = Number{1};
    }
}


/**
 * Computes the trace of a matrix A.
 *
 * Optionally, a bit mask can be passed to the function. If a mask is present,
 * then a diagonal element (i,i) is only considered for the trace computation
 * if the i-th value in mask is nonzero.
 *
 * @param[in] mask A bitmask of length nrow indicating which rows to ignore
 * (zero means ignore).
 */
template<typename Number, typename Real = typename get_real<Number>::type>
__global__ void
traceKernel(Number* a, size_t nrow, Real* trace, bool* mask, int mask_present) {
    ENFORCE(a);
    ENFORCE(trace);
    ENFORCE(mask || mask_present == 0);
    ENFORCE(mask_present == 0 || mask_present == 1);

    auto size = nrow * nrow;
    auto i = blockDim.x * blockIdx.x + threadIdx.x;

    if(mask_present == 0) {
        if(i < size) {
            if(i % (nrow + 1) == 0) {
                trace[i % nrow] = a[i].x;
            }
        }
    }
    if(mask_present == 1) {
        if(i < size) {
            if(i % (nrow + 1) == 0) {
                if(mask[i % nrow]) {
                    trace[i % nrow] = a[i].x;
                } else {
                    trace[i % nrow] = 0.0;
                }
            }
        }
    }
}


/*
 * DATA MOVEMENT ROUTINES
 */

extern "C" int cu_createMat(void** d_A, size_t bytecount) {
    ENFORCE(d_A);
    cudaError_t err = cudaMalloc(d_A, bytecount);
    //printf("create mat at GPU Address: %p \n",*d_A);
    return err;
}

extern "C" int cu_cudaFreeAsync(void** d_A) {
    if(!*d_A) {
        return 0;
    }
    int stat = cudaFreeAsync(*d_A, 0);
    *d_A = NULL;
    return stat;
}

extern "C" int cu_cudaMallocAsync(void** d_A, size_t bytecount) {
    ENFORCE(d_A);
    cudaError_t err = cudaMallocAsync(d_A, bytecount, 0);
    //printf("create mat at GPU Address: %p \n",*d_A);
    ENFORCE(err == cudaSuccess);
    return err;
}

extern "C" int cu_copyMatH2D(void* h_A, void* d_A, size_t bytecount) {
    ENFORCE(h_A);
    ENFORCE(d_A);
    // printf("copy %p to %p\n",h_A,d_A);
    cudaError_t err = cudaMemcpy(d_A, h_A, bytecount, cudaMemcpyHostToDevice);
    return err;
}

extern "C" int cu_copyMatH2D_async(void* h_A, void* d_A, size_t bytecount) {
    ENFORCE(h_A);
    ENFORCE(d_A);
    // printf("copy %p to %p\n",h_A,d_A);
    cudaError_t err =
        cudaMemcpyAsync(d_A, h_A, bytecount, cudaMemcpyHostToDevice);
    return err;
}

extern "C" int cu_copyMatD2H(void* h_A, void* d_A, size_t bytecount) {
    ENFORCE(h_A);
    ENFORCE(d_A);

    cudaError_t err = cudaMemcpy(h_A, d_A, bytecount, cudaMemcpyDeviceToHost);
    return err;
}

extern "C" int cu_deleteMat(void** d_A) {
    int stat = 0;
    if(*d_A != NULL) {
        stat = cudaFree(*d_A);
        *d_A = NULL;
    }
    return stat;
}


/*
 * INIT/FINAL ROUTINES
 */

extern "C" int cu_cudaGetDeviceCount(int* count) {
    ENFORCE(count);
    cudaError_t err = cudaGetDeviceCount(count);
    ENFORCE(err == cudaSuccess);
    return err;
}

extern "C" int cu_cudaGetDeviceProperties(int device) {
    cudaDeviceProp prop;
    cudaError_t err = cudaGetDeviceProperties(&prop, device);
    ENFORCE(err == cudaSuccess);

    printf(" Found GPU: Device Name: %s\n", prop.name);
    printf(" TotalMemory: %lu\n", (unsigned long)prop.totalGlobalMem);
    printf(" Shared per block: %lu\n", (unsigned long)prop.sharedMemPerBlock);

    return err;
}

extern "C" int cu_cudaSetDevice(int count) {
    cudaError_t err = cudaSetDevice(count);
    ENFORCE(err == cudaSuccess);
    return err;
}

extern "C" int cu_cublasInit(cublasHandle_t* hcublas) {
    ENFORCE(hcublas);
    cublasStatus_t err = cublasCreate(hcublas);
    ENFORCE(err == CUBLAS_STATUS_SUCCESS);
    if(err != CUBLAS_STATUS_SUCCESS) {
        printf("cublas create error: %d\n", err);
    }
    // printf("hcublas Addr: %p \n",*hcublas);
    return err;
}

extern "C" int cu_cublasFinalize(cublasHandle_t hcublas) {
    cublasStatus_t err = cublasDestroy(hcublas);
    ENFORCE(err == CUBLAS_STATUS_SUCCESS);
    return err;
}

extern "C" int cu_cusolverInit(cusolverDnHandle_t* hcusolver) {
    ENFORCE(hcusolver);
    cusolverStatus_t err = cusolverDnCreate(hcusolver);
    ENFORCE(err == cudaSuccess);
    if(err != 0) {
        printf("cusolver create error: %d\n", err);
    }
    // printf("hcusolver Addr: %p \n",*hcusolver);
    return err;
}

extern "C" int cu_cusolverFinalize(cusolverDnHandle_t hcusolver) {
    cusolverStatus_t err = cusolverDnDestroy(hcusolver);
    ENFORCE(err == cudaSuccess);
    return err;
}


/*
 * MATRIX ROUTINES
 */

/**
 * Multiplies the matrices A, B and adds the result to C.
 *
 * op(.) indicates if the matrix or its complex-conjugate is used. The allowed
 * values for dagger are:
 * * `dagger == 0`: compute C ≔ α A · B + β C
 * * `dagger == 1`: compute C ≔ α A^* · B + β C
 * * `dagger == 2`: compute C ≔ α A · B^* + β C
 *
 * @param[in] m The number of rows of C and op(A).
 * @param[in] n The number of columns of C and op(B).
 * @param[in] k The number of columns of op(A) and the number of rows of op(B).
 * @param[in] dagger A shorthand for various combinations of op(A), op(B).
 */
template<typename Number>
int cu_multMat(
    cublasHandle_t hcublas, size_t m, size_t n, size_t k, const Number* alpha,
    const Number* d_A, const Number* d_B, const Number* beta, Number* d_C,
    int dagger
) {
    ENFORCE(alpha);
    ENFORCE(d_A);
    ENFORCE(d_B);
    ENFORCE(beta);
    ENFORCE(dagger == 0 || dagger == 1 || dagger == 2);

    cublasStatus_t err;
    if(dagger == 0) {
        err = libnegf::cublasGemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, alpha, d_A, m, d_B, k,
            beta, d_C, m
        );
    } else if(dagger == 1) {
        err = libnegf::cublasGemm(
            hcublas, CUBLAS_OP_C, CUBLAS_OP_N, m, n, k, alpha, d_A, k, d_B, k,
            beta, d_C, m
        );
    } else if(dagger == 2) {
        err = libnegf::cublasGemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_C, m, n, k, alpha, d_A, m, d_B, n,
            beta, d_C, m
        );
    } else {
        std::fprintf(stderr, "expected dagger in [0, 1, 2], got %d\n", dagger);
        std::exit(EXIT_FAILURE);
    }
    ENFORCE(err == CUBLAS_STATUS_SUCCESS);
    return err;
}

extern "C" int cu_CmultMat(
    cublasHandle_t hcublas, size_t m, size_t n, size_t k,
    const cuComplex* alpha, const cuComplex* d_A, const cuComplex* d_B,
    const cuComplex* beta, cuComplex* d_C, int dagger
) {
    return cu_multMat(hcublas, m, n, k, alpha, d_A, d_B, beta, d_C, dagger);
}

extern "C" int cu_ZmultMat(
    cublasHandle_t hcublas, size_t m, size_t n, size_t k,
    const cuDoubleComplex* alpha, const cuDoubleComplex* d_A,
    const cuDoubleComplex* d_B, const cuDoubleComplex* beta,
    cuDoubleComplex* d_C, int dagger
) {
    return cu_multMat(hcublas, m, n, k, alpha, d_A, d_B, beta, d_C, dagger);
}


template<typename Number>
int inverse(
    cublasHandle_t hcublas, cusolverDnHandle_t hcusolver, Number* d_A,
    Number* d_Ainv, size_t n
) {
    ENFORCE(hcusolver);
    ENFORCE(d_A);
    ENFORCE(d_Ainv);

    // compute buffer size and prep . memory
    int lwork;
    cusolverStatus_t cusolverStatus =
        libnegf::cusolverDngetrf_bufferSize(hcusolver, n, n, d_A, n, &lwork);
    ENFORCE(cusolverStatus == CUSOLVER_STATUS_SUCCESS);

    // prepare memory on the device
    Number* d_LU;
    cudaError_t cudaStatus = cudaMalloc((void**)&d_LU, n * n * sizeof(Number));
    ENFORCE(cudaStatus == cudaSuccess);

    int* d_pivot;
    cudaStatus = cudaMalloc((void**)&d_pivot, n * sizeof(int));
    ENFORCE(cudaStatus == cudaSuccess);
    int* d_info;
    cudaStatus = cudaMalloc((void**)&d_info, sizeof(int));
    ENFORCE(cudaStatus == cudaSuccess);
    // copy d_LU <- pdA
    auto cublasStatus = libnegf::cublasCopy(hcublas, n * n, d_A, 1, d_LU, 1);
    ENFORCE(cublasStatus == CUBLAS_STATUS_SUCCESS);

    Number* d_work;
    cudaStatus = cudaMalloc((void**)&d_work, lwork * sizeof(Number));
    ENFORCE(cudaStatus == cudaSuccess);

    // LU factorization of d_A , with partial pivoting and row
    // interchanges ; row i is interchanged with row d_pivot ( i );
    cusolverStatus = libnegf::cusolverDngetrf(
        hcusolver, n, n, d_LU, n, d_work, d_pivot, d_info
    );

    // use the LU factorization to solve the system d_LU * x = d_Ainv ;
    // the solution overwrites d_Ainv
    cusolverStatus = libnegf::cusolverDngetrs(
        hcusolver, CUBLAS_OP_N, n, n, d_LU, n, d_pivot, d_Ainv, n, d_info
    );
    ENFORCE(cusolverStatus == CUSOLVER_STATUS_SUCCESS);

    int info_gpu;
    // d_info -> info_gpu
    cudaStatus =
        cudaMemcpy(&info_gpu, d_info, sizeof(int), cudaMemcpyDeviceToHost);
    ENFORCE(cudaStatus == cudaSuccess);

    cudaStatus = cudaFree(d_pivot);
    ENFORCE(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_info);
    ENFORCE(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_work);
    ENFORCE(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_LU);
    ENFORCE(cudaStatus == cudaSuccess);

    return cudaStatus;
}

extern "C" int cu_Cinverse(
    cublasHandle_t hcublas, cusolverDnHandle_t hcusolver, cuComplex* d_A,
    cuComplex* d_Ainv, size_t n
) {
    return inverse(hcublas, hcusolver, d_A, d_Ainv, n);
}

extern "C" int cu_Zinverse(
    cublasHandle_t hcublas, cusolverDnHandle_t hcusolver, cuDoubleComplex* d_A,
    cuDoubleComplex* d_Ainv, size_t n
) {
    return inverse(hcublas, hcusolver, d_A, d_Ainv, n);
}

template<typename Number>
int cu_kernelsum(
    Number* d_C, Number* alpha, Number* d_A, Number* beta, Number* d_B,
    size_t size
) {
    auto num_blocks = (size + BLOCK_SIZE - 1) / BLOCK_SIZE;

    addKernel<<<num_blocks, BLOCK_SIZE>>>(d_C, *alpha, d_A, *beta, d_B, size);
    ENFORCE(cudaGetLastError() == cudaSuccess);

    return 0;
}


extern "C" int cu_Ckernelsum(
    cuComplex* d_C, cuComplex* alpha, cuComplex* d_A, cuComplex* beta,
    cuComplex* d_B, size_t size
) {
    return cu_kernelsum(d_C, alpha, d_A, beta, d_B, size);
}

extern "C" int cu_Zkernelsum(
    cuDoubleComplex* d_C, cuDoubleComplex* alpha, cuDoubleComplex* d_A,
    cuDoubleComplex* beta, cuDoubleComplex* d_B, size_t size
) {
    return cu_kernelsum(d_C, alpha, d_A, beta, d_B, size);
}


/**
 * Computes the sum of α op(A) + β op(B), where op(.) indicates if the matrix
 * or its complex-conjugate is used.
 *
 * The possible options are:
 * * `dagger == 0`: compute C ≔ α A + β B
 * * `dagger == 1`: compute C ≔ α A^* + β B
 * * `dagger == 2`: compute C ≔ α A + B^*
 *
 * @param[in] m The number of rows of op(A) and op(B).
 * @param[in] n The number of columns of op(A) and op(B).
 * @param[in] dagger A shorthand for various combinations of op(A), op(B).
 */
template<typename Number>
int cu_matsum(
    cublasHandle_t hcublas, size_t m, size_t n, const Number* alpha,
    const Number* d_A, const Number* beta, const Number* d_B, Number* d_C,
    int dagger
) {
    ENFORCE(d_A);
    ENFORCE(d_B);
    ENFORCE(d_C);
    ENFORCE(d_A != d_C || dagger == 0 || dagger == 2);
    ENFORCE(d_B != d_C || dagger == 1);

    cublasStatus_t err;
    if(dagger == 0) {
        err = libnegf::cublasGeam(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, m, n, alpha, d_A, m, beta, d_B,
            m, d_C, m
        );
    } else if(dagger == 1) {
        err = libnegf::cublasGeam(
            hcublas, CUBLAS_OP_C, CUBLAS_OP_N, m, n, alpha, d_A, n, beta, d_B,
            m, d_C, m
        );
    } else if(dagger == 2) {
        err = libnegf::cublasGeam(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_C, m, n, alpha, d_A, m, beta, d_B,
            n, d_C, m
        );
    } else {
        std::fprintf(stderr, "expected dagger in [0, 1, 2], got %d\n", dagger);
        std::exit(EXIT_FAILURE);
    }

    return err;
}

extern "C" int cu_Cmatsum(
    cublasHandle_t hcublas, size_t m, size_t n, const cuComplex* alpha,
    const cuComplex* d_A, const cuComplex* beta, const cuComplex* d_B,
    cuComplex* d_C, int dagger
) {
    return cu_matsum(hcublas, m, n, alpha, d_A, beta, d_B, d_C, dagger);
}

extern "C" int cu_Zmatsum(
    cublasHandle_t hcublas, size_t m, size_t n, const cuDoubleComplex* alpha,
    const cuDoubleComplex* d_A, const cuDoubleComplex* beta,
    const cuDoubleComplex* d_B, cuDoubleComplex* d_C, int dagger
) {
    return cu_matsum(hcublas, m, n, alpha, d_A, beta, d_B, d_C, dagger);
}

extern "C" int cu_Cinitmat(cuComplex* d_A, size_t nrow) {
    ENFORCE(d_A);
    auto size = nrow * nrow;
    auto num_blocks = (size + BLOCK_SIZE - 1) / BLOCK_SIZE;

    initKernel<<<num_blocks, BLOCK_SIZE>>>(d_A, nrow);
    ENFORCE(cudaGetLastError() == cudaSuccess);

    return 0;
}

extern "C" int cu_Zinitmat(cuDoubleComplex* d_A, size_t nrow) {
    ENFORCE(d_A);
    auto size = nrow * nrow;
    auto num_blocks = (size + BLOCK_SIZE - 1) / BLOCK_SIZE;

    initKernel<<<num_blocks, BLOCK_SIZE>>>(d_A, nrow);
    ENFORCE(cudaGetLastError() == cudaSuccess);

    return 0;
}


template<typename Number, typename Real = typename get_real<Number>::type>
Real trace(
    cublasHandle_t hcublas, Number* d_A, size_t nrow, void* h_mask,
    int mask_present
) {
    auto size = nrow * nrow;
    auto num_blocks = (size + BLOCK_SIZE - 1) / BLOCK_SIZE;
    Real* d_work;
    cudaError_t cudaStatus = cudaMalloc((void**)&d_work, nrow * sizeof(Real));
    Real* d_iden;
    cudaStatus = cudaMalloc((void**)&d_iden, nrow * sizeof(Real));
    ENFORCE(cudaStatus == cudaSuccess);
    bool* d_mask;
    cudaStatus = cudaMalloc((void**)&d_mask, nrow * sizeof(bool));
    ENFORCE(cudaStatus == cudaSuccess);
    if(h_mask) {
        cudaStatus = cudaMemcpy(
            d_mask, h_mask, nrow * sizeof(bool), cudaMemcpyHostToDevice
        );
        ENFORCE(cudaStatus == cudaSuccess);
    }

    initArrayWithOnes<<<num_blocks, BLOCK_SIZE>>>(d_iden, nrow);
    traceKernel<<<num_blocks, BLOCK_SIZE>>>(
        d_A, nrow, d_work, d_mask, mask_present
    );

    Real result;
    cublasStatus_t err =
        libnegf::cublasDot(hcublas, nrow, d_iden, 1, d_work, 1, &result);
    ENFORCE(err == CUBLAS_STATUS_SUCCESS);

    cudaStatus = cudaFree(d_work);
    cudaStatus = cudaFree(d_iden);
    cudaStatus = cudaFree(d_mask);

    return result;
}

extern "C" float cu_Ctrace(
    cublasHandle_t hcublas, cuComplex* d_A, size_t nrow, void* h_mask,
    int mask_present
) {
    return trace(hcublas, d_A, nrow, h_mask, mask_present);
}


extern "C" double cu_Ztrace(
    cublasHandle_t hcublas, cuDoubleComplex* d_A, size_t nrow, void* h_mask,
    int mask_present
) {
    return trace(hcublas, d_A, nrow, h_mask, mask_present);
}


extern "C" int cu_Cmatcopy(
    cublasHandle_t hcublas, const cuComplex* d_A, cuComplex* d_B,
    size_t num_elements
) {
    ENFORCE(d_A);
    ENFORCE(d_B);
    ENFORCE(num_elements <= INT_MAX);

    auto err = cublasCcopy(hcublas, num_elements, d_A, 1, d_B, 1);
    ENFORCE(err == CUBLAS_STATUS_SUCCESS);
    return err;
}

extern "C" int cu_Zmatcopy(
    cublasHandle_t hcublas, cuDoubleComplex* d_A, cuDoubleComplex* d_B,
    size_t num_elements
) {
    ENFORCE(d_A);
    ENFORCE(d_B);
    ENFORCE(num_elements <= INT_MAX);

    auto err = cublasZcopy(hcublas, num_elements, d_A, 1, d_B, 1);
    ENFORCE(err == CUBLAS_STATUS_SUCCESS);
    return err;
}

extern "C" int
cu_Casum(cublasHandle_t hcublas, void* d_A, float* summ, size_t n) {
    cuComplex* pdA = (cuComplex*)d_A;

    auto err = cublasScasum(hcublas, n, pdA, 1, summ);
    ENFORCE(err == CUBLAS_STATUS_SUCCESS);
    return err;
}

extern "C" int
cu_Zasum(cublasHandle_t hcublas, void* d_A, double* summ, size_t n) {
    cuDoubleComplex* pdA = (cuDoubleComplex*)d_A;

    auto err = cublasDzasum(hcublas, n, pdA, 1, summ);
    ENFORCE(err == CUBLAS_STATUS_SUCCESS);
    return err;
}


template<typename Number>
int decimation(
    cublasHandle_t hcublas, cusolverDnHandle_t hcusolver, Number* h_Go_out,
    Number* h_Ao_in, Number* h_Bo_in, Number* h_Co_in, size_t n, int tf32,
    int* ncyc, float SGFACC
) {
    ENFORCE(h_Go_out);
    ENFORCE(h_Ao_in);
    ENFORCE(h_Bo_in);
    ENFORCE(h_Co_in);
    ENFORCE(tf32 == 0 || tf32 == 1);
    ENFORCE(ncyc);
    ENFORCE(SGFACC > 0.0);

    auto num_elements = n * n;
    auto num_blocks = (num_elements + BLOCK_SIZE - 1) / BLOCK_SIZE;

    using Real = typename get_real<Number>::type;
    auto one = Number{Real{1}};
    auto mone = Number{-Real{1}};
    auto zero = Number{Real{0}};

    Number* d_Ao;
    cudaError_t cudaStatus =
        cudaMalloc((void**)&d_Ao, num_elements * sizeof(Number));
    ENFORCE(cudaStatus == cudaSuccess);

    Number* d_Bo;
    cudaStatus = cudaMalloc((void**)&d_Bo, num_elements * sizeof(Number));
    ENFORCE(cudaStatus == cudaSuccess);

    Number* d_Co;
    cudaStatus = cudaMalloc((void**)&d_Co, num_elements * sizeof(Number));
    ENFORCE(cudaStatus == cudaSuccess);

    cudaStatus = cudaMemcpy(
        d_Ao, h_Ao_in, n * n * sizeof(Number), cudaMemcpyHostToDevice
    );
    ENFORCE(cudaStatus == cudaSuccess);
    cudaStatus = cudaMemcpy(
        d_Bo, h_Bo_in, n * n * sizeof(Number), cudaMemcpyHostToDevice
    );
    ENFORCE(cudaStatus == cudaSuccess);
    cudaStatus = cudaMemcpy(
        d_Co, h_Co_in, n * n * sizeof(Number), cudaMemcpyHostToDevice
    );
    ENFORCE(cudaStatus == cudaSuccess);

    cublasStatus_t cublasStatus =
        cublasSetPointerMode(hcublas, CUBLAS_POINTER_MODE_HOST);
    ENFORCE(cublasStatus == cudaSuccess);

    if(tf32 == 1) {
        cublasStatus = cublasSetMathMode(hcublas, CUBLAS_TENSOR_OP_MATH);
        ENFORCE(cublasStatus == cudaSuccess);
    }

    Number* d_Ao_s;
    cudaStatus = cudaMalloc((void**)&d_Ao_s, n * n * sizeof(Number));
    ENFORCE(cudaStatus == cudaSuccess);
    Number* d_C1;
    cudaStatus = cudaMalloc((void**)&d_C1, n * n * sizeof(Number));
    ENFORCE(cudaStatus == cudaSuccess);
    Number* d_Go;
    cudaStatus = cudaMalloc((void**)&d_Go, n * n * sizeof(Number));
    ENFORCE(cudaStatus == cudaSuccess);
    int* d_pivot;
    cudaStatus = cudaMalloc((void**)&d_pivot, n * sizeof(int));
    ENFORCE(cudaStatus == cudaSuccess);
    Number* d_T;
    cudaStatus = cudaMalloc((void**)&d_T, n * n * sizeof(Number));
    ENFORCE(cudaStatus == cudaSuccess);
    Number* d_Self;
    cudaStatus = cudaMalloc((void**)&d_Self, n * n * sizeof(Number));
    ENFORCE(cudaStatus == cudaSuccess);
    int* d_info;
    cudaStatus = cudaMalloc((void**)&d_info, sizeof(int));
    ENFORCE(cudaStatus == cudaSuccess);

    int lwork;
    cusolverStatus_t cusolverStatus =
        libnegf::cusolverDngetrf_bufferSize(hcusolver, n, n, d_Self, n, &lwork);
    ENFORCE(cusolverStatus == cudaSuccess);
    Number* d_work;
    cudaStatus = cudaMalloc((void**)&d_work, lwork * sizeof(Number));
    ENFORCE(cudaStatus == cudaSuccess);

    cublasStatus = libnegf::cublasCopy(hcublas, n * n, d_Ao, 1, d_Ao_s, 1);
    ENFORCE(cublasStatus == cudaSuccess);

    bool okCo = false;
    for(int i1 = 1; i1 <= 300; i1++) {
        *ncyc = i1;

        initKernel<<<num_blocks, BLOCK_SIZE>>>(d_Go, n);

        cublasStatus = libnegf::cublasCopy(hcublas, n * n, d_Ao, 1, d_Self, 1);
        ENFORCE(cublasStatus == cudaSuccess);

        cusolverStatus = libnegf::cusolverDngetrf(
            hcusolver, n, n, d_Self, n, d_work, d_pivot, d_info
        );
        ENFORCE(cusolverStatus == cudaSuccess);
        cusolverStatus = libnegf::cusolverDngetrs(
            hcusolver, CUBLAS_OP_N, n, n, d_Self, n, d_pivot, d_Go, n, d_info
        );
        ENFORCE(cusolverStatus == cudaSuccess);

        cublasStatus = libnegf::cublasGemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, &one, d_Go, n, d_Co, n,
            &zero, d_T, n
        );
        ENFORCE(cublasStatus == cudaSuccess);

        cublasStatus = libnegf::cublasGemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, &one, d_Co, n, d_T, n,
            &zero, d_C1, n
        );
        ENFORCE(cublasStatus == cudaSuccess);

        Real summ;
        cublasStatus = libnegf::cublasAsum(hcublas, n * n, d_C1, 1, &summ);
        ENFORCE(cublasStatus == cudaSuccess);
        // printf("loop it= %d , summ= %f \n ", i1, summ);

        if(summ <= SGFACC) {
            if(okCo) {
                break;
            } else {
                okCo = true;
            }
        } else {
            okCo = false;
        }

        cublasStatus = libnegf::cublasGemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, &one, d_Bo, n, d_T, n,
            &zero, d_Self, n
        );
        ENFORCE(cublasStatus == cudaSuccess);

        cublasStatus =
            libnegf::cublasAxpy(hcublas, n * n, &mone, d_Self, 1, d_Ao_s, 1);
        ENFORCE(cublasStatus == cudaSuccess);
        cublasStatus =
            libnegf::cublasAxpy(hcublas, n * n, &mone, d_Self, 1, d_Ao, 1);
        ENFORCE(cublasStatus == cudaSuccess);

        cublasStatus = libnegf::cublasGemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, &one, d_Go, n, d_Bo, n,
            &zero, d_T, n
        );
        ENFORCE(cublasStatus == cudaSuccess);
        cublasStatus = libnegf::cublasGemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, &mone, d_Co, n, d_T, n,
            &one, d_Ao, n
        );
        ENFORCE(cublasStatus == cudaSuccess);

        cublasStatus = libnegf::cublasCopy(hcublas, n * n, d_C1, 1, d_Co, 1);
        ENFORCE(cublasStatus == cudaSuccess);

        cublasStatus = libnegf::cublasGemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, &one, d_Bo, n, d_T, n,
            &zero, d_C1, n
        );
        ENFORCE(cublasStatus == cudaSuccess);

        cublasStatus = libnegf::cublasCopy(hcublas, n * n, d_C1, 1, d_Bo, 1);
        ENFORCE(cublasStatus == cudaSuccess);
    }

    initKernel<<<num_blocks, BLOCK_SIZE>>>(d_Go, n);
    cublasStatus = libnegf::cublasCopy(hcublas, n * n, d_Ao_s, 1, d_Self, 1);
    ENFORCE(cublasStatus == cudaSuccess);
    cusolverStatus = libnegf::cusolverDngetrf(
        hcusolver, n, n, d_Self, n, d_work, d_pivot, d_info
    );
    ENFORCE(cusolverStatus == cudaSuccess);
    cusolverStatus = libnegf::cusolverDngetrs(
        hcusolver, CUBLAS_OP_N, n, n, d_Self, n, d_pivot, d_Go, n, d_info
    );
    ENFORCE(cusolverStatus == cudaSuccess);

    cudaStatus = cudaMemcpy(
        h_Go_out, d_Go, n * n * sizeof(Number), cudaMemcpyDeviceToHost
    );
    ENFORCE(cudaStatus == cudaSuccess);

    cudaStatus = cudaFree(d_pivot);
    ENFORCE(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_info);
    ENFORCE(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_Ao);
    ENFORCE(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_Bo);
    ENFORCE(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_Co);
    ENFORCE(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_Go);
    ENFORCE(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_Ao_s);
    ENFORCE(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_C1);
    ENFORCE(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_T);
    ENFORCE(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_Self);
    ENFORCE(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_work);
    ENFORCE(cudaStatus == cudaSuccess);

    return cudaStatus;
}


extern "C" int cu_Cdecimation(
    cublasHandle_t hcublas, cusolverDnHandle_t hcusolver, cuComplex* h_Go_out,
    cuComplex* h_Ao_in, cuComplex* h_Bo_in, cuComplex* h_Co_in, size_t n,
    int tf32, int* ncyc, float SGFACC
) {
    return decimation(
        hcublas, hcusolver, h_Go_out, h_Ao_in, h_Bo_in, h_Co_in, n, tf32, ncyc,
        SGFACC
    );
}

extern "C" int cu_Zdecimation(
    cublasHandle_t hcublas, cusolverDnHandle_t hcusolver,
    cuDoubleComplex* h_Go_out, cuDoubleComplex* h_Ao_in,
    cuDoubleComplex* h_Bo_in, cuDoubleComplex* h_Co_in, size_t n, int tf32,
    int* ncyc, double SGFACC
) {
    return decimation(
        hcublas, hcusolver, h_Go_out, h_Ao_in, h_Bo_in, h_Co_in, n, tf32, ncyc,
        SGFACC
    );
}


extern "C" int cu_meminfo(size_t* freemem, size_t* totalmem) {
    cudaError_t cudaStatus;
    cudaStatus = cudaDeviceSynchronize();
    cudaStatus = cudaMemGetInfo(freemem, totalmem);
    ENFORCE(cudaStatus == cudaSuccess);
    return cudaStatus;
}
