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
#include <cusolverDn.h>

#include <cassert>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

#include <libnegf/cublas_v2.h>


constexpr auto BLOCK_SIZE = std::size_t{1024};


template<typename>
struct get_real {};

template<>
struct get_real<cuComplex> {
    using type = float;
};

template<>
struct get_real<cuDoubleComplex> {
    using type = double;
};


template<typename Number>
__global__ void addKernel(
    Number* c, Number alpha, const Number* a, Number beta, const Number* b,
    size_t size
) {
    assert(c);
    assert(a);
    assert(b);

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


template<typename Number>
__global__ void initKernel(Number* a, size_t nrow) {
    assert(a);

    auto size = nrow * nrow;
    auto i = blockDim.x * blockIdx.x + threadIdx.x;

    if(i < size) {
        if(i % (nrow + 1) == 0) {
            a[i].x = 1.0;
            a[i].y = 0.0;
        } else {
            a[i].x = 0.0;
            a[i].y = 0.0;
        }
    }
}


__global__ void DinitKernel(double* a, size_t nrow) {
    assert(a);

    auto i = blockDim.x * blockIdx.x + threadIdx.x;

    if(i < nrow) {
        a[i] = 1.0;
    }
}

__global__ void SinitKernel(float* a, size_t nrow) {
    assert(a);

    auto i = blockDim.x * blockIdx.x + threadIdx.x;

    if(i < nrow) {
        a[i] = 1.0;
    }
}


template<typename Number, typename Real = typename get_real<Number>::type>
__device__ void
traceKernel(Number* a, size_t nrow, Real* trace, bool* mask, int mask_present) {
    assert(a);
    assert(trace);
    assert(mask || mask_present == 0);
    assert(mask_present == 0 || mask_present == 1);

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

__global__ void CtraceKernel(
    cuComplex* a, size_t nrow, float* trace, bool* mask, int mask_present
) {
    traceKernel(a, nrow, trace, mask, mask_present);
}

__global__ void ZtraceKernel(
    cuDoubleComplex* a, size_t nrow, double* trace, bool* mask, int mask_present
) {
    traceKernel(a, nrow, trace, mask, mask_present);
}


/*~-~-~-~-~-~-~-~-~-~-~-~-~-~ DATA MOVEMENT  ROUTINES
 * -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-*/

extern "C" int cu_createMat(void** d_A, size_t bytecount) {
    assert(d_A);
    cudaError_t err = cudaMalloc(d_A, bytecount);
    //printf("create mat at GPU Address: %p \n",*d_A);
    return err;
}

extern "C" int cu_copyMatH2D(void* h_A, void* d_A, size_t bytecount) {
    assert(h_A);
    assert(d_A);
    // printf("copy %p to %p\n",h_A,d_A);
    cudaError_t err = cudaMemcpy(d_A, h_A, bytecount, cudaMemcpyHostToDevice);
    return err;
}

extern "C" int cu_copyMatD2H(void* h_A, void* d_A, size_t bytecount) {
    assert(h_A);
    assert(d_A);

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

/*~-~-~-~-~-~-~-~-~-~-~-~-~-~ INIT/FINAL ROUTINES
 * -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-*/

extern "C" int cu_cudaGetDeviceCount(int* count) {
    assert(count);
    cudaError_t err = cudaGetDeviceCount(count);
    assert(err == cudaSuccess);
    return err;
}

extern "C" int cu_cudaGetDeviceProperties(int device) {
    cudaDeviceProp prop;
    cudaError_t err = cudaGetDeviceProperties(&prop, device);
    assert(err == cudaSuccess);

    printf(" Found GPU: Device Name: %s\n", prop.name);
    printf(" TotalMemory: %lu\n", (unsigned long)prop.totalGlobalMem);
    printf(" Shared per block: %lu\n", (unsigned long)prop.sharedMemPerBlock);

    return err;
}

extern "C" int cu_cudaSetDevice(int count) {
    cudaError_t err = cudaSetDevice(count);
    assert(err == cudaSuccess);
    return err;
}

extern "C" int cu_cublasInit(cublasHandle_t* hcublas) {
    assert(hcublas);
    cublasStatus_t err = cublasCreate(hcublas);
    assert(err == CUBLAS_STATUS_SUCCESS);
    if(err != CUBLAS_STATUS_SUCCESS) {
        printf("cublas create error: %d\n", err);
    }
    // printf("hcublas Addr: %p \n",*hcublas);
    return err;
}

extern "C" int cu_cublasFinalize(cublasHandle_t hcublas) {
    cublasStatus_t err = cublasDestroy(hcublas);
    assert(err == CUBLAS_STATUS_SUCCESS);
    return err;
}

extern "C" int cu_cusolverInit(cusolverDnHandle_t* hcusolver) {
    assert(hcusolver);
    cusolverStatus_t err = cusolverDnCreate(hcusolver);
    assert(err == cudaSuccess);
    if(err != 0) {
        printf("cusolver create error: %d\n", err);
    }
    // printf("hcusolver Addr: %p \n",*hcusolver);
    return err;
}

extern "C" int cu_cusolverFinalize(cusolverDnHandle_t hcusolver) {
    cusolverStatus_t err = cusolverDnDestroy(hcusolver);
    assert(err == cudaSuccess);
    return err;
}


/*
 * MATRIX ROUTINES
 */

/* C = alpha op(A) op(B) + beta C
 * m: #rows of op(A)
 * n: #cols of op(B)
 * k: #cols of op(A) = #rows of op(B)
 */
template<typename Number>
int cu_multMat(
    cublasHandle_t hcublas, size_t m, size_t n, size_t k, const Number* alpha,
    const Number* d_A, const Number* d_B, const Number* beta, Number* d_C,
    int dagger
) {
    assert(alpha);
    assert(d_A);
    assert(d_B);
    assert(beta);
    assert(dagger == 0 || dagger == 1 || dagger == 2);

    cublasStatus_t err;
    if(dagger == 0) {
        err = cublasGemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, alpha, d_A, m, d_B, k,
            beta, d_C, m
        );
    } else if(dagger == 1) {
        err = cublasGemm(
            hcublas, CUBLAS_OP_C, CUBLAS_OP_N, m, n, k, alpha, d_A, k, d_B, k,
            beta, d_C, m
        );
    } else if(dagger == 2) {
        err = cublasGemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_C, m, n, k, alpha, d_A, m, d_B, n,
            beta, d_C, m
        );
    } else {
        std::fprintf(stderr, "expected dagger in [0, 1, 2], got %d\n", dagger);
        std::exit(EXIT_FAILURE);
    }
    assert(err == CUBLAS_STATUS_SUCCESS);
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


extern "C" int cu_Cinverse(
    cublasHandle_t hcublas, cusolverDnHandle_t hcusolver, void* d_A,
    void* d_Ainv, size_t n
) {
    assert(hcusolver);
    assert(d_A);
    assert(d_Ainv);

    // compute buffer size and prep . memory
    cuComplex* pdA = (cuComplex*)d_A;
    int lwork;
    cusolverStatus_t cusolverStatus =
        cusolverDnCgetrf_bufferSize(hcusolver, n, n, pdA, n, &lwork);
    assert(cusolverStatus == CUSOLVER_STATUS_SUCCESS);

    // prepare memory on the device
    cuComplex* d_LU;
    cudaError_t cudaStatus =
        cudaMalloc((void**)&d_LU, n * n * sizeof(cuComplex));
    assert(cudaStatus == cudaSuccess);

    int* d_pivot;
    cudaStatus = cudaMalloc((void**)&d_pivot, n * sizeof(int));
    assert(cudaStatus == cudaSuccess);
    int* d_info;
    cudaStatus = cudaMalloc((void**)&d_info, sizeof(int));
    assert(cudaStatus == cudaSuccess);
    // copy d_LU <- pdA
    cublasStatus_t cublasStatus = cublasCcopy(hcublas, n * n, pdA, 1, d_LU, 1);
    assert(cublasStatus == CUBLAS_STATUS_SUCCESS);

    cuComplex* d_work;
    cudaStatus = cudaMalloc((void**)&d_work, lwork * sizeof(cuComplex));
    assert(cudaStatus == cudaSuccess);

    // LU factorization of d_A , with partial pivoting and row
    // interchanges ; row i is interchanged with row d_pivot ( i );
    cusolverStatus =
        cusolverDnCgetrf(hcusolver, n, n, d_LU, n, d_work, d_pivot, d_info);

    // use the LU factorization to solve the system d_LU * x = d_Ainv ;
    // the solution overwrites d_Ainv
    cuComplex* pdAinv = (cuComplex*)d_Ainv;
    cusolverStatus = cusolverDnCgetrs(
        hcusolver, CUBLAS_OP_N, n, n, d_LU, n, d_pivot, pdAinv, n, d_info
    );
    assert(cusolverStatus == CUSOLVER_STATUS_SUCCESS);

    int info_gpu;
    // d_info -> info_gpu
    cudaStatus =
        cudaMemcpy(&info_gpu, d_info, sizeof(int), cudaMemcpyDeviceToHost);
    assert(cudaStatus == cudaSuccess);

    cudaStatus = cudaFree(d_pivot);
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_info);
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_work);
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_LU);
    assert(cudaStatus == cudaSuccess);

    return cudaStatus;
}

extern "C" int cu_Zinverse(
    cublasHandle_t hcublas, cusolverDnHandle_t hcusolver, void* d_A,
    void* d_Ainv, size_t n
) {
    assert(hcusolver);
    assert(d_A);
    assert(d_Ainv);

    // compute buffer size and prep . memory
    cuDoubleComplex* pdA = (cuDoubleComplex*)d_A;
    int lwork;
    cusolverStatus_t cusolverStatus =
        cusolverDnZgetrf_bufferSize(hcusolver, n, n, pdA, n, &lwork);
    assert(cusolverStatus == CUSOLVER_STATUS_SUCCESS);

    // prepare memory on the device
    cuDoubleComplex* d_LU;
    cudaError_t cudaStatus =
        cudaMalloc((void**)&d_LU, n * n * sizeof(cuDoubleComplex));
    assert(cudaStatus == cudaSuccess);

    int* d_pivot;
    cudaStatus = cudaMalloc((void**)&d_pivot, n * sizeof(int));
    assert(cudaStatus == cudaSuccess);
    int* d_info;
    cudaStatus = cudaMalloc((void**)&d_info, sizeof(int));
    assert(cudaStatus == cudaSuccess);
    // copy d_LU <- pdA
    cublasStatus_t cublasStatus = cublasZcopy(hcublas, n * n, pdA, 1, d_LU, 1);
    assert(cublasStatus == CUBLAS_STATUS_SUCCESS);

    cuDoubleComplex* d_work;
    cudaStatus = cudaMalloc((void**)&d_work, lwork * sizeof(cuDoubleComplex));
    assert(cudaStatus == cudaSuccess);

    // LU factorization of d_A , with partial pivoting and row
    // interchanges ; row i is interchanged with row d_pivot ( i );
    cusolverStatus =
        cusolverDnZgetrf(hcusolver, n, n, d_LU, n, d_work, d_pivot, d_info);

    // use the LU factorization to solve the system d_LU * x = d_Ainv ;
    // the solution overwrites d_Ainv
    cuDoubleComplex* pdAinv = (cuDoubleComplex*)d_Ainv;
    cusolverStatus = cusolverDnZgetrs(
        hcusolver, CUBLAS_OP_N, n, n, d_LU, n, d_pivot, pdAinv, n, d_info
    );
    assert(cusolverStatus == CUSOLVER_STATUS_SUCCESS);

    int info_gpu;
    // d_info -> info_gpu
    cudaStatus =
        cudaMemcpy(&info_gpu, d_info, sizeof(int), cudaMemcpyDeviceToHost);
    assert(cudaStatus == cudaSuccess);

    cudaStatus = cudaFree(d_pivot);
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_info);
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_work);
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_LU);
    assert(cudaStatus == cudaSuccess);

    return cudaStatus;
}

template<typename Number>
int cu_kernelsum(
    Number* d_C, Number* alpha, Number* d_A, Number* beta, Number* d_B,
    size_t size
) {
    auto num_blocks = (size + BLOCK_SIZE - 1) / BLOCK_SIZE;

    addKernel<<<num_blocks, BLOCK_SIZE>>>(d_C, *alpha, d_A, *beta, d_B, size);
    assert(cudaGetLastError() == cudaSuccess);

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


template<typename Number>
int cu_matsum(
    cublasHandle_t hcublas, size_t m, size_t n, const Number* alpha,
    const Number* d_A, const Number* beta, const Number* d_B, Number* d_C,
    int dagger
) {
    cublasStatus_t err;
    if(dagger == 0) {
        err = cublasGeam(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, m, n, alpha, d_A, m, beta, d_B,
            m, d_C, m
        );
    } else if(dagger == 1) {
        err = cublasGeam(
            hcublas, CUBLAS_OP_C, CUBLAS_OP_N, m, n, alpha, d_A, n, beta, d_B,
            m, d_C, m
        );
    } else if(dagger == 2) {
        err = cublasGeam(
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
    assert(d_A);
    auto size = nrow * nrow;
    auto num_blocks = (size + BLOCK_SIZE - 1) / BLOCK_SIZE;

    initKernel<<<num_blocks, BLOCK_SIZE>>>(d_A, nrow);
    assert(cudaGetLastError() == cudaSuccess);

    return 0;
}

extern "C" int cu_Zinitmat(cuDoubleComplex* d_A, size_t nrow) {
    assert(d_A);
    auto size = nrow * nrow;
    auto num_blocks = (size + BLOCK_SIZE - 1) / BLOCK_SIZE;

    initKernel<<<num_blocks, BLOCK_SIZE>>>(d_A, nrow);
    assert(cudaGetLastError() == cudaSuccess);

    return 0;
}


extern "C" float cu_Ctrace(
    cublasHandle_t hcublas, void* d_A, size_t nrow, void* h_mask,
    int mask_present
) {
    cuComplex* pdA = (cuComplex*)d_A;
    auto size = nrow * nrow;
    auto num_blocks = (size + BLOCK_SIZE - 1) / BLOCK_SIZE;
    float* d_work;
    cudaError_t cudaStatus = cudaMalloc((void**)&d_work, nrow * sizeof(float));
    float* d_iden;
    cudaStatus = cudaMalloc((void**)&d_iden, nrow * sizeof(float));
    assert(cudaStatus == cudaSuccess);
    bool* d_mask;
    cudaStatus = cudaMalloc((void**)&d_mask, nrow * sizeof(bool));
    assert(cudaStatus == cudaSuccess);
    if(h_mask) {
        cudaStatus = cudaMemcpy(
            d_mask, h_mask, nrow * sizeof(bool), cudaMemcpyHostToDevice
        );
        assert(cudaStatus == cudaSuccess);
    }

    SinitKernel<<<num_blocks, BLOCK_SIZE>>>(d_iden, nrow);
    CtraceKernel<<<num_blocks, BLOCK_SIZE>>>(
        pdA, nrow, d_work, d_mask, mask_present
    );

    float result;
    cublasStatus_t err =
        cublasSdot(hcublas, nrow, d_iden, 1, d_work, 1, &result);
    assert(err == CUBLAS_STATUS_SUCCESS);

    cudaStatus = cudaFree(d_work);
    cudaStatus = cudaFree(d_iden);
    cudaStatus = cudaFree(d_mask);

    return result;
}

extern "C" double cu_Ztrace(
    cublasHandle_t hcublas, void* d_A, size_t nrow, void* h_mask,
    int mask_present
) {
    assert(d_A);
    assert(h_mask || mask_present == 0);
    assert(mask_present == 0 || mask_present == 1);

    cuDoubleComplex* pdA = (cuDoubleComplex*)d_A;
    auto size = nrow * nrow;
    auto num_blocks = (size + BLOCK_SIZE - 1) / BLOCK_SIZE;
    double* d_work;
    cudaError_t cudaStatus = cudaMalloc((void**)&d_work, nrow * sizeof(double));
    assert(cudaStatus == cudaSuccess);
    double* d_iden;
    cudaStatus = cudaMalloc((void**)&d_iden, nrow * sizeof(double));
    assert(cudaStatus == cudaSuccess);
    bool* d_mask;
    cudaStatus = cudaMalloc((void**)&d_mask, nrow * sizeof(bool));
    assert(cudaStatus == cudaSuccess);
    if(h_mask) {
        cudaStatus = cudaMemcpy(
            d_mask, h_mask, nrow * sizeof(bool), cudaMemcpyHostToDevice
        );
        assert(cudaStatus == cudaSuccess);
    }

    DinitKernel<<<num_blocks, BLOCK_SIZE>>>(d_iden, nrow);
    assert(cudaGetLastError() == cudaSuccess);
    ZtraceKernel<<<num_blocks, BLOCK_SIZE>>>(
        pdA, nrow, d_work, d_mask, mask_present
    );
    assert(cudaGetLastError() == cudaSuccess);

    double result;
    cublasStatus_t err =
        cublasDdot(hcublas, nrow, d_iden, 1, d_work, 1, &result);
    assert(err == CUBLAS_STATUS_SUCCESS);

    cudaStatus = cudaFree(d_work);
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_iden);
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_mask);
    assert(cudaStatus == cudaSuccess);

    return result;
}

extern "C" int cu_Cmatcopy(
    cublasHandle_t hcublas, const cuComplex* d_A, cuComplex* d_B,
    size_t num_elements
) {
    assert(d_A);
    assert(d_B);
    assert(num_elements <= INT_MAX);

    auto err = cublasCcopy(hcublas, num_elements, d_A, 1, d_B, 1);
    assert(err == CUBLAS_STATUS_SUCCESS);
    return err;
}

extern "C" int cu_Zmatcopy(
    cublasHandle_t hcublas, cuDoubleComplex* d_A, cuDoubleComplex* d_B,
    size_t num_elements
) {
    assert(d_A);
    assert(d_B);
    assert(num_elements <= INT_MAX);

    auto err = cublasZcopy(hcublas, num_elements, d_A, 1, d_B, 1);
    assert(err == CUBLAS_STATUS_SUCCESS);
    return err;
}

extern "C" int
cu_Casum(cublasHandle_t hcublas, void* d_A, float* summ, size_t n) {
    cuComplex* pdA = (cuComplex*)d_A;

    auto err = cublasScasum(hcublas, n, pdA, 1, summ);
    assert(err == CUBLAS_STATUS_SUCCESS);
    return err;
}

extern "C" int
cu_Zasum(cublasHandle_t hcublas, void* d_A, double* summ, size_t n) {
    cuDoubleComplex* pdA = (cuDoubleComplex*)d_A;

    auto err = cublasDzasum(hcublas, n, pdA, 1, summ);
    assert(err == CUBLAS_STATUS_SUCCESS);
    return err;
}

extern "C" int cu_Cdecimation(
    cublasHandle_t hcublas, cusolverDnHandle_t hcusolver, void* h_Go_out,
    void* h_Ao_in, void* h_Bo_in, void* h_Co_in, size_t n, int tf32, int* ncyc,
    cuComplex* one, cuComplex* mone, cuComplex* zero, float SGFACC
) {
    assert(h_Go_out);
    assert(h_Ao_in);
    assert(h_Bo_in);
    assert(h_Co_in);
    assert(tf32 == 0 || tf32 == 1);
    assert(ncyc);
    assert(one);
    assert(mone);
    assert(zero);
    assert(SGFACC > 0.0);

    auto num_elements = n * n;
    auto num_blocks = (num_elements + BLOCK_SIZE - 1) / BLOCK_SIZE;

    cuComplex* d_Ao;
    cudaError_t cudaStatus =
        cudaMalloc((void**)&d_Ao, num_elements * sizeof(cuComplex));
    assert(cudaStatus == cudaSuccess);

    cuComplex* d_Bo;
    cudaStatus = cudaMalloc((void**)&d_Bo, num_elements * sizeof(cuComplex));
    assert(cudaStatus == cudaSuccess);

    cuComplex* d_Co;
    cudaStatus = cudaMalloc((void**)&d_Co, num_elements * sizeof(cuComplex));
    assert(cudaStatus == cudaSuccess);

    cudaStatus = cudaMemcpy(
        d_Ao, h_Ao_in, n * n * sizeof(cuComplex), cudaMemcpyHostToDevice
    );
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaMemcpy(
        d_Bo, h_Bo_in, n * n * sizeof(cuComplex), cudaMemcpyHostToDevice
    );
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaMemcpy(
        d_Co, h_Co_in, n * n * sizeof(cuComplex), cudaMemcpyHostToDevice
    );
    assert(cudaStatus == cudaSuccess);

    cublasStatus_t cublasStatus =
        cublasSetPointerMode(hcublas, CUBLAS_POINTER_MODE_HOST);
    assert(cublasStatus == cudaSuccess);

    if(tf32 == 1) {
        cublasStatus = cublasSetMathMode(hcublas, CUBLAS_TENSOR_OP_MATH);
        assert(cublasStatus == cudaSuccess);
    }

    cuComplex* d_Ao_s;
    cudaStatus = cudaMalloc((void**)&d_Ao_s, n * n * sizeof(cuComplex));
    assert(cudaStatus == cudaSuccess);
    cuComplex* d_C1;
    cudaStatus = cudaMalloc((void**)&d_C1, n * n * sizeof(cuComplex));
    assert(cudaStatus == cudaSuccess);
    cuComplex* d_Go;
    cudaStatus = cudaMalloc((void**)&d_Go, n * n * sizeof(cuComplex));
    assert(cudaStatus == cudaSuccess);
    int* d_pivot;
    cudaStatus = cudaMalloc((void**)&d_pivot, n * sizeof(int));
    assert(cudaStatus == cudaSuccess);
    cuComplex* d_T;
    cudaStatus = cudaMalloc((void**)&d_T, n * n * sizeof(cuComplex));
    assert(cudaStatus == cudaSuccess);
    cuComplex* d_Self;
    cudaStatus = cudaMalloc((void**)&d_Self, n * n * sizeof(cuComplex));
    assert(cudaStatus == cudaSuccess);
    int* d_info;
    cudaStatus = cudaMalloc((void**)&d_info, sizeof(int));
    assert(cudaStatus == cudaSuccess);

    int lwork;
    cusolverStatus_t cusolverStatus =
        cusolverDnCgetrf_bufferSize(hcusolver, n, n, d_Self, n, &lwork);
    assert(cusolverStatus == cudaSuccess);
    cuComplex* d_work;
    cudaStatus = cudaMalloc((void**)&d_work, lwork * sizeof(cuComplex));
    assert(cudaStatus == cudaSuccess);

    cublasStatus = cublasCcopy(hcublas, n * n, d_Ao, 1, d_Ao_s, 1);
    assert(cublasStatus == cudaSuccess);

    bool okCo = false;
    for(int i1 = 1; i1 <= 300; i1++) {
        *ncyc = i1;

        initKernel<<<num_blocks, BLOCK_SIZE>>>(d_Go, n);

        cublasStatus = cublasCcopy(hcublas, n * n, d_Ao, 1, d_Self, 1);
        assert(cublasStatus == cudaSuccess);

        cusolverStatus = cusolverDnCgetrf(
            hcusolver, n, n, d_Self, n, d_work, d_pivot, d_info
        );
        assert(cusolverStatus == cudaSuccess);
        cusolverStatus = cusolverDnCgetrs(
            hcusolver, CUBLAS_OP_N, n, n, d_Self, n, d_pivot, d_Go, n, d_info
        );
        assert(cusolverStatus == cudaSuccess);

        cublasStatus = cublasCgemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, one, d_Go, n, d_Co, n,
            zero, d_T, n
        );
        assert(cublasStatus == cudaSuccess);

        cublasStatus = cublasCgemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, one, d_Co, n, d_T, n,
            zero, d_C1, n
        );
        assert(cublasStatus == cudaSuccess);

        float summ;
        cublasStatus = cublasScasum(hcublas, n * n, d_C1, 1, &summ);
        assert(cublasStatus == cudaSuccess);
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

        cublasStatus = cublasCgemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, one, d_Bo, n, d_T, n,
            zero, d_Self, n
        );
        assert(cublasStatus == cudaSuccess);

        cublasStatus = cublasCaxpy(hcublas, n * n, mone, d_Self, 1, d_Ao_s, 1);
        assert(cublasStatus == cudaSuccess);
        cublasStatus = cublasCaxpy(hcublas, n * n, mone, d_Self, 1, d_Ao, 1);
        assert(cublasStatus == cudaSuccess);

        cublasStatus = cublasCgemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, one, d_Go, n, d_Bo, n,
            zero, d_T, n
        );
        assert(cublasStatus == cudaSuccess);
        cublasStatus = cublasCgemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, mone, d_Co, n, d_T, n,
            one, d_Ao, n
        );
        assert(cublasStatus == cudaSuccess);

        cublasStatus = cublasCcopy(hcublas, n * n, d_C1, 1, d_Co, 1);
        assert(cublasStatus == cudaSuccess);

        cublasStatus = cublasCgemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, one, d_Bo, n, d_T, n,
            zero, d_C1, n
        );
        assert(cublasStatus == cudaSuccess);

        cublasStatus = cublasCcopy(hcublas, n * n, d_C1, 1, d_Bo, 1);
        assert(cublasStatus == cudaSuccess);
    }

    initKernel<<<num_blocks, BLOCK_SIZE>>>(d_Go, n);
    cublasStatus = cublasCcopy(hcublas, n * n, d_Ao_s, 1, d_Self, 1);
    assert(cublasStatus == cudaSuccess);
    cusolverStatus =
        cusolverDnCgetrf(hcusolver, n, n, d_Self, n, d_work, d_pivot, d_info);
    assert(cusolverStatus == cudaSuccess);
    cusolverStatus = cusolverDnCgetrs(
        hcusolver, CUBLAS_OP_N, n, n, d_Self, n, d_pivot, d_Go, n, d_info
    );
    assert(cusolverStatus == cudaSuccess);

    cudaStatus = cudaMemcpy(
        h_Go_out, d_Go, n * n * sizeof(cuComplex), cudaMemcpyDeviceToHost
    );
    assert(cudaStatus == cudaSuccess);

    cudaStatus = cudaFree(d_pivot);
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_info);
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_Ao);
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_Bo);
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_Co);
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_Go);
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_Ao_s);
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_C1);
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_T);
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_Self);
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_work);
    assert(cudaStatus == cudaSuccess);

    return cudaStatus;
}

extern "C" int cu_Zdecimation(
    cublasHandle_t hcublas, cusolverDnHandle_t hcusolver, void* h_Go_out,
    void* h_Ao_in, void* h_Bo_in, void* h_Co_in, size_t n, int tf32, int* ncyc,
    cuDoubleComplex* one, cuDoubleComplex* mone, cuDoubleComplex* zero,
    double SGFACC
) {
    assert(h_Go_out);
    assert(h_Ao_in);
    assert(h_Bo_in);
    assert(h_Co_in);
    assert(tf32 == 0 || tf32 == 1);
    assert(ncyc);
    assert(one);
    assert(mone);
    assert(zero);
    assert(SGFACC > 0.0);

    auto num_elements = n * n;
    auto num_blocks = (num_elements + BLOCK_SIZE - 1) / BLOCK_SIZE;

    cuDoubleComplex* d_Ao;
    cudaError_t cudaStatus =
        cudaMalloc((void**)&d_Ao, num_elements * sizeof(cuDoubleComplex));
    assert(cudaStatus == cudaSuccess);

    cuDoubleComplex* d_Bo;
    cudaStatus =
        cudaMalloc((void**)&d_Bo, num_elements * sizeof(cuDoubleComplex));
    assert(cudaStatus == cudaSuccess);

    cuDoubleComplex* d_Co;
    cudaStatus =
        cudaMalloc((void**)&d_Co, num_elements * sizeof(cuDoubleComplex));
    assert(cudaStatus == cudaSuccess);

    cudaStatus = cudaMemcpy(
        d_Ao, h_Ao_in, n * n * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice
    );
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaMemcpy(
        d_Bo, h_Bo_in, n * n * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice
    );
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaMemcpy(
        d_Co, h_Co_in, n * n * sizeof(cuDoubleComplex), cudaMemcpyHostToDevice
    );
    assert(cudaStatus == cudaSuccess);

    cublasStatus_t cublasStatus =
        cublasSetPointerMode(hcublas, CUBLAS_POINTER_MODE_HOST);
    assert(cublasStatus == cudaSuccess);

    if(tf32 == 1) {
        cublasStatus = cublasSetMathMode(hcublas, CUBLAS_TENSOR_OP_MATH);
        assert(cublasStatus == cudaSuccess);
    }

    cuDoubleComplex* d_Ao_s;
    cudaStatus = cudaMalloc((void**)&d_Ao_s, n * n * sizeof(cuDoubleComplex));
    assert(cudaStatus == cudaSuccess);
    cuDoubleComplex* d_C1;
    cudaStatus = cudaMalloc((void**)&d_C1, n * n * sizeof(cuDoubleComplex));
    assert(cudaStatus == cudaSuccess);
    cuDoubleComplex* d_Go;
    cudaStatus = cudaMalloc((void**)&d_Go, n * n * sizeof(cuDoubleComplex));
    assert(cudaStatus == cudaSuccess);
    int* d_pivot;
    cudaStatus = cudaMalloc((void**)&d_pivot, n * sizeof(int));
    assert(cudaStatus == cudaSuccess);
    cuDoubleComplex* d_T;
    cudaStatus = cudaMalloc((void**)&d_T, n * n * sizeof(cuDoubleComplex));
    assert(cudaStatus == cudaSuccess);
    cuDoubleComplex* d_Self;
    cudaStatus = cudaMalloc((void**)&d_Self, n * n * sizeof(cuDoubleComplex));
    assert(cudaStatus == cudaSuccess);
    int* d_info;
    cudaStatus = cudaMalloc((void**)&d_info, sizeof(int));
    assert(cudaStatus == cudaSuccess);

    int lwork;
    cusolverStatus_t cusolverStatus =
        cusolverDnZgetrf_bufferSize(hcusolver, n, n, d_Self, n, &lwork);
    assert(cusolverStatus == cudaSuccess);
    cuDoubleComplex* d_work;
    cudaStatus = cudaMalloc((void**)&d_work, lwork * sizeof(cuDoubleComplex));
    assert(cudaStatus == cudaSuccess);

    cublasStatus = cublasZcopy(hcublas, n * n, d_Ao, 1, d_Ao_s, 1);
    assert(cublasStatus == cudaSuccess);

    bool okCo = false;
    for(int i1 = 1; i1 <= 300; i1++) {
        *ncyc = i1;

        initKernel<<<num_blocks, BLOCK_SIZE>>>(d_Go, n);

        cublasStatus = cublasZcopy(hcublas, n * n, d_Ao, 1, d_Self, 1);
        assert(cublasStatus == cudaSuccess);

        cusolverStatus = cusolverDnZgetrf(
            hcusolver, n, n, d_Self, n, d_work, d_pivot, d_info
        );
        assert(cusolverStatus == cudaSuccess);
        cusolverStatus = cusolverDnZgetrs(
            hcusolver, CUBLAS_OP_N, n, n, d_Self, n, d_pivot, d_Go, n, d_info
        );
        assert(cusolverStatus == cudaSuccess);

        cublasStatus = cublasZgemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, one, d_Go, n, d_Co, n,
            zero, d_T, n
        );
        assert(cublasStatus == cudaSuccess);

        cublasStatus = cublasZgemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, one, d_Co, n, d_T, n,
            zero, d_C1, n
        );
        assert(cublasStatus == cudaSuccess);

        double summ;
        cublasStatus = cublasDzasum(hcublas, n * n, d_C1, 1, &summ);
        assert(cublasStatus == cudaSuccess);
        //printf("loop it= %d , summ= %f \n ", i1, summ);

        if(summ <= SGFACC) {
            if(okCo) {
                break;
            } else {
                okCo = true;
            }
        } else {
            okCo = false;
        }

        cublasStatus = cublasZgemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, one, d_Bo, n, d_T, n,
            zero, d_Self, n
        );
        assert(cublasStatus == cudaSuccess);

        cublasStatus = cublasZaxpy(hcublas, n * n, mone, d_Self, 1, d_Ao_s, 1);
        assert(cublasStatus == cudaSuccess);
        cublasStatus = cublasZaxpy(hcublas, n * n, mone, d_Self, 1, d_Ao, 1);
        assert(cublasStatus == cudaSuccess);

        cublasStatus = cublasZgemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, one, d_Go, n, d_Bo, n,
            zero, d_T, n
        );
        assert(cublasStatus == cudaSuccess);
        cublasStatus = cublasZgemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, mone, d_Co, n, d_T, n,
            one, d_Ao, n
        );
        assert(cublasStatus == cudaSuccess);

        cublasStatus = cublasZcopy(hcublas, n * n, d_C1, 1, d_Co, 1);
        assert(cublasStatus == cudaSuccess);

        cublasStatus = cublasZgemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, one, d_Bo, n, d_T, n,
            zero, d_C1, n
        );
        assert(cublasStatus == cudaSuccess);

        cublasStatus = cublasZcopy(hcublas, n * n, d_C1, 1, d_Bo, 1);
        assert(cublasStatus == cudaSuccess);
    }

    initKernel<<<num_blocks, BLOCK_SIZE>>>(d_Go, n);
    cublasStatus = cublasZcopy(hcublas, n * n, d_Ao_s, 1, d_Self, 1);
    assert(cublasStatus == cudaSuccess);
    cusolverStatus =
        cusolverDnZgetrf(hcusolver, n, n, d_Self, n, d_work, d_pivot, d_info);
    assert(cusolverStatus == cudaSuccess);
    cusolverStatus = cusolverDnZgetrs(
        hcusolver, CUBLAS_OP_N, n, n, d_Self, n, d_pivot, d_Go, n, d_info
    );
    assert(cusolverStatus == cudaSuccess);

    cudaStatus = cudaMemcpy(
        h_Go_out, d_Go, n * n * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost
    );
    assert(cudaStatus == cudaSuccess);

    cudaStatus = cudaFree(d_pivot);
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_info);
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_Ao);
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_Bo);
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_Co);
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_Go);
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_Ao_s);
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_C1);
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_T);
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_Self);
    assert(cudaStatus == cudaSuccess);
    cudaStatus = cudaFree(d_work);
    assert(cudaStatus == cudaSuccess);

    return cudaStatus;
}

extern "C" int cu_meminfo(size_t* freemem, size_t* totalmem) {
    cudaError_t cudaStatus;
    cudaStatus = cudaDeviceSynchronize();
    cudaStatus = cudaMemGetInfo(freemem, totalmem);
    assert(cudaStatus == cudaSuccess);
    return cudaStatus;
}
