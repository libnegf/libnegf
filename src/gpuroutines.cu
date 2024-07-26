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

#include "cublas_v2.h"
#include "cusolverDn.h"
#include <assert.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>

#define BLOCK_SIZE 1024
#define TILE_DIM 32
#define BLOCK_ROWS 8

__global__ void CaddKernel(
    cuComplex *c, const cuComplex alpha, const cuComplex *a,
    const cuComplex beta, const cuComplex *b, int size
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if(i < size) {
        c[i].x = (alpha.x * a[i].x - alpha.y * a[i].y) +
                 (beta.x * b[i].x - beta.y * b[i].y);
        c[i].y = (alpha.x * a[i].y + alpha.y * a[i].x) +
                 (beta.x * b[i].y + beta.y * b[i].x);
    }
}

__global__ void ZaddKernel(
    cuDoubleComplex *c, const cuDoubleComplex alpha, const cuDoubleComplex *a,
    const cuDoubleComplex beta, const cuDoubleComplex *b, int size
) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
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
__global__ void CinitKernel(cuComplex *a, int nrow) {
    int size;
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    size = nrow * nrow;
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

__global__ void ZinitKernel(cuDoubleComplex *a, int nrow) {
    int size;
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    size = nrow * nrow;
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

__global__ void DinitKernel(double *a, int nrow) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if(i < nrow) {
        a[i] = 1.0;
    }
}

__global__ void SinitKernel(float *a, int nrow) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    if(i < nrow) {
        a[i] = 1.0;
    }
}

__global__ void CtraceKernel(
    cuComplex *a, int nrow, float *trace, bool *mask, int mask_present
) {
    int size;
    int i = blockDim.x * blockIdx.x + threadIdx.x;

    size = nrow * nrow;
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
                }
            }
        }
    }
}

__global__ void ZtraceKernel(
    cuDoubleComplex *a, int nrow, double *trace, bool *mask, int mask_present
) {
    int i = blockDim.x * blockIdx.x + threadIdx.x;
    int size = nrow * nrow;
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
                }
            }
        }
    }
}

/*~-~-~-~-~-~-~-~-~-~-~-~-~-~ DATA MOVEMENT  ROUTINES
 * -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-*/

extern "C" int cu_createMat(void **d_A, int bytecount) {
    cudaError_t err;
    err = cudaMalloc(d_A, bytecount);
    // printf("GPU Address: %p \n",*d_A);
    return err;
}

extern "C" int cu_copyMatH2D(void *h_A, void *d_A, int bytecount) {
    cudaError_t err;
    // printf("copy %p to %p\n",h_A,d_A);
    err = cudaMemcpy(d_A, h_A, bytecount, cudaMemcpyHostToDevice);
    return err;
}

extern "C" int cu_copyMatD2H(void *h_A, void *d_A, int bytecount) {
    cudaError_t err;
    // printf("copy %p to %p\n",d_A,h_A);
    err = cudaMemcpy(h_A, d_A, bytecount, cudaMemcpyDeviceToHost);
    return err;
}

extern "C" int cu_deleteMat(void *d_A) {
    cudaError_t err;
    // printf("add_free: %p",d_A);
    err = cudaFree(d_A);
    return err;
}

/*~-~-~-~-~-~-~-~-~-~-~-~-~-~ INIT/FINAL ROUTINES
 * -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-*/

extern "C" int cu_cudaGetDeviceCount(int *count) {
    cudaError_t err;

    err = cudaGetDeviceCount(count);
    return err;
}

extern "C" int cu_cudaGetDeviceProperties(int device) {
    cudaError_t err;
    cudaDeviceProp prop;

    err = cudaGetDeviceProperties(&prop, device);

    printf(" Found GPU: Device Name: %s\n", prop.name);
    printf(" TotalMemory: %lu\n", (unsigned long)prop.totalGlobalMem);
    printf(" Shared per block: %lu\n", (unsigned long)prop.sharedMemPerBlock);
    return err;
}

extern "C" int cu_cublasInit(cublasHandle_t *hcublas) {
    cublasStatus_t err;
    err = cublasCreate(hcublas);
    if(err != 0) {
        printf("cublas create error: %d\n", err);
    }
    // printf("hcublas Addr: %p \n",*hcublas);
    return err;
}

extern "C" int cu_cublasFinalize(cublasHandle_t hcublas) {
    cublasStatus_t err;
    err = cublasDestroy(hcublas);
    return err;
}

extern "C" int cu_cusolverInit(cusolverDnHandle_t *hcusolver) {
    cusolverStatus_t err;
    err = cusolverDnCreate(hcusolver);
    if(err != 0) {
        printf("cusolver create error: %d\n", err);
    }
    // printf("hcusolver Addr: %p \n",*hcusolver);
    return err;
}

extern "C" int cu_cusolverFinalize(cusolverDnHandle_t hcusolver) {
    cusolverStatus_t err;
    err = cusolverDnDestroy(hcusolver);
    return err;
}

/*~-~-~-~-~-~-~-~-~-~-~-~-~-~ MATRIX ROUTINES
 * -~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-~-*/
/* C = alpha op(A) op(B) + beta C
 * m: #rows of op(A)
 * n: #cols of op(B)
 * k: #cols of op(A) = #rows of op(B)
 */
extern "C" int cu_CmultMat(
    cublasHandle_t hcublas, int m, int n, int k, cuComplex *alpha, void *d_A,
    void *d_B, cuComplex *beta, void *d_C, int dagger
) {
    cuComplex *pdA, *pdB, *pdC;

    // printf("A: %p B: %p C: %p\n",d_A,d_B,d_C);
    pdA = (cuComplex *)d_A;
    pdB = (cuComplex *)d_B;
    pdC = (cuComplex *)d_C;
    cublasStatus_t err;
    if(dagger == 0) {
        err = cublasCgemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, alpha, pdA, m, pdB, k,
            beta, pdC, m
        );
    }
    if(dagger == 1) {
        err = cublasCgemm(
            hcublas, CUBLAS_OP_C, CUBLAS_OP_N, m, n, k, alpha, pdA, k, pdB, k,
            beta, pdC, m
        );
    }
    if(dagger == 2) {
        err = cublasCgemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_C, m, n, k, alpha, pdA, m, pdB, n,
            beta, pdC, m
        );
    }
    return err;
}

//  C = alpha op(A) op(B) + beta C
// op(A):  m x k
// op(B):  k x n
//     C:  m x n
extern "C" int cu_ZmultMat(
    cublasHandle_t hcublas, int m, int n, int k, cuDoubleComplex *alpha,
    void *d_A, void *d_B, cuDoubleComplex *beta, void *d_C, int dagger
) {
    cuDoubleComplex *pdA, *pdB, *pdC;

    pdA = (cuDoubleComplex *)d_A;
    pdB = (cuDoubleComplex *)d_B;
    pdC = (cuDoubleComplex *)d_C;
    cublasStatus_t err;
    if(dagger == 0) {
        err = cublasZgemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, m, n, k, alpha, pdA, m, pdB, k,
            beta, pdC, m
        );
    }
    if(dagger == 1) {
        err = cublasZgemm(
            hcublas, CUBLAS_OP_C, CUBLAS_OP_N, m, n, k, alpha, pdA, k, pdB, k,
            beta, pdC, m
        );
    }
    if(dagger == 2) {
        err = cublasZgemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_C, m, n, k, alpha, pdA, m, pdB, n,
            beta, pdC, m
        );
    }
    return err;
}

extern "C" int cu_Cinverse(
    cublasHandle_t hcublas, cusolverDnHandle_t hcusolver, void *d_A,
    void *d_Ainv, int N
) {
    cudaError_t cudaStatus;
    cusolverStatus_t cusolverStatus;
    cublasStatus_t cublasStatus;
    // declare arrays on the device
    cuComplex *pdA, *pdAinv, *d_LU, *d_Work;

    pdA = (cuComplex *)d_A;
    pdAinv = (cuComplex *)d_Ainv;
    // coeff . matrix , rhs , workspace
    int *d_pivot, *d_info, lwork; // pivots , info , worksp . size
    int info_gpu = 0;

    // compute buffer size and prep . memory
    cusolverStatus =
        cusolverDnCgetrf_bufferSize(hcusolver, N, N, pdA, N, &lwork);
    // prepare memory on the device

    cudaStatus = cudaMalloc((void **)&d_LU, N * N * sizeof(cuComplex));
    cudaStatus = cudaMalloc((void **)&d_pivot, N * sizeof(int));
    cudaStatus = cudaMalloc((void **)&d_info, sizeof(int));
    // copy d_LU <- pdA
    cublasStatus = cublasCcopy(hcublas, N * N, pdA, 1, d_LU, 1);

    cudaStatus = cudaMalloc((void **)&d_Work, lwork * sizeof(cuComplex));

    // LU factorization of d_A , with partial pivoting and row
    // interchanges ; row i is interchanged with row d_pivot ( i );
    cusolverStatus =
        cusolverDnCgetrf(hcusolver, N, N, d_LU, N, d_Work, d_pivot, d_info);

    // use the LU factorization to solve the system d_LU * x = d_Ainv ;
    // the solution overwrites d_Ainv
    cusolverStatus = cusolverDnCgetrs(
        hcusolver, CUBLAS_OP_N, N, N, d_LU, N, d_pivot, pdAinv, N, d_info
    );

    cudaStatus =
        cudaMemcpy(&info_gpu, d_info, sizeof(int), cudaMemcpyDeviceToHost);
    // d_info -> info_gpu
    cudaStatus = cudaFree(d_pivot);
    cudaStatus = cudaFree(d_info);
    cudaStatus = cudaFree(d_Work);
    cudaStatus = cudaFree(d_LU);
    return cudaStatus;
}

extern "C" int cu_Zinverse(
    cublasHandle_t hcublas, cusolverDnHandle_t hcusolver, void *d_A,
    void *d_Ainv, int N
) {
    cudaError_t cudaStatus;
    cusolverStatus_t cusolverStatus;
    cublasStatus_t cublasStatus;
    // declare arrays on the device
    cuDoubleComplex *pdA, *pdAinv, *d_LU, *d_Work;

    pdA = (cuDoubleComplex *)d_A;
    pdAinv = (cuDoubleComplex *)d_Ainv;
    // coeff . matrix , rhs , workspace
    int *d_pivot, *d_info, lwork; // pivots , info , worksp . size
    int info_gpu = 0;

    // compute buffer size and prep . memory
    cusolverStatus =
        cusolverDnZgetrf_bufferSize(hcusolver, N, N, pdA, N, &lwork);
    // prepare memory on the device

    cudaStatus = cudaMalloc((void **)&d_LU, N * N * sizeof(cuDoubleComplex));
    cudaStatus = cudaMalloc((void **)&d_pivot, N * sizeof(int));
    cudaStatus = cudaMalloc((void **)&d_info, sizeof(int));
    // copy d_LU <- pdA
    cublasStatus = cublasZcopy(hcublas, N * N, pdA, 1, d_LU, 1);

    cudaStatus = cudaMalloc((void **)&d_Work, lwork * sizeof(cuDoubleComplex));

    // LU factorization of d_A , with partial pivoting and row
    // interchanges ; row i is interchanged with row d_pivot ( i );
    cusolverStatus =
        cusolverDnZgetrf(hcusolver, N, N, d_LU, N, d_Work, d_pivot, d_info);

    // use the LU factorization to solve the system d_LU * x = d_Ainv ;
    // the solution overwrites d_Ainv
    cusolverStatus = cusolverDnZgetrs(
        hcusolver, CUBLAS_OP_N, N, N, d_LU, N, d_pivot, pdAinv, N, d_info
    );

    cudaStatus =
        cudaMemcpy(&info_gpu, d_info, sizeof(int), cudaMemcpyDeviceToHost);
    // d_info -> info_gpu
    cudaStatus = cudaFree(d_pivot);
    cudaStatus = cudaFree(d_info);
    cudaStatus = cudaFree(d_Work);
    cudaStatus = cudaFree(d_LU);
    return cudaStatus;
}

extern "C" int cu_Ckernelsum(
    void *d_C, cuComplex *alpha, void *d_A, cuComplex *beta, void *d_B, int size
) {
    int num_blocks;
    cuComplex *pdA = (cuComplex *)d_A;
    cuComplex *pdB = (cuComplex *)d_B;
    cuComplex *pdC = (cuComplex *)d_C;

    num_blocks = (size / BLOCK_SIZE) + 1;

    CaddKernel<<<num_blocks, BLOCK_SIZE>>>(pdC, *alpha, pdA, *beta, pdB, size);

    return 0;
}

extern "C" int cu_Zkernelsum(
    void *d_C, cuDoubleComplex *alpha, void *d_A, cuDoubleComplex *beta,
    void *d_B, int size
) {
    int num_blocks;
    cuDoubleComplex *pdA = (cuDoubleComplex *)d_A;
    cuDoubleComplex *pdB = (cuDoubleComplex *)d_B;
    cuDoubleComplex *pdC = (cuDoubleComplex *)d_C;

    num_blocks = (size / BLOCK_SIZE) + 1;

    ZaddKernel<<<num_blocks, BLOCK_SIZE>>>(pdC, *alpha, pdA, *beta, pdB, size);

    return 0;
}

extern "C" int cu_Cmatsum(
    cublasHandle_t hcublas, int m, int n, cuComplex *alpha, void *d_A,
    cuComplex *beta, void *d_B, void *d_C, int dagger
) {
    // m number of rows of matrix op(A) and C
    // n number of columns of matrix op(B) and C
    cuComplex *pdA = (cuComplex *)d_A;
    cuComplex *pdB = (cuComplex *)d_B;
    cuComplex *pdC = (cuComplex *)d_C;

    cublasStatus_t err;
    if(dagger == 0) {
        err = cublasCgeam(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, m, n, alpha, pdA, m, beta, pdB,
            m, pdC, m
        );
    }
    if(dagger == 1) {
        err = cublasCgeam(
            hcublas, CUBLAS_OP_C, CUBLAS_OP_N, m, n, alpha, pdA, n, beta, pdB,
            m, pdC, m
        );
    }
    if(dagger == 2) {
        err = cublasCgeam(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_C, m, n, alpha, pdA, m, beta, pdB,
            n, pdC, m
        );
    }
    return err;
}

extern "C" int cu_Zmatsum(
    cublasHandle_t hcublas, int m, int n, cuDoubleComplex *alpha, void *d_A,
    cuDoubleComplex *beta, void *d_B, void *d_C, int dagger
) {
    // m number of rows of matrix op(A) and C
    // n number of columns of matrix op(B) and C
    cuDoubleComplex *pdA = (cuDoubleComplex *)d_A;
    cuDoubleComplex *pdB = (cuDoubleComplex *)d_B;
    cuDoubleComplex *pdC = (cuDoubleComplex *)d_C;

    cublasStatus_t err;
    if(dagger == 0) {
        err = cublasZgeam(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, m, n, alpha, pdA, m, beta, pdB,
            m, pdC, m
        );
    }
    if(dagger == 1) {
        err = cublasZgeam(
            hcublas, CUBLAS_OP_C, CUBLAS_OP_N, m, n, alpha, pdA, n, beta, pdB,
            m, pdC, m
        );
    }
    if(dagger == 2) {
        err = cublasZgeam(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_C, m, n, alpha, pdA, m, beta, pdB,
            n, pdC, m
        );
    }
    return err;
}

extern "C" int cu_Cinitmat(void *d_A, int nrow) {
    int num_blocks;
    int size = nrow * nrow;
    cuComplex *pdA = (cuComplex *)d_A;

    num_blocks = (size / BLOCK_SIZE) + 1;

    CinitKernel<<<num_blocks, BLOCK_SIZE>>>(pdA, nrow);

    return 0;
}

extern "C" int cu_Zinitmat(void *d_A, int nrow) {
    int num_blocks;
    int size = nrow * nrow;
    cuDoubleComplex *pdA = (cuDoubleComplex *)d_A;

    num_blocks = (size / BLOCK_SIZE) + 1;

    ZinitKernel<<<num_blocks, BLOCK_SIZE>>>(pdA, nrow);

    return 0;
}

extern "C" float cu_Ctrace(
    cublasHandle_t hcublas, void *d_A, int nrow, void *h_mask, int mask_present
) {
    cuComplex *pdA = (cuComplex *)d_A;
    int size = nrow * nrow;
    int num_blocks = (size / BLOCK_SIZE) + 1;
    float *d_work;
    cudaError_t cudaStatus = cudaMalloc((void **)&d_work, nrow * sizeof(float));
    float *d_iden;
    cudaStatus = cudaMalloc((void **)&d_iden, nrow * sizeof(float));
    assert(cudaStatus == cudaSuccess);
    bool *d_mask;
    cudaStatus = cudaMalloc((void **)&d_mask, nrow * sizeof(bool));
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
    assert(err == cudaSuccess);

    cudaStatus = cudaFree(d_work);
    cudaStatus = cudaFree(d_iden);
    cudaStatus = cudaFree(d_mask);

    return result;
}

extern "C" double cu_Ztrace(
    cublasHandle_t hcublas, void *d_A, int nrow, void *h_mask, int mask_present
) {
    cuDoubleComplex *pdA = (cuDoubleComplex *)d_A;
    int size = nrow * nrow;
    int num_blocks = (size / BLOCK_SIZE) + 1;
    double *d_work;
    cudaError_t cudaStatus =
        cudaMalloc((void **)&d_work, nrow * sizeof(double));
    double *d_iden;
    cudaStatus = cudaMalloc((void **)&d_iden, nrow * sizeof(double));
    assert(cudaStatus == cudaSuccess);
    bool *d_mask;
    cudaStatus = cudaMalloc((void **)&d_mask, nrow * sizeof(bool));
    assert(cudaStatus == cudaSuccess);
    if(h_mask) {
        cudaStatus = cudaMemcpy(
            d_mask, h_mask, nrow * sizeof(bool), cudaMemcpyHostToDevice
        );
        assert(cudaStatus == cudaSuccess);
    }

    DinitKernel<<<num_blocks, BLOCK_SIZE>>>(d_iden, nrow);
    ZtraceKernel<<<num_blocks, BLOCK_SIZE>>>(
        pdA, nrow, d_work, d_mask, mask_present
    );

    double result;
    cublasStatus_t err =
        cublasDdot(hcublas, nrow, d_iden, 1, d_work, 1, &result);
    assert(err == cudaSuccess);

    cudaStatus = cudaFree(d_work);
    cudaStatus = cudaFree(d_iden);
    cudaStatus = cudaFree(d_mask);

    return result;
}

extern "C" int
cu_Cmatcopy(cublasHandle_t hcublas, void *d_A, void *d_B, int N) {
    cuComplex *pdA = (cuComplex *)d_A;
    cuComplex *pdB = (cuComplex *)d_B;

    cublasStatus_t err;

    err = cublasCcopy(hcublas, N * N, pdA, 1, pdB, 1);
    return err;
}

extern "C" int
cu_Zmatcopy(cublasHandle_t hcublas, void *d_A, void *d_B, int size) {
    cuDoubleComplex *pdA = (cuDoubleComplex *)d_A;
    cuDoubleComplex *pdB = (cuDoubleComplex *)d_B;

    cublasStatus_t err;

    err = cublasZcopy(hcublas, size, pdA, 1, pdB, 1);
    return err;
}

extern "C" int cu_Casum(cublasHandle_t hcublas, void *d_A, float *summ, int N) {
    cuComplex *pdA = (cuComplex *)d_A;

    cublasStatus_t err;

    err = cublasScasum(hcublas, N, pdA, 1, summ);
    return err;
}

extern "C" int
cu_Zasum(cublasHandle_t hcublas, void *d_A, double *summ, int N) {
    cuDoubleComplex *pdA = (cuDoubleComplex *)d_A;

    cublasStatus_t err;

    err = cublasDzasum(hcublas, N, pdA, 1, summ);
    return err;
}

extern "C" int cu_Cdecimation(
    cublasHandle_t hcublas, cusolverDnHandle_t hcusolver, void *h_Go_out,
    void *h_Ao_in, void *h_Bo_in, void *h_Co_in, int n, int tf32, int *ncyc,
    cuComplex *one, cuComplex *mone, cuComplex *zero, float SGFACC
) {
    cudaError_t cudaStatus;
    cusolverStatus_t cusolverStatus;
    cublasStatus_t cublasStatus;

    float summ;
    bool okCo = false;
    cuComplex *phGo_out = (cuComplex *)h_Go_out;
    cuComplex *phAo_in = (cuComplex *)h_Ao_in;
    cuComplex *phBo_in = (cuComplex *)h_Bo_in;
    cuComplex *phCo_in = (cuComplex *)h_Co_in;

    cuComplex *d_Ao, *d_Bo, *d_Co, *d_Go, *d_Ao_s, *d_C1, *d_T, *d_Self,
        *d_work;
    int *d_pivot, *d_info, lwork; // pivots , info , worksp . size
    int i1, num_blocks;

    num_blocks = ((n * n) / BLOCK_SIZE) + 1;

    cudaStatus = cudaMalloc((void **)&d_Ao, n * n * sizeof(cuComplex));
    cudaStatus = cudaMalloc((void **)&d_Bo, n * n * sizeof(cuComplex));
    cudaStatus = cudaMalloc((void **)&d_Co, n * n * sizeof(cuComplex));
    cudaStatus = cudaMemcpy(
        d_Ao, phAo_in, n * n * sizeof(cuComplex), cudaMemcpyHostToDevice
    );
    cudaStatus = cudaMemcpy(
        d_Bo, phBo_in, n * n * sizeof(cuComplex), cudaMemcpyHostToDevice
    );
    cudaStatus = cudaMemcpy(
        d_Co, phCo_in, n * n * sizeof(cuComplex), cudaMemcpyHostToDevice
    );

    cublasStatus = cublasSetPointerMode(hcublas, CUBLAS_POINTER_MODE_HOST);
    if(tf32 == 1) {
        cublasStatus = cublasSetMathMode(hcublas, CUBLAS_TENSOR_OP_MATH);
    }

    cudaStatus = cudaMalloc((void **)&d_Ao_s, n * n * sizeof(cuComplex));
    cudaStatus = cudaMalloc((void **)&d_C1, n * n * sizeof(cuComplex));
    cudaStatus = cudaMalloc((void **)&d_Go, n * n * sizeof(cuComplex));
    cudaStatus = cudaMalloc((void **)&d_pivot, n * sizeof(int));
    cudaStatus = cudaMalloc((void **)&d_T, n * n * sizeof(cuComplex));
    cudaStatus = cudaMalloc((void **)&d_Self, n * n * sizeof(cuComplex));
    cudaStatus = cudaMalloc((void **)&d_info, sizeof(int));

    cusolverStatus =
        cusolverDnCgetrf_bufferSize(hcusolver, n, n, d_Self, n, &lwork);
    cudaStatus = cudaMalloc((void **)&d_work, lwork * sizeof(cuComplex));

    cublasStatus = cublasCcopy(hcublas, n * n, d_Ao, 1, d_Ao_s, 1);

    for(i1 = 1; i1 <= 300; i1++) {
        *ncyc = i1;

        CinitKernel<<<num_blocks, BLOCK_SIZE>>>(d_Go, n);

        cublasStatus = cublasCcopy(hcublas, n * n, d_Ao, 1, d_Self, 1);

        cusolverStatus = cusolverDnCgetrf(
            hcusolver, n, n, d_Self, n, d_work, d_pivot, d_info
        );
        cusolverStatus = cusolverDnCgetrs(
            hcusolver, CUBLAS_OP_N, n, n, d_Self, n, d_pivot, d_Go, n, d_info
        );

        cublasStatus = cublasCgemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, one, d_Go, n, d_Co, n,
            zero, d_T, n
        );

        cublasStatus = cublasCgemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, one, d_Co, n, d_T, n,
            zero, d_C1, n
        );

        cublasStatus = cublasScasum(hcublas, n * n, d_C1, 1, &summ);
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

        cublasStatus = cublasCaxpy(hcublas, n * n, mone, d_Self, 1, d_Ao_s, 1);
        cublasStatus = cublasCaxpy(hcublas, n * n, mone, d_Self, 1, d_Ao, 1);

        cublasStatus = cublasCgemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, one, d_Go, n, d_Bo, n,
            zero, d_T, n
        );
        cublasStatus = cublasCgemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, mone, d_Co, n, d_T, n,
            one, d_Ao, n
        );

        cublasStatus = cublasCcopy(hcublas, n * n, d_C1, 1, d_Co, 1);

        cublasStatus = cublasCgemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, one, d_Bo, n, d_T, n,
            zero, d_C1, n
        );

        cublasStatus = cublasCcopy(hcublas, n * n, d_C1, 1, d_Bo, 1);
    }

    CinitKernel<<<num_blocks, BLOCK_SIZE>>>(d_Go, n);
    cublasStatus = cublasCcopy(hcublas, n * n, d_Ao_s, 1, d_Self, 1);
    cusolverStatus =
        cusolverDnCgetrf(hcusolver, n, n, d_Self, n, d_work, d_pivot, d_info);
    cusolverStatus = cusolverDnCgetrs(
        hcusolver, CUBLAS_OP_N, n, n, d_Self, n, d_pivot, d_Go, n, d_info
    );

    // cublasStatus = cublasCcopy(hcublas, n*n, d_Go, 1, d_Go_out, 1);

    cudaStatus = cudaMemcpy(
        phGo_out, d_Go, n * n * sizeof(cuComplex), cudaMemcpyDeviceToHost
    );

    cudaStatus = cudaFree(d_Ao);
    cudaStatus = cudaFree(d_Bo);
    cudaStatus = cudaFree(d_Co);
    cudaStatus = cudaFree(d_Go);
    cudaStatus = cudaFree(d_Ao_s);
    cudaStatus = cudaFree(d_C1);
    cudaStatus = cudaFree(d_pivot);
    cudaStatus = cudaFree(d_T);
    cudaStatus = cudaFree(d_Self);
    cudaStatus = cudaFree(d_info);
    cudaStatus = cudaFree(d_work);

    return cudaStatus;
}

extern "C" int cu_Zdecimation(
    cublasHandle_t hcublas, cusolverDnHandle_t hcusolver, void *h_Go_out,
    void *h_Ao_in, void *h_Bo_in, void *h_Co_in, int n, int tf32, int *ncyc,
    cuDoubleComplex *one, cuDoubleComplex *mone, cuDoubleComplex *zero,
    double SGFACC
) {
    assert(h_Go_out);
    assert(h_Ao_in);
    assert(h_Bo_in);
    assert(h_Co_in);
    assert(n > 0);
    assert(tf32 == 0 || tf32 == 1);
    assert(ncyc);
    assert(one);
    assert(mone);
    assert(zero);
    assert(SGFACC > 0.0);

    int num_elements = n * n;
    int num_blocks = (num_elements / BLOCK_SIZE) + 1;

    cuDoubleComplex *d_Ao;
    cudaError_t cudaStatus =
        cudaMalloc((void **)&d_Ao, num_elements * sizeof(cuDoubleComplex));
    assert(cudaStatus == cudaSuccess);

    cuDoubleComplex *d_Bo;
    cudaStatus =
        cudaMalloc((void **)&d_Bo, num_elements * sizeof(cuDoubleComplex));
    assert(cudaStatus == cudaSuccess);

    cuDoubleComplex *d_Co;
    cudaStatus =
        cudaMalloc((void **)&d_Co, num_elements * sizeof(cuDoubleComplex));
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

    cuDoubleComplex *d_Ao_s;
    cudaStatus = cudaMalloc((void **)&d_Ao_s, n * n * sizeof(cuDoubleComplex));
    assert(cudaStatus == cudaSuccess);
    cuDoubleComplex *d_C1;
    cudaStatus = cudaMalloc((void **)&d_C1, n * n * sizeof(cuDoubleComplex));
    assert(cudaStatus == cudaSuccess);
    cuDoubleComplex *d_Go;
    cudaStatus = cudaMalloc((void **)&d_Go, n * n * sizeof(cuDoubleComplex));
    assert(cudaStatus == cudaSuccess);
    int *d_pivot;
    cudaStatus = cudaMalloc((void **)&d_pivot, n * sizeof(int));
    assert(cudaStatus == cudaSuccess);
    cuDoubleComplex *d_T;
    cudaStatus = cudaMalloc((void **)&d_T, n * n * sizeof(cuDoubleComplex));
    assert(cudaStatus == cudaSuccess);
    cuDoubleComplex *d_Self;
    cudaStatus = cudaMalloc((void **)&d_Self, n * n * sizeof(cuDoubleComplex));
    assert(cudaStatus == cudaSuccess);
    int *d_info;
    cudaStatus = cudaMalloc((void **)&d_info, sizeof(int));
    assert(cudaStatus == cudaSuccess);

    int lwork;
    cusolverStatus_t cusolverStatus =
        cusolverDnZgetrf_bufferSize(hcusolver, n, n, d_Self, n, &lwork);
    assert(cusolverStatus == cudaSuccess);
    cuDoubleComplex *d_work;
    cudaStatus = cudaMalloc((void **)&d_work, lwork * sizeof(cuDoubleComplex));
    assert(cudaStatus == cudaSuccess);

    cublasStatus = cublasZcopy(hcublas, n * n, d_Ao, 1, d_Ao_s, 1);
    assert(cublasStatus == cudaSuccess);

    bool okCo = false;
    for(int i1 = 1; i1 <= 300; i1++) {
        *ncyc = i1;

        ZinitKernel<<<num_blocks, BLOCK_SIZE>>>(d_Go, n);

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

    ZinitKernel<<<num_blocks, BLOCK_SIZE>>>(d_Go, n);
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
