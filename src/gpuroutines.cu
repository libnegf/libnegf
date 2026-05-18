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
#include <cub/cub.cuh>

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
#define CUDA_CHECK(error) \
    do { \
        if(!(error == cudaSuccess)) { \
            std::fprintf( \
                stderr, "CUDA error in %s line %u \n", __FILE__, \
                __LINE__ \
		    ); \
            std::fprintf( \
                stderr, "ERROR \n %s \n", cudaGetErrorString(error) \
		    ); \
            std::exit(EXIT_FAILURE); \
        } \
    } while(false)


#define LIBNEGF_CHECK_INFO(info_var) \
    do { \
        check_info<<<1, 1>>>(__FILE__, __LINE__, info_var); \
    } while(false)


/// Check if the preceeding BLAS-style routine terminated successfully.
__global__ void check_info(const char* file, unsigned line, int* d_info) {
    if(!d_info) {
        std::printf("%s:%u: info parameter is NULL\n", file, line);
        __trap();
    }
    if(d_info && *d_info != 0) {
        std::printf(
            "%s:%u: got info=%d, expected info=0\n", file, line, *d_info
        );
        __trap();
    }
}


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
    size_t num_elements
) {
    ENFORCE(c);
    ENFORCE(a);
    ENFORCE(b);

    size_t i = blockIdx.x * blockDim.x + threadIdx.x;

    if(i < num_elements) {
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
__global__ void initKernel(Number* a, int32_t nrow) {
    using Real = typename get_real<Number>::type;

    ENFORCE(a);

    size_t n64 = static_cast<size_t>(nrow);
    size_t num_elements = n64 * n64;
    size_t i = blockDim.x * blockIdx.x + threadIdx.x;
    auto one = Number{Real{1}};
    auto zero = Number{Real{0}};

    if(i < num_elements) {
        if(i % (nrow + 1) == 0) {
            a[i] = one;
        } else {
            a[i] = zero;
        }
    }
}


template<typename Number>
__global__ void initArrayWithOnes(Number* a, int32_t nrow) {
    static_assert(std::is_same<Number, typename get_real<Number>::type>::value);

    assert(a);

    auto i = blockDim.x * blockIdx.x + threadIdx.x;

    if(i < nrow) {
        a[i] = Number{1};
    }
}

/* ***************************************************************
 * Extract diagonal from a square matrix
 * ***************************************************************
 */

template<typename Number, typename Real = typename get_real<Number>::type>
__global__ void getdiagKernel(const Number* A, int n, Real* out,
                              const bool* mask, int mask_present)
{
    ENFORCE(A);
    ENFORCE(out);

    int i = blockIdx.x * blockDim.x + threadIdx.x;

    if (i < n) {
        if (mask_present == 0 || mask[i]) {
            out[i] = A[i * n + i].x;
        } else {
            out[i] = 0;
        }
    }
}

/* ***************************************************************
 * Simple trace kernel with a mask
 * ***************************************************************
 */
/*
template<typename Number>
__global__ void traceKernel(const Number* A, int32_t nrow, Number* result,
                           bool* mask, int mask_present)
{
    ENFORCE(A);
    ENFORCE(result);
    ENFORCE(mask || mask_present == 0);
    ENFORCE(mask_present == 0 || mask_present == 1);

    size_t i = blockIdx.x * blockDim.x + threadIdx.x;

    if(mask_present == 0) {
      if (i < nrow) {
        Number val = A[i * nrow + i];
        atomicAdd(&((*result).x), val.x);
        atomicAdd(&((*result).y), val.y);
      }
    }
    else
    {
      if (i < nrow && mask[i]) {
        Number val = A[i * nrow + i];
        atomicAdd(&((*result).x), val.x);
        atomicAdd(&((*result).y), val.y);
      }
    }
}
*/

__device__ inline float absval(cuComplex z)  { return cuCabsf(z); }
__device__ inline double absval(cuDoubleComplex z){ return cuCabs(z); }

__device__ __forceinline__ float atomicMaxAbsFloat(float* addr, float value) {
    // assumes value >= 0
    return __int_as_float(atomicMax((int*)addr, __float_as_int(value)));
}

__device__ __forceinline__ double atomicMaxAbsDouble(double* addr, double value)
{
    // assumes value >= 0
    return __longlong_as_double(
        atomicMax((unsigned long long*)addr, __double_as_longlong(value))
    );
}


/* ------------------------------------------------------------------------
 * Kernel to find maxval(abs(A))
 * ------------------------------------------------------------------------
 */
template<typename Number, typename Real = typename get_real<Number>::type>
__global__ void maxabs_kernel(const Number* __restrict__ data,
                              size_t n,
                              Real* __restrict__ norm)
{
    extern __shared__ unsigned char smem[];
    Real* sdata = reinterpret_cast<Real*>(smem);

    unsigned int tid = threadIdx.x;
    size_t i = blockDim.x * blockIdx.x + tid;

    Real val = 0;
    if (i < n) val = absval(data[i]);

    sdata[tid] = val;

    __syncthreads();

    // reduction to find block maximum
    for (unsigned int s = blockDim.x / 2; s > 0; s >>= 1)
    {
        if (tid < s)
            sdata[tid] = fmax(sdata[tid], sdata[tid + s]);
        __syncthreads();
    }

    // final block max reduced atomically to global result
    if (tid == 0) {
        if constexpr (std::is_same<Real, float>::value)
            atomicMaxAbsFloat(norm, sdata[0]);
        else
            atomicMaxAbsDouble(norm, sdata[0]);
    }
}

/*
 * DATA MOVEMENT ROUTINES
 */

extern "C" int cu_createMat(void** d_A, size_t bytecount) {
    ENFORCE(d_A);
    cudaError_t err = cudaMalloc(d_A, bytecount);
    CUDA_CHECK(err);
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
    CUDA_CHECK(err);
    return err;
}

extern "C" int cu_copyMatH2D(void* h_A, void* d_A, size_t bytecount) {
    ENFORCE(h_A);
    ENFORCE(d_A);
    cudaError_t err = cudaMemcpy(d_A, h_A, bytecount, cudaMemcpyHostToDevice);
    CUDA_CHECK(err);
    return err;
}

extern "C" int cu_copyMatH2D_async(void* h_A, void* d_A, size_t bytecount) {
    ENFORCE(h_A);
    ENFORCE(d_A);
    cudaError_t err =
        cudaMemcpyAsync(d_A, h_A, bytecount, cudaMemcpyHostToDevice);
    return err;
}

extern "C" int cu_copyMatD2H(void* h_A, void* d_A, size_t bytecount) {
    ENFORCE(h_A);
    ENFORCE(d_A);
    cudaError_t err = cudaMemcpy(h_A, d_A, bytecount, cudaMemcpyDeviceToHost);
    CUDA_CHECK(err);
    return err;
}

extern "C" int cu_copyMatD2D_async(void* d_A, void* d_A_in, size_t bytecount) {
    ENFORCE(d_A);
    ENFORCE(d_A_in);
    cudaError_t err =
        cudaMemcpyAsync(d_A, d_A_in, bytecount, cudaMemcpyDeviceToDevice);
    return err;
}

extern "C" int cu_copyMatD2D(void* d_A, void* d_A_in, size_t bytecount) {
    ENFORCE(d_A);
    ENFORCE(d_A_in);
    cudaError_t err = cudaMemcpy(d_A, d_A_in, bytecount, cudaMemcpyDeviceToDevice);
    CUDA_CHECK(err);
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
    CUDA_CHECK(err);
    return err;
}

extern "C" int cu_cudaGetDeviceProperties(int device) {
    cudaDeviceProp prop;
    cudaError_t err = cudaGetDeviceProperties(&prop, device);
    CUDA_CHECK(err);

    printf(" Found GPU: Device Name: %s\n", prop.name);
    printf(" TotalMemory: %lu\n", (unsigned long)prop.totalGlobalMem);
    printf(" Shared per block: %lu\n", (unsigned long)prop.sharedMemPerBlock);

    return err;
}

extern "C" int cu_cudaGetDevice(int* device) {
    cudaError_t err = cudaGetDevice(device);
    CUDA_CHECK(err);
    return err;
}

extern "C" int cu_cudaSetDevice(int count) {
    cudaError_t err = cudaSetDevice(count);
    CUDA_CHECK(err);
    return err;
}

extern "C" int cu_cublasInit(void** hpblas) {
    ENFORCE(hpblas);

    cublasHandle_t hb;
    cublasStatus_t err = cublasCreate(&hb);
    ENFORCE(err == CUBLAS_STATUS_SUCCESS);
    if(err != CUBLAS_STATUS_SUCCESS) {
        printf("cublas create error: %d\n", err);
    }
    // printf("hcublas Addr: %p \n",*hcublas);
    *hpblas = hb;
    return err;
}

extern "C" int cu_cublasFinalize(void *hpblas) {
    ENFORCE(hpblas);
    cublasHandle_t hcublas = static_cast<cublasHandle_t>(hpblas);
    cublasStatus_t err = cublasDestroy(hcublas);
    ENFORCE(err == CUBLAS_STATUS_SUCCESS);
    return err;
}

extern "C" int cu_cusolverInit(void** hpsolver) {
    ENFORCE(hpsolver);
    cusolverDnHandle_t hs;
    cusolverStatus_t err = cusolverDnCreate(&hs);
    ENFORCE(err == CUSOLVER_STATUS_SUCCESS);
    if(err != 0) {
        printf("cusolver create error: %d\n", err);
    }
    // printf("hcusolver Addr: %p \n",*hcusolver);
    *hpsolver = hs;
    return err;
}

extern "C" int cu_cusolverFinalize(void *hpsolver) {
    ENFORCE(hpsolver);
    cusolverDnHandle_t hcusolver = static_cast<cusolverDnHandle_t>(hpsolver);
    cusolverStatus_t err = cusolverDnDestroy(hcusolver);
    ENFORCE(err == CUSOLVER_STATUS_SUCCESS);
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
 * * `dagger == 0`: compute C = α A · B + β C
 * * `dagger == 1`: compute C = α A^* · B + β C
 * * `dagger == 2`: compute C = α A · B^* + β C
 *
 * @param[in] m The number of rows of C and op(A).
 * @param[in] n The number of columns of C and op(B).
 * @param[in] k The number of columns of op(A) and the number of rows of op(B).
 * @param[in] dagger A shorthand for various combinations of op(A), op(B).
 */
template<typename Number>
int cu_multMat(
    void *hpblas, int32_t m, int32_t n, int32_t k, const Number* alpha,
    const Number* d_A, const Number* d_B, const Number* beta, Number* d_C,
    int dagger
) {
    ENFORCE(alpha);
    ENFORCE(d_A);
    ENFORCE(d_B);
    ENFORCE(beta);
    ENFORCE(dagger == 0 || dagger == 1 || dagger == 2);

    cublasHandle_t hcublas = static_cast<cublasHandle_t>(hpblas);
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
    void* hpblas, int32_t m, int32_t n, int32_t k,
    const cuComplex* alpha, const cuComplex* d_A, const cuComplex* d_B,
    const cuComplex* beta, cuComplex* d_C, int dagger
) {
    return cu_multMat(hpblas, m, n, k, alpha, d_A, d_B, beta, d_C, dagger);
}

extern "C" int cu_ZmultMat(
    void* hpblas, int32_t m, int32_t n, int32_t k,
    const cuDoubleComplex* alpha, const cuDoubleComplex* d_A,
    const cuDoubleComplex* d_B, const cuDoubleComplex* beta,
    cuDoubleComplex* d_C, int dagger
) {
    return cu_multMat(hpblas, m, n, k, alpha, d_A, d_B, beta, d_C, dagger);
}



template<typename Number>
int cu_kernelsum(
    Number* d_C, Number* alpha, Number* d_A, Number* beta, Number* d_B,
    size_t num_elements
) {
    auto num_blocks = (num_elements + BLOCK_SIZE - 1) / BLOCK_SIZE;

    addKernel<<<num_blocks, BLOCK_SIZE>>>(d_C, *alpha, d_A, *beta, d_B, num_elements);
    CUDA_CHECK(cudaGetLastError());

    return 0;
}


extern "C" int cu_Ckernelsum(
    cuComplex* d_C, cuComplex* alpha, cuComplex* d_A, cuComplex* beta,
    cuComplex* d_B, size_t num_elements
) {
    return cu_kernelsum(d_C, alpha, d_A, beta, d_B, num_elements);
}

extern "C" int cu_Zkernelsum(
    cuDoubleComplex* d_C, cuDoubleComplex* alpha, cuDoubleComplex* d_A,
    cuDoubleComplex* beta, cuDoubleComplex* d_B, size_t num_elements
) {
    return cu_kernelsum(d_C, alpha, d_A, beta, d_B, num_elements);
}


/**
 * Computes the sum of α op(A) + β op(B), where op(.) indicates if the matrix
 * or its complex-conjugate is used.
 *
 * The possible options are:
 * * `dagger == 0`: compute C = α A + β B
 * * `dagger == 1`: compute C = α A^* + β B
 * * `dagger == 2`: compute C = α A + B^*
 *
 * @param[in] m The number of rows of op(A) and op(B).
 * @param[in] n The number of columns of op(A) and op(B).
 * @param[in] dagger A shorthand for various combinations of op(A), op(B).
 */
template<typename Number>
int cu_matsum(
    void *hpblas, int32_t m, int32_t n, const Number* alpha,
    const Number* d_A, const Number* beta, const Number* d_B, Number* d_C,
    int dagger
) {
    ENFORCE(d_A);
    ENFORCE(d_B);
    ENFORCE(d_C);
    ENFORCE(d_A != d_C || dagger == 0 || dagger == 2);
    ENFORCE(d_B != d_C || dagger == 1);
    cublasHandle_t hcublas = static_cast<cublasHandle_t>(hpblas);
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
    void *hpblas, int32_t m, int32_t n, const cuComplex* alpha,
    const cuComplex* d_A, const cuComplex* beta, const cuComplex* d_B,
    cuComplex* d_C, int dagger
) {
    return cu_matsum(hpblas, m, n, alpha, d_A, beta, d_B, d_C, dagger);
}

extern "C" int cu_Zmatsum(
    void *hpblas, int32_t m, int32_t n, const cuDoubleComplex* alpha,
    const cuDoubleComplex* d_A, const cuDoubleComplex* beta,
    const cuDoubleComplex* d_B, cuDoubleComplex* d_C, int dagger
) {
    return cu_matsum(hpblas, m, n, alpha, d_A, beta, d_B, d_C, dagger);
}

extern "C" int cu_Cinitmat(cuComplex* d_A, int32_t nrow) {
    ENFORCE(d_A);
    size_t n64 = static_cast<size_t>(nrow);
    size_t num_elements = n64 * n64;
    size_t num_blocks = (num_elements + BLOCK_SIZE - 1) / BLOCK_SIZE;

    initKernel<<<num_blocks, BLOCK_SIZE>>>(d_A, nrow);
    CUDA_CHECK(cudaGetLastError());

    return 0;
}

extern "C" int cu_Zinitmat(cuDoubleComplex* d_A, int32_t nrow) {
    ENFORCE(d_A);
    size_t n64 = static_cast<size_t>(nrow);
    size_t num_elements = n64 * n64;
    size_t num_blocks = (num_elements + BLOCK_SIZE - 1) / BLOCK_SIZE;

    initKernel<<<num_blocks, BLOCK_SIZE>>>(d_A, nrow);
    CUDA_CHECK(cudaGetLastError());

    return 0;
}



template<typename Number, typename Real = typename get_real<Number>::type>
Real trace(
    Number* d_A, int32_t nrow, void* h_mask,
    int mask_present
) {

    size_t n64 = static_cast<size_t>(nrow);
    size_t num_elements = n64 * n64;
    size_t num_blocks = (num_elements + BLOCK_SIZE - 1) / BLOCK_SIZE;

    cudaError_t cudaStatus;

    bool* d_mask;
    cudaStatus = cudaMalloc((void**)&d_mask, n64 * sizeof(bool));
    CUDA_CHECK(cudaStatus);
    if (h_mask) {
       cu_copyMatH2D(h_mask, d_mask, n64*sizeof(bool));
    }

    Real h_result;
    Real *d_result;
    Real *d_real_diag;
    void *d_temp = nullptr;
    size_t temp_bytes = 0;
    cudaStatus = cudaMalloc(&d_result, sizeof(Real));
    cudaStatus = cudaMemset(static_cast<void*>(d_result), 0, sizeof(Real));
    CUDA_CHECK(cudaStatus);
    cudaStatus = cudaMalloc(&d_real_diag, n64*sizeof(Real));
    CUDA_CHECK(cudaStatus);

    getdiagKernel<<<num_blocks, BLOCK_SIZE>>>(d_A, nrow, d_real_diag, d_mask, mask_present);

    cub::DeviceReduce::Sum(nullptr, temp_bytes, d_real_diag, d_result, nrow);

    cudaStatus = cudaMalloc(&d_temp, temp_bytes);
    CUDA_CHECK(cudaStatus);

    cub::DeviceReduce::Sum(d_temp, temp_bytes, d_real_diag, d_result, nrow);

    cudaStatus = cudaMemcpy(&h_result, d_result, sizeof(Real), cudaMemcpyDeviceToHost);

    cudaStatus = cudaFree(d_mask);
    cudaStatus = cudaFree(d_result);
    cudaStatus = cudaFree(d_temp);
    cudaStatus = cudaFree(d_real_diag);

    return h_result;
}

extern "C" float cu_Ctrace(
    void *hpblas, cuComplex* d_A, int32_t nrow, void* h_mask,
    int mask_present
) {
    return trace(d_A, nrow, h_mask, mask_present);
}


extern "C" double cu_Ztrace(
    void *hpblas, cuDoubleComplex* d_A, int32_t nrow, void* h_mask,
    int mask_present
) {
    return trace(d_A, nrow, h_mask, mask_present);
}


extern "C" int cu_Cmatcopy(
    void *hpblas, const cuComplex* d_A, cuComplex* d_B,
    size_t num_elements
) {
    ENFORCE(d_A);
    ENFORCE(d_B);
    ENFORCE(num_elements <= INT_MAX);
    cublasHandle_t hcublas = static_cast<cublasHandle_t>(hpblas);

    auto err = cublasCcopy(hcublas, num_elements, d_A, 1, d_B, 1);
    ENFORCE(err == CUBLAS_STATUS_SUCCESS);
    return err;
}

extern "C" int cu_Zmatcopy(
    void *hpblas, cuDoubleComplex* d_A, cuDoubleComplex* d_B,
    size_t num_elements
) {
    ENFORCE(d_A);
    ENFORCE(d_B);
    ENFORCE(num_elements <= INT_MAX);
    cublasHandle_t hcublas = static_cast<cublasHandle_t>(hpblas);

    auto err = cublasZcopy(hcublas, num_elements, d_A, 1, d_B, 1);
    ENFORCE(err == CUBLAS_STATUS_SUCCESS);
    return err;
}

extern "C" int
cu_Casum(void *hpblas, void* d_A, float* summ, int32_t n) {
    cuComplex* pdA = (cuComplex*)d_A;

    cublasHandle_t hcublas = static_cast<cublasHandle_t>(hpblas);
    auto err = cublasScasum(hcublas, n, pdA, 1, summ);
    ENFORCE(err == CUBLAS_STATUS_SUCCESS);
    return err;
}

extern "C" int
cu_Zasum(void *hpblas, void* d_A, double* summ, int32_t n) {
    cuDoubleComplex* pdA = (cuDoubleComplex*)d_A;

    cublasHandle_t hcublas = static_cast<cublasHandle_t>(hpblas);
    auto err = cublasDzasum(hcublas, n, pdA, 1, summ);
    ENFORCE(err == CUBLAS_STATUS_SUCCESS);
    return err;
}

/* **************************************************************************
 * Matrix inversion  A^-1
 * **************************************************************************
 */

template<typename Number>
int inverse(void* hpsolver, Number* d_A, Number* d_Ainv, int32_t n
) {
    ENFORCE(hpsolver);
    ENFORCE(d_A);
    ENFORCE(d_Ainv);

    size_t n64 = static_cast<size_t>(n);
    size_t num_elements = n64 * n64;
    size_t bytecount = num_elements * sizeof(Number);

    cusolverDnHandle_t hcusolver = static_cast<cusolverDnHandle_t>(hpsolver);

    cusolverStatus_t cusolverStatus;
    cudaError_t cudaStatus;

    // compute buffer msize and prep . memory
    int lwork;
    cusolverStatus =
        libnegf::cusolverDngetrf_bufferSize(hcusolver, n, n, d_A, n, &lwork);
    ENFORCE(cusolverStatus == CUSOLVER_STATUS_SUCCESS);

    // prepare memory on the device
    Number* d_LU;
    cudaStatus = cudaMalloc((void**)&d_LU, bytecount);
    CUDA_CHECK(cudaStatus);

    int* d_pivot;
    cudaStatus = cudaMalloc((void**)&d_pivot, n64 * sizeof(int));
    CUDA_CHECK(cudaStatus);
    int* d_info;
    cudaStatus = cudaMalloc((void**)&d_info, sizeof(int));
    CUDA_CHECK(cudaStatus);
    // copy d_LU <- pdA
    cu_copyMatD2D(d_LU, d_A, bytecount);

    Number* d_work;
    size_t lwork_size = static_cast<size_t>(lwork) * sizeof(Number);
    cudaStatus = cudaMalloc((void**)&d_work, lwork_size);
    CUDA_CHECK(cudaStatus);

    // LU factorization of d_A , with partial pivoting and row
    // interchanges ; row i is interchanged with row d_pivot ( i );
    cusolverStatus = libnegf::cusolverDngetrf(
        hcusolver, n, n, d_LU, n, d_work, d_pivot, d_info
    );
    ENFORCE(cusolverStatus == CUSOLVER_STATUS_SUCCESS);
    LIBNEGF_CHECK_INFO(d_info);

    // use the LU factorization to solve the system d_LU * x = d_Ainv ;
    // the solution overwrites d_Ainv
    cusolverStatus = libnegf::cusolverDngetrs(
        hcusolver, CUBLAS_OP_N, n, n, d_LU, n, d_pivot, d_Ainv, n, d_info
    );
    ENFORCE(cusolverStatus == CUSOLVER_STATUS_SUCCESS);
    LIBNEGF_CHECK_INFO(d_info);

    cudaStatus = cudaFree(d_pivot);
    CUDA_CHECK(cudaStatus);
    cudaStatus = cudaFree(d_info);
    CUDA_CHECK(cudaStatus);
    cudaStatus = cudaFree(d_work);
    CUDA_CHECK(cudaStatus);
    cudaStatus = cudaFree(d_LU);
    CUDA_CHECK(cudaStatus);

    return cudaStatus;
}

extern "C" int cu_Cinverse(
    void* hpsolver, cuComplex* d_A,
    cuComplex* d_Ainv, int32_t n
) {
    return inverse(hpsolver, d_A, d_Ainv, n);
}

extern "C" int cu_Zinverse(
    void * hpsolver, cuDoubleComplex* d_A,
    cuDoubleComplex* d_Ainv, int32_t n
) {
    return inverse(hpsolver, d_A, d_Ainv, n);
}

/* ***************************************************************
 *  Sanchez-Sanchez algorithm
 * ****************************************************************
 */
template<typename Number, typename Real = typename get_real<Number>::type>
int decimation(
    void* hpblas, void* hpsolver, Number* h_Go_out,
    Number* h_Ao_in, Number* h_Bo_in, Number* h_Co_in, int32_t n,
    int* ncyc, Real accuracy, bool log
) {
    ENFORCE(h_Go_out);
    ENFORCE(h_Ao_in);
    ENFORCE(h_Bo_in);
    ENFORCE(h_Co_in);
    ENFORCE(ncyc);
    ENFORCE(accuracy > 0.0);

    cublasHandle_t hcublas = static_cast<cublasHandle_t>(hpblas);
    cusolverDnHandle_t hcusolver = static_cast<cusolverDnHandle_t>(hpsolver);

    cudaError_t cudaStatus;
    cublasStatus_t err;

    size_t n64 = static_cast<size_t>(n);
    size_t num_elements = n64 * n64;
    size_t bytecount = num_elements * sizeof(Number);
    size_t num_blocks = (num_elements + BLOCK_SIZE - 1) / BLOCK_SIZE;

    auto one = Number{Real{1}};
    auto mone = Number{-Real{1}};
    auto zero = Number{Real{0}};

    Number* d_Ao;
    cudaStatus = cudaMalloc((void**)&d_Ao, bytecount);
    CUDA_CHECK(cudaStatus);

    Number* d_Bo;
    cudaStatus = cudaMalloc((void**)&d_Bo, bytecount);
    CUDA_CHECK(cudaStatus);

    Number* d_Co;
    cudaStatus = cudaMalloc((void**)&d_Co, bytecount);
    CUDA_CHECK(cudaStatus);

    cu_copyMatH2D(h_Ao_in, d_Ao, bytecount);
    cu_copyMatH2D(h_Bo_in, d_Bo, bytecount);
    cu_copyMatH2D(h_Co_in, d_Co, bytecount);

    Number* d_Ao_s;
    cudaStatus = cudaMalloc((void**)&d_Ao_s, bytecount);
    CUDA_CHECK(cudaStatus);
    Number* d_C1;
    cudaStatus = cudaMalloc((void**)&d_C1, bytecount);
    CUDA_CHECK(cudaStatus);
    Number* d_Go;
    cudaStatus = cudaMalloc((void**)&d_Go, bytecount);
    CUDA_CHECK(cudaStatus);
    int* d_pivot;
    cudaStatus = cudaMalloc((void**)&d_pivot, n64 * sizeof(int));
    CUDA_CHECK(cudaStatus);
    Number* d_T;
    cudaStatus = cudaMalloc((void**)&d_T, bytecount);
    CUDA_CHECK(cudaStatus);
    Number* d_Self;
    cudaStatus = cudaMalloc((void**)&d_Self, bytecount);
    CUDA_CHECK(cudaStatus);
    int* d_info;
    cudaStatus = cudaMalloc((void**)&d_info, sizeof(int));
    CUDA_CHECK(cudaStatus);
    Real* d_norm;
    cudaStatus = cudaMalloc((void**)&d_norm, sizeof(Real));
    CUDA_CHECK(cudaStatus);

    int lwork;
    cusolverStatus_t cusolverStatus =
            libnegf::cusolverDngetrf_bufferSize(hcusolver, n, n, d_Self, n, &lwork);
    ENFORCE(cusolverStatus == CUSOLVER_STATUS_SUCCESS);
    Number* d_work;
    cudaStatus = cudaMalloc((void**)&d_work, lwork * sizeof(Number));
    CUDA_CHECK(cudaStatus);

    cu_copyMatD2D(d_Ao_s, d_Ao, bytecount);

    constexpr auto MAX_NUM_ITERATIONS = 300;
    bool okCo = false;
    for(int i1 = 1; i1 <= MAX_NUM_ITERATIONS; i1++) {
        *ncyc = i1;

        // Compute Go = Ao^-1
        initKernel<<<num_blocks, BLOCK_SIZE>>>(d_Go, n);

        cu_copyMatD2D(d_Self, d_Ao, bytecount);

        cusolverStatus = libnegf::cusolverDngetrf(
            hcusolver, n, n, d_Self, n, d_work, d_pivot, d_info
        );
        ENFORCE(cusolverStatus == CUSOLVER_STATUS_SUCCESS);
        LIBNEGF_CHECK_INFO(d_info);
        cusolverStatus = libnegf::cusolverDngetrs(
            hcusolver, CUBLAS_OP_N, n, n, d_Self, n, d_pivot, d_Go, n, d_info
        );
        ENFORCE(cusolverStatus == CUSOLVER_STATUS_SUCCESS);
        LIBNEGF_CHECK_INFO(d_info);

        // compute C1 = Co Go Co
        err = libnegf::cublasGemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, &one, d_Go, n, d_Co, n,
            &zero, d_T, n
        );
        ENFORCE(err == CUBLAS_STATUS_SUCCESS);

        err = libnegf::cublasGemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, &one, d_Co, n, d_T, n,
            &zero, d_C1, n
        );
        ENFORCE(err == CUBLAS_STATUS_SUCCESS);

        // convergence check on |C1| < tol
        cudaStatus = cudaMemset(d_norm, 0, sizeof(Real));
        CUDA_CHECK(cudaStatus);

        maxabs_kernel<<<num_blocks, BLOCK_SIZE, BLOCK_SIZE*sizeof(Real)>>>
                                                       (d_C1, n * n, d_norm);
        CUDA_CHECK(cudaGetLastError());
        Real h_norm;
        cudaStatus = cudaMemcpy(&h_norm, d_norm, sizeof(Real), cudaMemcpyDeviceToHost);
        if (log){printf("loop it= %d , maxval|C1|= %g \n ", i1, h_norm);}

        if(h_norm <= accuracy) {
            if(okCo) {
                break;
            } else {
                okCo = true;
            }
        } else {
            okCo = false;
        }

        // compute Self = Bo Go Co
        err = libnegf::cublasGemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, &one, d_Bo, n, d_T, n,
            &zero, d_Self, n
        );
        ENFORCE(err == CUBLAS_STATUS_SUCCESS);

        // compute Ao_s = Ao_s - Self
        err = libnegf::cublasAxpy(hcublas, n * n, &mone, d_Self, 1, d_Ao_s, 1);
        ENFORCE(err == CUBLAS_STATUS_SUCCESS);

	  // compute Ao = Ao - Self
        err = libnegf::cublasAxpy(hcublas, n * n, &mone, d_Self, 1, d_Ao, 1);
        ENFORCE(err == CUBLAS_STATUS_SUCCESS);

        // compute Ao = Ao - Co Go Bo
        err = libnegf::cublasGemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, &one, d_Go, n, d_Bo, n,
            &zero, d_T, n
        );
        ENFORCE(err == CUBLAS_STATUS_SUCCESS);

        err = libnegf::cublasGemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, &mone, d_Co, n, d_T, n,
            &one, d_Ao, n
        );
        ENFORCE(err == CUBLAS_STATUS_SUCCESS);

        // Co = C1
        cu_copyMatD2D(d_Co, d_C1, bytecount);

        // compute C1 = Bo Go Bo
        err = libnegf::cublasGemm(
            hcublas, CUBLAS_OP_N, CUBLAS_OP_N, n, n, n, &one, d_Bo, n, d_T, n,
            &zero, d_C1, n
        );
        ENFORCE(err == CUBLAS_STATUS_SUCCESS);

        // compute Bo = C1
        cu_copyMatD2D(d_Bo, d_C1, bytecount);
    }

    ENFORCE(*ncyc < MAX_NUM_ITERATIONS);

    // compute Go = Ao_s^-1
    initKernel<<<num_blocks, BLOCK_SIZE>>>(d_Go, n);
    cu_copyMatD2D(d_Self, d_Ao_s, bytecount);
    cusolverStatus = libnegf::cusolverDngetrf(
        hcusolver, n, n, d_Self, n, d_work, d_pivot, d_info
    );
    ENFORCE(cusolverStatus == CUSOLVER_STATUS_SUCCESS);
    LIBNEGF_CHECK_INFO(d_info);
    cusolverStatus = libnegf::cusolverDngetrs(
        hcusolver, CUBLAS_OP_N, n, n, d_Self, n, d_pivot, d_Go, n, d_info
    );
    ENFORCE(cusolverStatus == CUSOLVER_STATUS_SUCCESS);
    LIBNEGF_CHECK_INFO(d_info);

    // copy to output to host
    cu_copyMatD2H(h_Go_out, d_Go, bytecount);

    cudaStatus = cudaFree(d_pivot);
    CUDA_CHECK(cudaStatus);
    cudaStatus = cudaFree(d_info);
    CUDA_CHECK(cudaStatus);
    cudaStatus = cudaFree(d_Ao);
    CUDA_CHECK(cudaStatus);
    cudaStatus = cudaFree(d_Bo);
    CUDA_CHECK(cudaStatus);
    cudaStatus = cudaFree(d_Co);
    CUDA_CHECK(cudaStatus);
    cudaStatus = cudaFree(d_Go);
    CUDA_CHECK(cudaStatus);
    cudaStatus = cudaFree(d_Ao_s);
    CUDA_CHECK(cudaStatus);
    cudaStatus = cudaFree(d_C1);
    CUDA_CHECK(cudaStatus);
    cudaStatus = cudaFree(d_T);
    CUDA_CHECK(cudaStatus);
    cudaStatus = cudaFree(d_Self);
    CUDA_CHECK(cudaStatus);
    cudaStatus = cudaFree(d_work);
    CUDA_CHECK(cudaStatus);
    cudaStatus = cudaFree(d_norm);
    CUDA_CHECK(cudaStatus);

    return cudaStatus;
}


extern "C" int cu_Cdecimation(
    void* hcublas, void* hcusolver, cuComplex* h_Go_out,
    cuComplex* h_Ao_in, cuComplex* h_Bo_in, cuComplex* h_Co_in, int32_t n,
    int* ncyc, float accuracy, bool log = false
) {
    return decimation(
        hcublas, hcusolver, h_Go_out, h_Ao_in, h_Bo_in, h_Co_in, n, ncyc,
        accuracy, log
    );
}

extern "C" int cu_Zdecimation(
    void* hcublas, void* hcusolver,
    cuDoubleComplex* h_Go_out, cuDoubleComplex* h_Ao_in,
    cuDoubleComplex* h_Bo_in, cuDoubleComplex* h_Co_in, int32_t n,
    int* ncyc, double accuracy, bool log = false
) {
    return decimation(
        hcublas, hcusolver, h_Go_out, h_Ao_in, h_Bo_in, h_Co_in, n, ncyc,
        accuracy, log
    );
}


extern "C" int cu_meminfo(size_t* freemem, size_t* totalmem) {
    //auto cudaStatus = cudaDeviceSynchronize();
    //printf("STATUS: %d\n",cudaStatus);
    //CUDA_CHECK(cudaStatus);
    auto cudaStatus = cudaMemGetInfo(freemem, totalmem);
    CUDA_CHECK(cudaStatus);
    return cudaStatus;
}
