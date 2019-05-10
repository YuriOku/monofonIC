#pragma once

#include <logger.hh>

#include <complex>

#if defined(USE_MPI)
#include <mpi.h>
  #include <fftw3-mpi.h>
#else
  #include <fftw3.h>
#endif

#ifdef USE_SINGLEPRECISION
using real_t = float;
using complex_t = fftwf_complex;
#define FFTW_PREFIX fftwf
#else
using real_t = double;
using complex_t = fftw_complex;
#define FFTW_PREFIX fftw
#endif

using ccomplex_t = std::complex<real_t>;

#define FFTW_GEN_NAME_PRIM(a, b) a##_##b
#define FFTW_GEN_NAME(a, b) FFTW_GEN_NAME_PRIM(a, b)
#define FFTW_API(x) FFTW_GEN_NAME(FFTW_PREFIX, x)

using fftw_plan_t = FFTW_GEN_NAME(FFTW_PREFIX, plan);

#if defined(FFTW_MODE_PATIENT)
#define FFTW_RUNMODE FFTW_PATIENT
#elif defined(FFTW_MODE_MEASURE)
#define FFTW_RUNMODE FFTW_MEASURE
#else
#define FFTW_RUNMODE FFTW_ESTIMATE
#endif

#if defined(USE_MPI)
inline double get_wtime()
{
    return MPI_Wtime();
}
#else
  #if defined(_OPENMP)
    #include <omp.h>
    inline double get_wtime()
    {
      return omp_get_wtime();
    }
  #else
    #include <ctime>
    inline double get_wtime()
    {
      return std::clock() / double(CLOCKS_PER_SEC);
    }
  #endif
#endif

namespace CONFIG
{
extern int MPI_thread_support;
extern int MPI_task_rank;
extern int MPI_task_size;
extern bool MPI_ok;
extern bool MPI_threads_ok;
extern bool FFTW_threads_ok;
}; // namespace CONFIG