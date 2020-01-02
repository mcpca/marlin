include(CheckCXXCompilerFlag)

check_cxx_compiler_flag("-march=native" MARLIN_COMPILER_SUPPORTS_MARCH_NATIVE)
check_cxx_compiler_flag("-mtune=native" MARLIN_COMPILER_SUPPORTS_MTUNE_NATIVE)

if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  if(MARLIN_INTEL_TARGET STREQUAL "SKX")
    check_cxx_compiler_flag("-xMIC-AVX512"
                            MARLIN_COMPILER_SUPPORTS_INTEL_MIC_AVX512)

    if(MARLIN_COMPILER_SUPPORTS_INTEL_MIC_AVX512)
      set(MARLIN_USE_INTEL_MIC_AVX512 ON)
    endif()

  elseif(MARLIN_INTEL_TARGET STREQUAL "KNL")
    check_cxx_compiler_flag("-xCORE-AVX512"
                            MARLIN_COMPILER_SUPPORTS_INTEL_CORE_AVX512)

    if(MARLIN_COMPILER_SUPPORTS_INTEL_CORE_AVX512)
      set(MARLIN_USE_INTEL_CORE_AVX512 ON)
    endif()
  endif()

  if(MARLIN_USE_INTEL_MKL)
    check_cxx_compiler_flag("-mkl" MARLIN_COMPILER_SUPPORTS_INTEL_MKL)

    if(NOT MARLIN_COMPILER_SUPPORTS_INTEL_MKL)
      set(MARLIN_USE_INTEL_MKL OFF)
    endif()
  endif()
else()
  set(MARLIN_USE_INTEL_MKL OFF)
endif()
