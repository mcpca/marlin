find_package(h5cpp REQUIRED)
find_package(OpenMP REQUIRED)

include(cmake/Compiler.cmake)

add_library(marlin src/solver.cpp src/solver_internals.cpp src/io.cpp)

target_compile_features(marlin PUBLIC cxx_std_14)

target_compile_options(marlin PRIVATE -Wall -Wextra -Werror -pedantic)

# Intel compiler options
if(MARLIN_USE_INTEL_MIC_AVX512)
  target_compile_options(marlin PRIVATE -xMIC-AVX512)
elseif(MARLIN_USE_INTEL_CORE_AVX512)
  target_compile_options(marlin PRIVATE -xCORE-AVX512)
else()
  if(MARLIN_COMPILER_SUPPORTS_MARCH_NATIVE)
    target_compile_options(marlin PRIVATE -march=native)
  endif()

  if(MARLIN_COMPILER_SUPPORTS_MTUNE_NATIVE)
    target_compile_options(marlin PRIVATE -mtune=native)
  endif()
endif()

if(MARLIN_USE_INTEL_MKL)
  target_compile_options(marlin PRIVATE -mkl)
endif()

target_include_directories(
  marlin
  PUBLIC include
  PRIVATE src)

target_link_libraries(
  marlin
  PRIVATE h5cpp
  PUBLIC OpenMP::OpenMP_CXX)
