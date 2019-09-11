include(CheckCXXCompilerFlag)

check_cxx_compiler_flag("-march=native" MARLIN_COMPILER_SUPPORTS_MARCH_NATIVE)
check_cxx_compiler_flag("-mtune=native" MARLIN_COMPILER_SUPPORTS_MTUNE_NATIVE)
