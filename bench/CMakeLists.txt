add_library(bench_lib OBJECT EXCLUDE_FROM_ALL bench.cpp)

add_executable(bench_eikonal3d EXCLUDE_FROM_ALL eikonal3d.cpp)
target_compile_features(bench_eikonal3d PRIVATE cxx_std_14)
target_compile_options(bench_eikonal3d PRIVATE -Wall -Wextra -pedantic)
target_link_libraries(bench_eikonal3d PRIVATE marlin bench_lib)

add_custom_target(bench)
add_dependencies(bench bench_eikonal3d)
