add_executable(eikonal2d EXCLUDE_FROM_ALL eikonal2d.cpp)
target_compile_features(eikonal2d PUBLIC cxx_std_14)
target_compile_options(eikonal2d PRIVATE -Wall -Wextra -pedantic)
target_link_libraries(eikonal2d PUBLIC fsm)

add_executable(eikonal3d EXCLUDE_FROM_ALL eikonal3d.cpp)
target_compile_features(eikonal3d PUBLIC cxx_std_14)
target_compile_options(eikonal3d PRIVATE -Wall -Wextra -pedantic)
target_link_libraries(eikonal3d PUBLIC fsm)

add_executable(constant_field EXCLUDE_FROM_ALL constant_field.cpp)
target_compile_features(constant_field PUBLIC cxx_std_14)
target_compile_options(constant_field PRIVATE -Wall -Wextra -pedantic)
target_link_libraries(constant_field PUBLIC fsm)

add_executable(geodesic EXCLUDE_FROM_ALL geodesic.cpp)
target_compile_features(geodesic PUBLIC cxx_std_14)
target_compile_options(geodesic PRIVATE -Wall -Wextra -pedantic)
target_link_libraries(geodesic PUBLIC fsm)

add_custom_target(examples)
add_dependencies(examples
                 eikonal2d
                 eikonal3d
                 constant_field
                 geodesic)
