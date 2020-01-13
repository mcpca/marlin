find_package(h5cpp)

if(h5cpp_FOUND)
  add_library(hdf5_helper EXCLUDE_FROM_ALL OBJECT hdf5.cpp)
  target_compile_features(hdf5_helper PRIVATE cxx_std_14)
  target_compile_options(hdf5_helper PRIVATE -Wall -Wextra -pedantic)
  target_link_libraries(hdf5_helper PUBLIC marlin h5cpp)

  add_executable(eikonal2d EXCLUDE_FROM_ALL eikonal2d.cpp)
  target_compile_features(eikonal2d PRIVATE cxx_std_14)
  target_compile_options(eikonal2d PRIVATE -Wall -Wextra -pedantic)
  target_link_libraries(eikonal2d PRIVATE marlin hdf5_helper)

  add_executable(eikonal3d EXCLUDE_FROM_ALL eikonal3d.cpp)
  target_compile_features(eikonal3d PRIVATE cxx_std_14)
  target_compile_options(eikonal3d PRIVATE -Wall -Wextra -pedantic)
  target_link_libraries(eikonal3d PRIVATE marlin hdf5_helper)

  add_executable(constant_field EXCLUDE_FROM_ALL constant_field.cpp)
  target_compile_features(constant_field PRIVATE cxx_std_14)
  target_compile_options(constant_field PRIVATE -Wall -Wextra -pedantic)
  target_link_libraries(constant_field PRIVATE marlin hdf5_helper)

  add_executable(geodesic EXCLUDE_FROM_ALL geodesic.cpp)
  target_compile_features(geodesic PRIVATE cxx_std_14)
  target_compile_options(geodesic PRIVATE -Wall -Wextra -pedantic)
  target_link_libraries(geodesic PRIVATE marlin hdf5_helper)

  add_custom_target(examples)
  add_dependencies(examples eikonal2d eikonal3d constant_field geodesic)
endif()
