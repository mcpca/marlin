add_subdirectory(external)

find_package(h5cpp REQUIRED)

add_library(fsm
            src/solver.cpp
            src/solver_internals.cpp
            src/grid.cpp
            src/data.cpp
            src/io.cpp)

target_compile_features(fsm PUBLIC cxx_std_14)

target_compile_options(fsm
                       PRIVATE -Wall
                               -Wextra
                               -Werror
                               -pedantic)

target_include_directories(fsm PUBLIC include PRIVATE src)
target_link_libraries(fsm
                      PRIVATE h5cpp
                      PUBLIC thread_pool
                             concurrent_queue)
