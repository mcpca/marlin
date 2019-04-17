set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)

add_library(concurrent_queue INTERFACE)

target_link_libraries(concurrent_queue INTERFACE Threads::Threads)
target_include_directories(concurrent_queue INTERFACE include)
