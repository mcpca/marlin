set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
set(THREADS_PREFER_PTHREAD_FLAG TRUE)
find_package(Threads REQUIRED)

add_library(thread_pool INTERFACE)

target_link_libraries(thread_pool INTERFACE Threads::Threads)
target_include_directories(thread_pool INTERFACE include)
