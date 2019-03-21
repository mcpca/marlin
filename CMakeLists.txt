cmake_minimum_required(VERSION 3.2)

# export a compilation database in JSON format.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(fsm LANGUAGES CXX C)

find_package(h5cpp REQUIRED)

add_subdirectory(fsm)
add_subdirectory(examples)
