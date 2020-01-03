cmake_minimum_required(VERSION 3.9)

# export a compilation database in JSON format.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(marlin LANGUAGES CXX C)

add_subdirectory(marlin)
add_subdirectory(examples)
