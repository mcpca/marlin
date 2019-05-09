cmake_minimum_required(VERSION 3.2)

# export a compilation database in JSON format.
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(fsm LANGUAGES CXX C)

add_subdirectory(fsm)
add_subdirectory(examples)
