cmake_minimum_required(VERSION 3.21)
project(Godunov_solver)

set(CMAKE_CXX_STANDARD 17)


file(GLOB_RECURSE SOURCES RELATIVE ${CMAKE_SOURCE_DIR} "src/*.cpp")

message(${SOURCES})

add_executable(Godunov_solver main.cpp ${SOURCES} ./src/system/system_helper.cpp)

#add_executable(Godunov_solver main.cpp )
