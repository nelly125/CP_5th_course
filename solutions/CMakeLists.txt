cmake_minimum_required(VERSION 3.21)
project(solutions)

set(CMAKE_CXX_STANDARD 17)

cmake_minimum_required(VERSION 2.6)

# Locate GTest
find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

file(GLOB_RECURSE SOURCES RELATIVE ${CMAKE_SOURCE_DIR} "../src/*.cpp")


# Link runTests with what we want to test and the GTest and pthread library
add_executable(runTests limiters.cpp sod_test.cpp)
find_package(Threads)
target_link_libraries(runTests ${GTEST_LIBRARIES} pthread Godunov_solver)