cmake_minimum_required(VERSION 3.9)
project(integration_gsl_wrapper
        VERSION 1.0
        LANGUAGES CXX)

# standard of c++ specification
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED YES)
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra -Wno-psabi -Wno-deprecated-declarations -Wold-style-cast")

# add gsl
find_package(GSL REQUIRED)

# Set the output folder where program will be created
 set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR})

# The following folder will be included
include_directories("${PROJECT_SOURCE_DIR}")
include_directories(include)

add_executable(test_contour ${PROJECT_SOURCE_DIR}/src/test_contour.cpp)
target_link_libraries(test_contour GSL::gsl GSL::gslcblas)

install ( TARGETS test_contour DESTINATION ${CMAKE_INSTALL_PREFIX} )