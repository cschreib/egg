cmake_minimum_required(VERSION 2.6)
project(gencat)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")
find_package(phypp)

set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}")

add_phypp_target(gencat.cpp)
add_phypp_target(genmf.cpp)
add_phypp_target(make_skymaker.cpp)
