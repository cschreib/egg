cmake_minimum_required(VERSION 2.6)
project(ifni)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")
find_package(phypp)

set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}")

add_phypp_target(ifni-gencat.cpp)
add_phypp_target(ifni-genmf.cpp)
add_phypp_target(ifni-2skymaker.cpp)
