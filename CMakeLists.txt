cmake_minimum_required(VERSION 2.6)
project(egg)

if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "No build type selected, default to Release (-DCMAKE_BUILD_TYPE=...)")
    set(CMAKE_BUILD_TYPE "Release")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")
find_package(phypp)

set(EGG_SHARE_DIR ${CMAKE_INSTALL_PREFIX}/share/egg)
install(FILES
    share/ir_lib_ce01.fits share/ir_lib_cs15.fits
    share/mass_func_candels.fits share/opt_lib_fast.fits
    DESTINATION ${EGG_SHARE_DIR} COMPONENT data)
install(DIRECTORY share/filter-db
    DESTINATION ${EGG_SHARE_DIR} COMPONENT data)
install(DIRECTORY share/skymaker-templates
    DESTINATION ${EGG_SHARE_DIR} COMPONENT data)
install(DIRECTORY share/psfs
    DESTINATION ${EGG_SHARE_DIR} COMPONENT data)

add_definitions(-DEGG_SHARE_DIR="${EGG_SHARE_DIR}")
add_definitions(-DFILTER_DB_DIR="${EGG_SHARE_DIR}/filter-db")

include_directories(${PHYPP_INCLUDE_DIRS})

message(STATUS ${PHYPP_INCLUDE_DIRS})

add_executable(egg-gencat src/egg-gencat.cpp)
target_link_libraries(egg-gencat ${PHYPP_LIBRARIES})
install(TARGETS egg-gencat DESTINATION bin)

add_executable(egg-genmap src/egg-genmap.cpp)
target_link_libraries(egg-genmap ${PHYPP_LIBRARIES})
install(TARGETS egg-genmap DESTINATION bin)

add_executable(egg-gennoise src/egg-gennoise.cpp)
target_link_libraries(egg-gennoise ${PHYPP_LIBRARIES})
install(TARGETS egg-gennoise DESTINATION bin)

add_executable(egg-buildmf src/egg-buildmf.cpp)
target_link_libraries(egg-buildmf ${PHYPP_LIBRARIES})
install(TARGETS egg-buildmf DESTINATION bin)

add_executable(egg-2skymaker src/egg-2skymaker.cpp)
target_link_libraries(egg-2skymaker ${PHYPP_LIBRARIES})
install(TARGETS egg-2skymaker DESTINATION bin)

add_executable(egg-getsed src/egg-getsed.cpp)
target_link_libraries(egg-getsed ${PHYPP_LIBRARIES})
install(TARGETS egg-getsed DESTINATION bin)
