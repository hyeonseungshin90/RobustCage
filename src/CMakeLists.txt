cmake_minimum_required(VERSION 3.14)
# Add .cmake files in directory 'cmake'
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
cmake_policy(SET CMP0074 NEW)

project(CageGenerator CXX)

# C++ standard 17.
# in C++17 the standard allocators have been updated to respect type's alignment.
set(CMAKE_CXX_STANDARD 17)

# OpenMP
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC" OR
    "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    find_package(OpenMP)
    if(NOT OPENMP_FOUND)
        message(FATAL_ERROR "OpenMP Not Found!")
    endif()
endif()

# Third Party Libraries
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/ThirdPartyLib)

# CGAL
set(CGAL_Boost_USE_STATIC_LIBS ON)
find_package(CGAL)
if(CGAL_FOUND)
    MESSAGE(STATUS "CGAL INCLUDE DIRS:${CGAL_INCLUDE_DIRS}")
else()
    MESSAGE(FATAL "Can't find CGAL.")
endif()

# Boost
# set(Boost_DEBUG ON)
set(Boost_USE_STATIC_LIBS ON)
find_package(Boost 1.79.0
    REQUIRED
    COMPONENTS filesystem json
)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    MESSAGE(STATUS "-------------------------------------------")
    MESSAGE(STATUS "Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}.")
    MESSAGE(STATUS "Boost_LIBRARY_DIRS = ${Boost_LIBRARY_DIRS}.")
    MESSAGE(STATUS "Boost_LIBRARIES = ${Boost_LIBRARIES}.")
    MESSAGE(STATUS "Boost_LIB_VERSION = ${Boost_LIB_VERSION}.")
    MESSAGE(STATUS "-------------------------------------------")
endif()

# Our source files
file(GLOB_RECURSE SRCS
  "CageInitializer/*"
  "CageSimplifier/*"
  "Geometry/*"
  "Mesh/*"
  "Graph/*"
  "Utils/*"
)

# build library
add_library(lib${PROJECT_NAME} STATIC
    ${SRCS}
    Config.hh
    CageGenerator.hh
    CageGenerator.cc
)
target_include_directories(lib${PROJECT_NAME} PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${EIGEN3_INCLUDE_DIR}
    ${AUTODIFF_INCLUDE_DIR}
    ${OpenMesh_INCLUDE_DIR}
)
target_link_libraries(lib${PROJECT_NAME} PUBLIC
    mesh_generator_lib
    spdlog::spdlog
    autodiff::autodiff
    indirectPredicates
    geometry-central
    ${OpenMesh_LIBRARIES}
    CGAL::CGAL
    OpenMP::OpenMP_CXX
)

target_compile_definitions(lib${PROJECT_NAME} PUBLIC _USE_MATH_DEFINES NOMINMAX)

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
    target_compile_options(lib${PROJECT_NAME} PRIVATE "/wd4819")
endif()

# options to control program behavior
target_compile_definitions(lib${PROJECT_NAME} PUBLIC
    #CAGE_DEBUG
    #OUTPUT_VMESH
    #OUTPUT_MIDDLE_RESULT
    #CHECK_MANIFOLD
    #USE_TREE_SEARCH
)

add_executable(exe${PROJECT_NAME} main.cpp)

target_include_directories(exe${PROJECT_NAME} PUBLIC
    ${Boost_INCLUDE_DIRS}
)

target_link_libraries(exe${PROJECT_NAME} PUBLIC 
    lib${PROJECT_NAME}
    ${Boost_LIBRARIES}
)

option(CAGE_BUILD_TESTS "Build cage quality regression tests" OFF)
if(CAGE_BUILD_TESTS)
    enable_testing()
    add_executable(cageQualityPolishTests tests/QualityPolishTests.cpp)
    target_link_libraries(cageQualityPolishTests PRIVATE lib${PROJECT_NAME} ${Boost_LIBRARIES})
    if(WIN32)
        add_custom_command(TARGET cageQualityPolishTests POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${CMAKE_CURRENT_SOURCE_DIR}/ThirdPartyLib/gmp/lib/libgmp-10.dll"
                "${CMAKE_CURRENT_SOURCE_DIR}/ThirdPartyLib/mpfr/lib/libmpfr-4.dll"
                "$<TARGET_FILE_DIR:cageQualityPolishTests>")
    endif()
    add_test(NAME quality_polish COMMAND cageQualityPolishTests)
    set_tests_properties(quality_polish PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
    add_executable(cageBoundaryRailGeodesicTests tests/BoundaryRailGeodesicTests.cpp)
    target_link_libraries(cageBoundaryRailGeodesicTests PRIVATE lib${PROJECT_NAME} ${Boost_LIBRARIES})
    if(WIN32)
        add_custom_command(TARGET cageBoundaryRailGeodesicTests POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${CMAKE_CURRENT_SOURCE_DIR}/ThirdPartyLib/gmp/lib/libgmp-10.dll"
                "${CMAKE_CURRENT_SOURCE_DIR}/ThirdPartyLib/mpfr/lib/libmpfr-4.dll"
                "$<TARGET_FILE_DIR:cageBoundaryRailGeodesicTests>")
    endif()
    add_test(NAME boundary_rail_geodesics COMMAND cageBoundaryRailGeodesicTests)
    set_tests_properties(boundary_rail_geodesics PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
    add_executable(cageBoundaryRailNumericsTests tests/BoundaryRailNumericsTests.cpp)
    target_link_libraries(cageBoundaryRailNumericsTests PRIVATE lib${PROJECT_NAME} ${Boost_LIBRARIES})
    if(WIN32)
        add_custom_command(TARGET cageBoundaryRailNumericsTests POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${CMAKE_CURRENT_SOURCE_DIR}/ThirdPartyLib/gmp/lib/libgmp-10.dll"
                "${CMAKE_CURRENT_SOURCE_DIR}/ThirdPartyLib/mpfr/lib/libmpfr-4.dll"
                "$<TARGET_FILE_DIR:cageBoundaryRailNumericsTests>")
    endif()
    add_test(NAME boundary_rail_numerics COMMAND cageBoundaryRailNumericsTests)
    set_tests_properties(boundary_rail_numerics PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
    add_executable(cageBoundaryRailSearchTests tests/BoundaryRailSearchTests.cpp)
    target_link_libraries(cageBoundaryRailSearchTests PRIVATE lib${PROJECT_NAME} ${Boost_LIBRARIES})
    if(WIN32)
        add_custom_command(TARGET cageBoundaryRailSearchTests POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different
                "${CMAKE_CURRENT_SOURCE_DIR}/ThirdPartyLib/gmp/lib/libgmp-10.dll"
                "${CMAKE_CURRENT_SOURCE_DIR}/ThirdPartyLib/mpfr/lib/libmpfr-4.dll"
                "$<TARGET_FILE_DIR:cageBoundaryRailSearchTests>")
    endif()
    add_test(NAME boundary_rail_search COMMAND cageBoundaryRailSearchTests)
    set_tests_properties(boundary_rail_search PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=1")
endif()
