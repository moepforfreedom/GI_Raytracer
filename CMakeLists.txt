cmake_minimum_required(VERSION 3.2)

project(global-illu)

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_STANDARD 14)

set(CMAKE_PREFIX_PATH /media/philip/Data2/qt_linux/5.7/gcc_64)

find_package(Qt5Widgets)
find_package(Threads)

set(SOURCES main.cpp include/octree.cpp include/photonMap.cpp include/meshLoader.cpp include/sceneLoader.cpp include/util.cpp)
set(HEADERS include/gui.h include/image.h include/ray.h include/atmosphere.h include/entities.h include/camera.h include/halton_enum.h include/halton_sampler.h include/raytracer.h include/viewer.h include/octree.h include/photonMap.h include/photon.h include/bbox.h include/material.h include/util.h include/sceneLoader.h include/meshLoader.h include/light.h)

add_executable(global-illu ${SOURCES} ${HEADERS})
include_directories(global-illu include 3rd_party)

option(CLANG_FORMAT_TARGET "automatic clang format" OFF)
if(CLANG_FORMAT_TARGET)
    include(ClangFormat)
    add_clang_format(TARGET global-illu STYLE "LLVM" FILES ${SOURCES} ${HEADERS})
endif(CLANG_FORMAT_TARGET)

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#target_compile_options(global-illu PRIVATE -g)
target_link_libraries(global-illu Qt5::Widgets ${CMAKE_THREAD_LIBS_INIT})


