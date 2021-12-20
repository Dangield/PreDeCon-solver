cmake_minimum_required(VERSION 3.10)
set (CMAKE_CXX_STANDARD 14)

# set the project name
project(Clustering)

# output dirs
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# for matplotlib-cpp
find_package(Python3 COMPONENTS Interpreter Development NumPy REQUIRED)

add_library(matplotlib_cpp INTERFACE)
target_include_directories(matplotlib_cpp
  INTERFACE
    ../matplotlib-cpp
)
target_link_libraries(matplotlib_cpp INTERFACE Python3::Python Python3::Module Python3::NumPy)

# build libraries
set(SAMPLE_SOURCES src/sample.cpp)
add_library(samples ${SAMPLE_SOURCES})
target_include_directories(samples PUBLIC include)

set(PREDECON_SOURCES src/predecon.cpp)
add_library(predecon ${PREDECON_SOURCES})
target_include_directories(predecon PUBLIC include)
target_link_libraries(predecon samples matplotlib_cpp)

# add the executable
set(SOURCES src/main.cpp)
add_executable(main ${SOURCES})
target_include_directories(main PUBLIC include)
target_link_libraries(main predecon)
