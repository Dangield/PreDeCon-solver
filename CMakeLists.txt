cmake_minimum_required(VERSION 3.10)
set (CMAKE_CXX_STANDARD 14)

# set the project name
project(Clustering)

# output dirs
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# build libraries
set(SAMPLE_SOURCES src/sample.cpp)
add_library(samples ${SAMPLE_SOURCES})
target_include_directories(samples PUBLIC include)

set(DBSCAN_SOURCES src/dbscan.cpp)
add_library(dbscan ${DBSCAN_SOURCES})
target_include_directories(dbscan PUBLIC include)
target_link_libraries(dbscan samples)

set(PREDECON_SOURCES src/predecon.cpp)
add_library(predecon ${PREDECON_SOURCES})
target_include_directories(predecon PUBLIC include)
target_link_libraries(predecon samples)

set(RUN_DBSCAN_SOURCES src/run_dbscan.cpp)
add_executable(run_dbscan ${RUN_DBSCAN_SOURCES})
target_include_directories(run_dbscan PUBLIC include)
target_link_libraries(run_dbscan dbscan)

set(RUN_PREDECON_SOURCES src/run_predecon.cpp)
add_executable(run_predecon ${RUN_PREDECON_SOURCES})
target_include_directories(run_predecon PUBLIC include)
target_link_libraries(run_predecon predecon)

add_executable(DIM512 src/prepare_dim512.cpp)
