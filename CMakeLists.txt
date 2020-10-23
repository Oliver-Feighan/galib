cmake_minimum_required(VERSION 3.15)

add_definitions(-DARMA_DONT_USE_WRAPPER)

find_package(Armadillo REQUIRED)

if (${ARMADILLO_FOUND})
    message("found armadillo " ${ARMADILLO_VERSION_STRING})
endif ()

message("")

find_package(Boost REQUIRED)
if (${BOOST_FOUND})
    message("found boost " ${Boost_VERSION_STRING})
endif ()

message("")


set(CMAKE_CXX_STANDARD 17)

Project("GALIB")

add_library(LIBGAUS

            src/genes.cpp
            src/genetic_algorithm.cpp

)