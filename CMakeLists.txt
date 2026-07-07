cmake_minimum_required(VERSION 3.15)
project(celestial C)

set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

add_library(celestial STATIC
    src/roots_solver.c
    src/ivp_solver.c
)

target_include_directories(celestial PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_link_libraries(celestial PUBLIC m)