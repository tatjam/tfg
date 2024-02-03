cmake_minimum_required (VERSION 3.2.3)
set (CMAKE_CXX_STANDARD 20)
set (CMAKE_EXPORT_COMPILE_COMMANDS OFF)

project ("panels-cpp")
add_executable(${PROJECT_NAME})
target_sources(${PROJECT_NAME} PRIVATE src/Main.cpp src/Geometry.cpp)

target_compile_options(${PROJECT_NAME} PRIVATE -msse2)

include_directories(Eigen/)
