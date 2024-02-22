cmake_minimum_required (VERSION 3.2.3)
set (CMAKE_CXX_STANDARD 20)
set (CMAKE_EXPORT_COMPILE_COMMANDS OFF)

project ("panels-cpp")

set( SOURCES src/PlanarGeometry.cpp src/PanelMethod.cpp src/ThinWing.cpp src/Wake.cpp)

add_executable(verify_theodorsen ${SOURCES} src/cases/Theodorsen.cpp)
target_compile_options(verify_theodorsen PRIVATE -msse2)

include_directories(Eigen/)
