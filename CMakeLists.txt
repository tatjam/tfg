cmake_minimum_required (VERSION 3.2.3)
set (CMAKE_CXX_STANDARD 20)
set (CMAKE_EXPORT_COMPILE_COMMANDS OFF)

project ("panels-cpp")

set( SOURCES src/SheetGeometry.cpp src/PanelMethod.cpp src/ThinWing.cpp src/Wake.cpp)

add_executable(verify_theodorsen ${SOURCES} src/cases/Theodorsen.cpp)
add_executable(verify_wagner ${SOURCES} src/cases/Wagner.cpp)
add_executable(verify_wake ${SOURCES} src/cases/WakeGeneration.cpp)
add_executable(verify_steady_rect ${SOURCES} src/cases/SteadyRectangle.cpp)
add_executable(verify_clalpha ${SOURCES} src/cases/clAlpha.cpp)
target_compile_options(verify_theodorsen PRIVATE -msse2)
target_compile_options(verify_wagner PRIVATE -msse2)
target_compile_options(verify_wake PRIVATE -msse2)

include_directories(Eigen/)
