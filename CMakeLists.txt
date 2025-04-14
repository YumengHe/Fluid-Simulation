cmake_minimum_required(VERSION 3.10)
project(FluidSimulation LANGUAGES CXX)

# Use C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set source files
file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS
    src/*.cpp
)

# Include headers if needed
include_directories(${CMAKE_SOURCE_DIR}/include)

# Create executable
add_executable(FluidSimulation ${SOURCES} "src/particle.h" "src/particleSystem.h" "src/particleSystem.cpp")
