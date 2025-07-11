cmake_minimum_required(VERSION 3.10)
project(CaWoSched)

# Specify C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")


# Add the source files
add_library(MultiProcessor STATIC lib/multipleProcessors/multiple_processors.cpp)

# Add executables
add_executable(multi_machine_scheduler app/main_multiple_processors.cpp)
target_link_libraries(multi_machine_scheduler PRIVATE MultiProcessor)

# Include the current directory for headers
target_include_directories(MultiProcessor PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_include_directories(MultiProcessor PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/lib/tools)

# Optionally, add the include directory for executables as well
target_include_directories(multi_machine_scheduler PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/lib/tools)
