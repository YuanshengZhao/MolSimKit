cmake_minimum_required(VERSION 3.12)
set(CMAKE_CXX_STANDARD 17)

add_compile_options(-Wall -O3)

project(md)
if(APPLE)
    link_directories(/usr/local/Cellar/openblas/0.3.21/lib/)
elseif(UNIX)
    link_directories(./lib/OpenBLAS/lib)
endif()

file(GLOB SRCFILES CONFIGURE_DEPENDS "src/*.h" "src/*.cpp")
add_executable(msk.exe ${SRCFILES} main.cpp)
target_link_libraries(msk.exe openblas)

