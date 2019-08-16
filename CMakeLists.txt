cmake_minimum_required(VERSION 3.7)
project(tcods)
# version
set(serial "1.0")
# basically, we use Release
if (NOT CMAKE_BUILD_TYPE)
   set(CMAKE_BUILD_TYPE "Release")
endif()

####
# need to install following libraries
#   (they need to be found via find_package)
# - SuiteSparseQR
####
# for windows (find packages installed via vcpkg)
# add 
# "-DCMAKE_TOOLCHAIN_FILE=C:/src/vcpkg/scripts/buildsystems/vcpkg.cmake"
# to argument passed to cmake
####

set(SuiteSparse_DIR "C:/src/vcpkg/installed/x64-windows/share/suitesparse")
find_package(SuiteSparse CONFIG REQUIRED)

# Include directories for this project
set(INCLUDE_DIR
)

file(GLOB sources *.cpp)
set(CPP_FILES
    ${sources}
)

set(LINK_LIBS
    ${SuiteSparse_LIBRARIES}
)

# Target for this build
add_executable(tcods ${CPP_FILES})

# Add a include files
target_include_directories(tcods PUBLIC ${INCLUDE_DIR})
# Add external libraries
target_link_libraries(tcods PRIVATE ${LINK_LIBS}
)

set(
    MY_DEBUG_OPTIONS
    -g
)
set(
    MY_RELEASE_OPTIONS
    -O2
)
target_compile_options(tcods PUBLIC "$<$<CONFIG:DEBUG>:${MY_DEBUG_OPTIONS}>")
target_compile_options(tcods PUBLIC "$<$<CONFIG:RELEASE>:${MY_RELEASE_OPTIONS}>")
target_compile_features(tcods PUBLIC cxx_std_11)
target_compile_definitions(tcods PUBLIC _USE_MATH_DEFINES)
