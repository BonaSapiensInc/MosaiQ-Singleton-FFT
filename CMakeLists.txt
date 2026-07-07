cmake_minimum_required(VERSION 3.20)

project(MosaiQ_Singleton_FFT VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

add_library(mosaiq_singleton_fft STATIC
    src/SingletonFFT.cpp
)

add_library(mosaiq::fft ALIAS mosaiq_singleton_fft)

target_include_directories(mosaiq_singleton_fft
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
        $<INSTALL_INTERFACE:include>
)

target_compile_features(mosaiq_singleton_fft PUBLIC cxx_std_20)

if(MSVC)
    target_compile_options(mosaiq_singleton_fft PRIVATE /W4 /WX)
else()
    target_compile_options(mosaiq_singleton_fft PRIVATE
        -Wall -Wextra -Wpedantic -Werror
    )
endif()

include(FetchContent)

FetchContent_Declare(
    Catch2
    GIT_REPOSITORY https://github.com/catchorg/Catch2.git
    GIT_TAG        v3.5.4
)

FetchContent_MakeAvailable(Catch2)

enable_testing()

add_executable(mosaiq_singleton_fft_tests tests/test_fft.cpp)
target_link_libraries(mosaiq_singleton_fft_tests PRIVATE mosaiq_singleton_fft Catch2::Catch2WithMain)

if(MSVC)
    target_compile_options(mosaiq_singleton_fft_tests PRIVATE /W4 /WX)
else()
    target_compile_options(mosaiq_singleton_fft_tests PRIVATE
        -Wall -Wextra -Wpedantic -Werror
    )
endif()

include(CTest)
include(Catch)
catch_discover_tests(mosaiq_singleton_fft_tests)
