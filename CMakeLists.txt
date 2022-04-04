cmake_minimum_required(VERSION 3.2 FATAL_ERROR)

# ソリューションエクスプローラーにCMake用のフォルダを別に用意する
set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set_property(GLOBAL PROPERTY PREDEFINED_TARGETS_FOLDER "CMake")

# edubpt本体プロジェクト
project(edubpt VERSION 1.0)
set(CMAKE_CXX_STANDARD 11)
add_executable(edubpt main.cpp)

# インストール
install (TARGETS edubpt DESTINATION ${CMAKE_SOURCE_DIR}/bin)