cmake_minimum_required(VERSION 3.1)
project(make-it-stand)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl settings
option(LIBIGL_WITH_OPENGL              "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW         "Use GLFW"           ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI   "Use ImGui"          ON)

find_package(LIBIGL REQUIRED QUIET)

# include directories
include_directories("${PROJECT_SOURCE_DIR}/include/")
include_directories(${LIBIGL_INCLUDE_DIR})

# include source code
file(GLOB H_FILES ${PROJECT_SOURCE_DIR}/include/*.h)
file(GLOB CPP_FILES ${PROJECT_SOURCE_DIR}/src/*.cpp)

add_executable(${PROJECT_NAME} ${PROJECT_SOURCE_DIR}/main.cpp ${H_FILES} ${CPP_FILES})
target_link_libraries(${PROJECT_NAME} igl::core igl::opengl igl::opengl_glfw igl::opengl_glfw_imgui)
