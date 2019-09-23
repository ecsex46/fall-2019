cmake_minimum_required(VERSION 3.5)
project(tinyrender)

# SDL 2
find_package(SDL2 REQUIRED)
if(SDL_FOUND)
    message("SDL Found")
endif()
include_directories(${SDL2_INCLUDE_DIRS})
string(STRIP ${SDL2_LIBRARIES} SDL2_LIBRARIES)

if(NOT WIN32)
    # OpenGL
    find_package(OpenGL REQUIRED)
    if(OPENGL_FOUND)
        message("OpenGL Found")
    endif()
    include_directories(${OPENGL_INCLUDE_DIRS})

    # Threads
    find_package(Threads REQUIRED)
endif()

if(APPLE)
    # Boost (macOS only) for filesystem
    # TODO: Remove this dependency with something lighter in the future
    find_package(Boost 1.65.0 COMPONENTS system filesystem REQUIRED)
    if(BOOST_FOUND)
        message("Boost Found")
    endif()
    include_directories(${Boost_INCLUDE_DIRS})
    link_directories(${Boost_LIBRARY_DIR})
endif()

# GLEW
find_package(GLEW REQUIRED)
if(GLEW_FOUND)
    message("GLEW Found")
endif()
include_directories(${GLEW_INCLUDE_DIRS})

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include_directories("src")
include_directories("externals/")
include_directories("externals/glm/")

file(GLOB_RECURSE srcs src/*.cpp src/*.h)
add_executable(tinyrender ${srcs})

if(WIN32)
    target_link_libraries(tinyrender ${OPENGL_LIBRARIES} ${GLEW_LIBRARIES} SDL2::SDL2 SDL2::SDL2main)
elseif(APPLE)
    target_link_libraries(tinyrender ${Boost_LIBRARIES} ${OPENGL_LIBRARIES} ${GLEW_LIBRARIES} ${SDL2_LIBRARIES})
else()
    target_link_libraries(tinyrender stdc++fs ${OPENGL_LIBRARIES} ${GLEW_LIBRARIES} ${SDL2_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
endif()