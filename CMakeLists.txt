cmake_minimum_required(VERSION 3.15.2)

# 设置工程名称和语言
project(LIBLBFGS VERSION 1.1.0)

if(${CMAKE_HOST_SYSTEM_NAME} STREQUAL "Linux")
	message(STATUS "Platform: " ${CMAKE_HOST_SYSTEM_NAME})
	set(CMAKE_INSTALL_PREFIX "/usr/local")
elseif (${CMAKE_HOST_SYSTEM_NAME} STREQUAL "Darwin")
	message(STATUS "Platform: " ${CMAKE_HOST_SYSTEM_NAME})
	set(CMAKE_INSTALL_PREFIX "/usr/local")
elseif (${CMAKE_HOST_SYSTEM_NAME} STREQUAL "Windows")
	message(STATUS "Platform: " ${CMAKE_HOST_SYSTEM_NAME})
	set(CMAKE_C_COMPILER gcc)
	set(CMAKE_CXX_COMPILER g++)
	# CMake中的路径分隔符统一为UNIX格式 不需要使用反斜杠 否则会报警
	set(CMAKE_INSTALL_PREFIX "D:/Library")
else()
	message(STATUS "Platform: " ${CMAKE_HOST_SYSTEM_NAME})
	set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
message(STATUS "Install prefix: " ${CMAKE_INSTALL_PREFIX})

# 添加源文件地址
add_subdirectory(src/)

# 添加编译选项
option(HAVE_CONFIG_H "use config.h" ON)
option(_MSC_VER "use __msc_ver" OFF)
option(USE_SSE "use see" OFF)
# 加入一个头文件配置，让cmake对源码进行操作
configure_file(
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_SOURCE_DIR}/src/lib/config.h"
  )
