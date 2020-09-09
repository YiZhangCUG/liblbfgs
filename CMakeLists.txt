cmake_minimum_required(VERSION 3.15.2)

# 设置工程名称和语言
project(LIBLBFGS VERSION 1.1.0)
# 设置安装地址 改用homebrew安装后这一句就不需要了
set(CMAKE_INSTALL_PREFIX /usr/local)
# 添加编译选项
option(HAVE_CONFIG_H "use config.h" ON)
option(_MSC_VER "use __msc_ver" OFF)
option(USE_SSE "use see" OFF)
# 加入一个头文件配置，让cmake对源码进行操作
configure_file(
  "${PROJECT_SOURCE_DIR}/config.h.in"
  "${PROJECT_SOURCE_DIR}/src/lib/config.h"
  )
# 添加源文件地址
add_subdirectory(src/)
