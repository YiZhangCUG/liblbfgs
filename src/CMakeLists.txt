# 设定源文件文件夹
aux_source_directory(lib/ LIBLBFGS_SRC)

# 以下部分为库的编译
# 注意目标名必须唯一 所以不能直接生成相同名称的动态库与静态库
# 注意此处不必为目标名称添加lib前缀和相应后缀，cmake会自行添加
add_library(lbfgs SHARED ${LIBLBFGS_SRC})
# 首先添加静态库的生成命令
add_library(lbfgs_static STATIC ${LIBLBFGS_SRC})
# 设置静态库的输出名称从而获得与动态库名称相同的静态库
set_target_properties(lbfgs_static PROPERTIES OUTPUT_NAME "lbfgs")
# 设置输出目标属性以同时输出动态库与静态库
set_target_properties(lbfgs PROPERTIES CLEAN_DIRECT_OUTPUT 1)
set_target_properties(lbfgs_static PROPERTIES CLEAN_DIRECT_OUTPUT 1)
# 设置动态库的版本号
set_target_properties(lbfgs PROPERTIES VERSION 1.0 SOVERSION 1.0)
# 设置库文件的输出地址
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
# 库的安装命令
install(TARGETS lbfgs lbfgs_static
	LIBRARY DESTINATION lib
	ARCHIVE DESTINATION lib)
# 头文件安装命令
install(FILES lib/lbfgs.h DESTINATION include)

# 以下部分为例子程序的编译
# 设置可执行文件的输出地址
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

# 添加可执行文件 命令行
add_executable(lbfgs_sample sample/sample.c)
# 为安装文件添加动态库的搜索地址
set_target_properties(lbfgs_sample PROPERTIES INSTALL_RPATH "/usr/local/lib")
# 链接动态库
target_link_libraries(lbfgs_sample PUBLIC lbfgs)

# 添加可执行文件 命令行
add_executable(lbfgs_sample2 sample/sample2.cpp)
# 为安装文件添加动态库的搜索地址
set_target_properties(lbfgs_sample2 PROPERTIES INSTALL_RPATH "/usr/local/lib")
# 链接动态库
target_link_libraries(lbfgs_sample2 PUBLIC lbfgs)