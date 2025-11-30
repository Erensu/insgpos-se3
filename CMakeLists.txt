cmake_minimum_required(VERSION 3.16)
project(insgpos-se3)

add_library(insgpos-se3 insgpos-se3.c insgpos-se3.h)

add_executable(test_insgpos_se3 test_insgpos_se3.c)
target_link_libraries(test_insgpos_se3 insgpos-se3)








