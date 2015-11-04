cmake_minimum_required(VERSION 2.6)
project(ex1)

add_definitions(-std=gnu++0x)

add_executable(ex1 main.cpp)

install(TARGETS ex1 RUNTIME DESTINATION bin)
