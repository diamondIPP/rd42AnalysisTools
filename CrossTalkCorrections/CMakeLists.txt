#cmake_minimum_required(VERSION 3.8)
cmake_minimum_required(VERSION 3.5)
project(CrossTalkCorrections)

set(CMAKE_CXX_STANDARD 11)

set(SOURCE_FILES FeedThroughCorrection.cpp)
add_executable(CrossTalkCorrections ${SOURCE_FILES})