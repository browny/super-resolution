cmake_minimum_required(VERSION 2.8)

INCLUDE_DIRECTORIES(../include/)
ADD_EXECUTABLE(main main.cpp)

TARGET_LINK_LIBRARIES(main /usr/lib/libX11.so
    /usr/lib/libpthread.so)

FILE(COPY ${PROJECT_SOURCE_DIR}/rc/raw.bmp DESTINATION ${PROJECT_BINARY_DIR}/bin)

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
