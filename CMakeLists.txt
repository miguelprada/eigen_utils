cmake_minimum_required(VERSION 2.8.3)
project(eigen_utils)

find_package(catkin REQUIRED COMPONENTS cmake_modules tf tf_conversions roscpp)
find_package(Eigen REQUIRED)

include_directories(SYSTEM ${Eigen_INCLUDE_DIRS})
include_directories(include ${catkin_INCLUDE_DIRS})

catkin_package(
  INCLUDE_DIRS include ${Eigen_INCLUDE_DIRS}
  LIBRARIES ${PROJECT_NAME}
  CATKIN_DEPENDS
   tf
   tf_conversions
   roscpp
)

add_library(eigen_utils
  src/eigen_utils.cpp
)

target_link_libraries(${PROJECT_NAME}
   ${catkin_LIBRARIES}
)

install(DIRECTORY include/${PROJECT_NAME}/
  DESTINATION ${CATKIN_PACKAGE_INCLUDE_DESTINATION})

install(TARGETS ${PROJECT_NAME}
  DESTINATION ${CATKIN_PACKAGE_LIB_DESTINATION})
 
