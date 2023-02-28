cmake_minimum_required(VERSION 3.9.0)

project(JointRec VERSION 1.0)
set(CMAKE_INSTALL_PREFIX "${PROJECT_SOURCE_DIR}/compiled")
set(EXECUTABLE_OUTPUT_PATH "${CMAKE_INSTALL_PREFIX}/bin")
set(LIBRARY_OUTPUT_PATH "${CMAKE_INSTALL_PREFIX}/lib")
set(CMAKE_EXPORT_COMPILE_COMMANDS 1 CACHE BOOL "for clang" FORCE)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_COMPILER "/usr/bin/g++")
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "  -O0  -g ")

find_package(ROOT  CONFIG REQUIRED COMPONENTS Minuit)
include("${ROOT_USE_FILE}")
root_generate_dictionary(LactRecEvent ${PROJECT_SOURCE_DIR}/include/LACTRecEvent.h 
		LINKDEF ${PROJECT_SOURCE_DIR}/include/LinkDef.h)
root_generate_dictionary(KM2AEvent ${PROJECT_SOURCE_DIR}/include/LHEvent.h
		LINKDEF ${PROJECT_SOURCE_DIR}/include/LHEventLinkDef.h)
root_generate_dictionary(KM2ARecEvent ${PROJECT_SOURCE_DIR}/include/KM2ARecEvent.h
		LINKDEF ${PROJECT_SOURCE_DIR}/include/LHRecEventLinkDef.h)

add_library(KM2A SHARED)
target_sources(KM2A PRIVATE KM2AEvent.cxx KM2ARecEvent.cxx ${PROJECT_SOURCE_DIR}/src/G4KM2A_Geometry.cc ${PROJECT_SOURCE_DIR}/src/G4KM2A_RecE.cpp ${PROJECT_SOURCE_DIR}/src/G4KM2A_Reconstruction.cpp ${PROJECT_SOURCE_DIR}/src/KM2ARecEvent.cc ${PROJECT_SOURCE_DIR}/src/LHEvent.cc)
target_include_directories(KM2A PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(KM2A PRIVATE  ${ROOT_LIBRARIES} )


add_executable(JointRec)
target_sources(JointRec PUBLIC ${PROJECT_SOURCE_DIR}/src/JointRec.cpp ${PROJECT_SOURCE_DIR}/src/LACT_Reconstruction.cpp ${PROJECT_SOURCE_DIR}/src/LACT_RunPara.cpp ${PROJECT_SOURCE_DIR}/src/LACT_Telconfig.cpp ${PROJECT_SOURCE_DIR}/src/LACT_TelData.cpp ${PROJECT_SOURCE_DIR}/src/LACTEvent.cpp  ${PROJECT_SOURCE_DIR}/src/LACTRecEvent.cpp LactRecEvent.cxx ${PROJECT_SOURCE_DIR}/src/straux.c )
target_include_directories(JointRec PUBLIC ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(JointRec PRIVATE  ${ROOT_LIBRARIES} KM2A)
