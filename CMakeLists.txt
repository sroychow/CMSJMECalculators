cmake_minimum_required(VERSION 3.11)

if(NOT DEFINED CMSJMECALCULATORS_VERSION)
  file(READ "${CMAKE_CURRENT_SOURCE_DIR}/version.txt" CMSJMECALCULATORS_VERSION)
endif()
string(REPLACE "." ";" VERSION_SPLIT ${CMSJMECALCULATORS_VERSION})
list(GET VERSION_SPLIT 0 SPLIT_VERSION_MAJOR)
list(GET VERSION_SPLIT 1 SPLIT_VERSION_MINOR)

project(CMSJMECalculators VERSION ${SPLIT_VERSION_MAJOR}.${SPLIT_VERSION_MINOR} LANGUAGES CXX)

if(SKBUILD) # add python package name, for scikit-build
  set(PKG_LIB_INSTALL "CMSJMECalculators/lib")
  set(PKG_HDR_INSTALL "CMSJMECalculators/include")
  set(PKG_CMK_INSTALL "CMSJMECalculators/cmake")
else()
  set(PKG_LIB_INSTALL "lib")
  set(PKG_HDR_INSTALL "include")
  set(PKG_CMK_INSTALL "cmake")
endif()

find_package(ROOT REQUIRED COMPONENTS MathCore Hist ROOTVecOps Physics)

install(DIRECTORY interface/ DESTINATION ${PKG_HDR_INSTALL})

file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/jetclasses_filenames.txt FILENAMES_CMSJET_INC REGEX ".*/interface/.*")
file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/jetclasses_filenames.txt FILENAMES_CMSJET_SRC REGEX ".*/src/.*")
string(REGEX REPLACE "[^;]+/interface/([^;]+)" "${CMAKE_CURRENT_BINARY_DIR}/CMSJet/include/\\1" INC_CMSJET "${FILENAMES_CMSJET_INC}")
string(REGEX REPLACE "[^;]+/src/([^;]+)" "${CMAKE_CURRENT_BINARY_DIR}/CMSJet/src/\\1" SRC_CMSJET "${FILENAMES_CMSJET_SRC}")
add_custom_target(CMSJetSources COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/getjetclasses.sh COMMENT "Copying jet classes from CMSSW" BYPRODUCTS ${INC_CMSJET} ${SRC_CMSJET})
install(FILES ${INC_CMSJET} DESTINATION ${PKG_HDR_INSTALL})

add_library(CMSJMECalculators SHARED ${SRC_CMSJET} src/JMESystematicsCalculators.cc)
add_dependencies(CMSJMECalculators CMSJetSources)
target_include_directories(CMSJMECalculators PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/interface> $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/CMSJet/include> $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/CMSJet/src> $<INSTALL_INTERFACE:${PKG_HDR_INSTALL}>)
target_link_libraries(CMSJMECalculators ROOT::Hist ROOT::ROOTVecOps ROOT::Physics)
install(TARGETS CMSJMECalculators EXPORT CMSJMECalculators-targets LIBRARY DESTINATION ${PKG_LIB_INSTALL})

install(EXPORT CMSJMECalculators-targets FILE CMSJMECalculators-targets.cmake DESTINATION ${PKG_CMK_INSTALL})
include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/Config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/CMSJMECalculatorsConfig.cmake"
  INSTALL_DESTINATION "${PKG_CMK_INSTALL}"
  NO_SET_AND_CHECK_MACRO NO_CHECK_REQUIRED_COMPONENTS_MACRO
  )
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/CMSJMECalculatorsConfigVersion.cmake"
  VERSION "${CMSJMECalculators_VERSION_MAJOR}.${CMSJMECalculators_VERSION_MINOR}"
  COMPATIBILITY AnyNewerVersion)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/CMSJMECalculatorsConfig.cmake DESTINATION ${PKG_CMK_INSTALL})
