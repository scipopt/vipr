cmake_minimum_required(VERSION 3.11.0)

# set the project name and version
set(VIPR_VERSION_MAJOR 1)
set(VIPR_VERSION_MINOR 1)
project(VIPR VERSION ${VIPR_VERSION_MAJOR}.${VIPR_VERSION_MINOR})

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# set function visibility default to hidden
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN 1)

# path to find modules
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

# FindBoost is removed in version cmake 3.30
if(POLICY CMP0167)
	cmake_policy(SET CMP0167 NEW)
endif()

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# find GMP
find_package(GMP REQUIRED)
include_directories(${GMP_INCLUDE_DIRS})
include_directories(${PROJECT_BINARY_DIR})
set(libs ${libs} ${GMP_LIBRARIES})

# add executables
add_executable(viprttn viprttn.cpp)
add_executable(vipr2html vipr2html.cpp)
add_executable(viprchk viprchk.cpp)
add_executable(viprchk_parallel viprchk_parallel.cpp)

# find TBB
if(WIN32)
      find_package(TBB 2021 COMPONENTS tbb tbbmalloc REQUIRED)
   else()
      # Utilities to automatically download missing dependencies
      include(cmake/Dependencies.cmake)

      find_or_download_package(
         NAME TBB
         VERSION 2021
        URL https://github.com/oneapi-src/oneTBB/archive/refs/tags/v2021.7.0.tar.gz
         COMPONENTS tbb tbbmalloc
         URL_HASH SHA256=2cae2a80cda7d45dc7c072e4295c675fff5ad8316691f26f40539f7e7e54c0cc
         CONFIGURE_ARGS -DTBB_TEST=OFF -DTBB_EXAMPLES=OFF -DTBB4PY_BUILD=OFF)
         set(libs ${libs} ${TBB_LIBRARIES})
endif()



# option to install viprcomp
option(VIPRCOMP "Use viprcomp" ON)



target_link_libraries(viprchk ${libs})
target_link_libraries(viprchk_parallel ${libs})
target_link_libraries(viprchk_parallel TBB::tbb)


set(BOOST_MIN_VERSION 1.71)
if(VIPRCOMP)
	find_package(Boost ${BOOST_MIN_VERSION} REQUIRED)
	if(Boost_FOUND)
		include_directories(${Boost_INCLUDE_DIRS})
	endif()
	# Only install viprcomp if working SoPlex is found
	find_package(ZLIB)
	if(ZLIB_FOUND)
		find_package(SOPLEX)
		if(SOPLEX_FOUND)
			set(VIPR_HAVE_SOPLEX 1)
			# include ZLIB
			include_directories(${ZLIB_INCLUDE_DIRS})
			set(libs ${libs} ${ZLIB_LIBRARIES})

			# include SoPlex
			include_directories(${SOPLEX_INCLUDE_DIRS})
			set(libs ${libs} ${SOPLEX_LIBRARIES})

			# add viprcomp target and link
			add_executable(viprcompseq viprcompseq.cpp)
			add_executable(viprincomp incompletify.cpp)
			add_definitions(-DSOPLEX_WITH_GMP)
			target_link_libraries(viprcompseq ${libs})
			target_link_libraries(viprincomp ${libs})
         target_link_libraries(viprcompseq ${libs})
#         target_link_libraries(viprcompseq  TBB::tbb)
         message(STATUS "Soplex found.")
		else()
			message(STATUS "viprcomp not installed, because SoPlex could not be found.")
			set(VIPR_HAVE_SOPLEX 0)
		endif()
	else()
		message(STATUS "viprcomp not installed, because ZLIB could not be found (required by SoPlex).")
	endif()
endif()


configure_file("${PROJECT_SOURCE_DIR}/CMakeConfig.hpp.in"
               "${PROJECT_BINARY_DIR}/CMakeConfig.hpp")
