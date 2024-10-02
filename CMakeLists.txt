CMAKE_MINIMUM_REQUIRED(VERSION 3.15.0)
PROJECT(ibImplicit)

set(MPI_DIR /home/gkluhana/autoibamr/packages/petsc-3.21-2)
#set(MPI_DIR /home/gkluhana/autoibamr/packages/petsc-3.21-2-opt)
set(MPI_INCLUDES ${MPI_DIR}/include)
set(MPI_LIBRARY_DIRS ${MPI_DIR}/lib)

message(STATUS "MPI_DIR: ${MPI_DIR}")
message(STATUS "MPI_INCLUDES: ${MPI_INCLUDES}")
message(STATUS "MPI_LIBRARY_DIRS: ${MPI_LIBRARY_DIRS}")

find_library(MPI_LIBRARY
	NAMES mpi 
	HINTS ${MPI_LIBRARY_DIRS}
)

if(NOT MPI_LIBRARY)
	message(FATAL_ERROR "MPI library not found")
else() 
	message("MPI library found")
endif()


#set(MPI_HOME "/home/gkluhana/autoibamr/packages/petsc-3.21-2")
set(MPI_HOME "/home/gkluhana/autoibamr/packages/petsc-3.21-2-opt")
find_package(MPI REQUIRED)
if (MPI_C_FOUND AND MPI_CXX_FOUND)
	set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
	set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})
else()
	message(FATAL_ERROR "MPI C or C++ compiler not found")
endif()

set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g")
set(CMAKE_BUILD_TYPE Debug)

include_directories(${ibImplicit_SOURCE_DIR} ${ibImplicit_SOURCE_DIR}/include)
set(PETSc_DIR /home/gkluhana/autoibamr/packages/petsc-3.21-2)
#set(PETSc_DIR /home/gkluhana/autoibamr/packages/petsc-3.21-2-opt)

set(PETSc_INCLUDES ${PETSc_DIR}/include)
set(PETSc_LIBRARY_DIRS ${PETSc_DIR}/lib)
message(STATUS "PETSc_DIR: ${PETSc_DIR}")
message(STATUS "PETSc_LIBRARY_DIRS: ${PETSc_LIBRARY_DIRS}")
find_library(PETSc_LIBRARY
	NAMES petsc 
	HINTS ${PETSc_LIBRARY_DIRS}
)

if(NOT PETSc_LIBRARY)
	message(FATAL_ERROR "PETSc library not found")
endif()


include_directories(${PETSc_INCLUDES} )

add_library(ibImplicit SHARED
			src/initialize.cpp src/staggered.cpp src/utilities.cpp src/boundary.cpp
			src/lagrangian.cpp src/stokes.cpp src/multigrid.cpp
) 

target_link_libraries(ibImplicit ${MPI_LIBRARY} ${PETSc_LIBRARY}) 
target_include_directories(ibImplicit PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_compile_options(ibImplicit PUBLIC ${COMPILE_FLAGS})

add_executable(stokes ./examples/stokes.cpp)
add_dependencies(stokes ibImplicit)
target_link_libraries(stokes ${MPI_LIBRARY} ${PETSc_LIBRARY} ibImplicit)



