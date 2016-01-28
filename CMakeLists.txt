CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )
PROJECT( heat )

SET( SRCS

   heat.cpp
   # more *.cpp files follow here
)

OPTION( BIG_HELLO   "enables/disables big letter mode (off by default)"                                          OFF )
OPTION( ENABLE_ANSI "enables/disables addition of '-ansi -pedantic' flags for gnu c++ compiler (off by default)" OFF )

SET(CMAKE_CXX_COMPILER mpicxx)
target_link_libraries(mpi)

IF( BIG_HELLO )
   ADD_DEFINITIONS( -DBIG_HELLO )
ENDIF()

IF( CMAKE_COMPILER_IS_GNUCC )
      SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wextra -Wshadow -std=c++11 " )
ENDIF()

ADD_EXECUTABLE( heat ${SRCS} )