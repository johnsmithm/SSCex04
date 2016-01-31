CMAKE_MINIMUM_REQUIRED( VERSION 2.6 )
PROJECT( heat )

SET( SRCS

   heat.cpp
   # more *.cpp files follow here
)

OPTION( ARCHITECTURE_OPTIMIZATION "enables/disables addition of '-march=native' flags for gnu c++ compiler (off by default)" OFF )
OPTION( ENABLE_WERROR "enables/disables addition of '-Werror' flags for gnu c++ compiler (off by default)" OFF )

SET(CMAKE_CXX_COMPILER mpicxx)
target_link_libraries(mpi)


IF( CMAKE_COMPILER_IS_GNUCC )
      SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Wall -Wextra -Wshadow -std=c++11 " )
      IF ( ARCHITECTURE_OPTIMIZATION )
      	SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native " )
      ENDIF()
      IF ( ENABLE_WERROR )
      	SET( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror " )
      ENDIF()
ENDIF()

ADD_EXECUTABLE( heat ${SRCS} )