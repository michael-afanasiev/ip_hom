INCLUDE_DIRECTORIES(
  "${CMAKE_CURRENT_SOURCE_DIR}/include"
  "/usr/local/Cellar/ipopt/3.12.1_3/include/coin"
  "/usr/local/include"
  )
LINK_DIRECTORIES(
  "/usr/local/lib"
  )

SET(cxxbind_SRCS
  src/hs071_main.cpp
  src/hs071_nlp.cpp
  )

ADD_EXECUTABLE(OptTest ${cxxbind_SRCS})
TARGET_LINK_LIBRARIES(OptTest ipopt)
