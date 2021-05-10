noinst_PROGRAMS = testGeodTest testGeodTime

#TESTS =

#check_PROGRAMS =

MCXXFLAGS = \
-std=c++17 \
-g \
-pg \
-Wno-unused-variable \
-Wall \
-Wextra \
-Werror \
-pedantic \
-W \
-Wshadow \
-Wdisabled-optimization \
-DDEBUG

testGeodTest_SOURCES   = test_geodtest.cpp
testGeodTest_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src
testGeodTest_LDADD     = $(top_srcdir)/src/libggeodesy.la

testGeodTime_SOURCES   = test_geodtime.cpp
testGeodTime_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src
testGeodTime_LDADD     = $(top_srcdir)/src/libggeodesy.la
