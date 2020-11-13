noinst_PROGRAMS = testGeodTest

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
