noinst_PROGRAMS = testGeodesy \
	transform \
	vincenty

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

testGeodesy_SOURCES   = test_geodesy.cpp
testGeodesy_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src -L$(top_srcdir)/src
testGeodesy_LDADD     = $(top_srcdir)/src/libggeodesy.la

transform_SOURCES   = crd_transform_ex.cpp
transform_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src -L$(top_srcdir)/src
transform_LDADD     = $(top_srcdir)/src/libggeodesy.la 

vincenty_SOURCES   = vincenty_ex.cpp
vincenty_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src -L$(top_srcdir)/src
vincenty_LDADD     = $(top_srcdir)/src/libggeodesy.la 
