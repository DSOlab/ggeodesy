noinst_PROGRAMS = testGeodesy \
	transform \
	vincenty

TESTS = testEllipsoid \
        testUnits \
	testAngleNormalization \
	testMeridianArc

check_PROGRAMS = testEllipsoid \
        testUnits \
	testAngleNormalization \
        testTransformations \
	testMeridianArc

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
	-DDEBUG \
        -DCHECK_PRECISION

testGeodesy_SOURCES   = test_geodesy.cpp
testGeodesy_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src 
testGeodesy_LDADD     = $(top_srcdir)/src/libggeodesy.la

transform_SOURCES   = crd_transform_ex.cpp
transform_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src
transform_LDADD     = $(top_srcdir)/src/libggeodesy.la 

vincenty_SOURCES   = vincenty_ex.cpp
vincenty_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src
vincenty_LDADD     = $(top_srcdir)/src/libggeodesy.la 

testEllipsoid_SOURCES   = test_ellipsoid.cpp
testEllipsoid_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src
testEllipsoid_LDADD     = $(top_srcdir)/src/libggeodesy.la

testUnits_SOURCES   = test_units.cpp
testUnits_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src -I.
testUnits_LDADD     = $(top_srcdir)/src/libggeodesy.la

testAngleNormalization_SOURCES   = test_angle_normalization.cc
testAngleNormalization_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src
testAngleNormalization_LDADD     = $(top_srcdir)/src/libggeodesy.la

testTransformations_SOURCES   = test_transformations.cpp
testTransformations_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src -I.
testTransformations_LDADD     = $(top_srcdir)/src/libggeodesy.la

testMeridianArc_SOURCES   = test_meridian_arc.cpp
testMeridianArc_CXXFLAGS  = $(MCXXFLAGS) -I$(top_srcdir)/src -I.
testMeridianArc_LDADD     = $(top_srcdir)/src/libggeodesy.la
