##
##  The library will be called: (lib)"geodesy".
## ----------------------------------------------
##
lib_LTLIBRARIES = libgeodesy.la

##
##  Only use the following flags for debuging purposes.
## ----------------------------------------------------
##
libgeodesy_la_CXXFLAGS = \
-std=c++17 \
-g \
-pg \
-Wall \
-Wextra \
-Werror \
-pedantic \
-W \
-Wshadow \
-Winline \
-Wdisabled-optimization \
-march=native \
-DDEBUG

##
##  Header files (distributed) installed in /$(includedir)/$(package name).
## ------------------------------------------------------------------------
##
pkginclude_HEADERS = \
geoconst.hpp \
ellipsoid_core.hpp \
ellipsoid.hpp \
geodesy.hpp \
units.hpp \
meridian_arc_length.hpp \
car2ell.hpp \
car2top.hpp \
ell2car.hpp \
vincenty.hpp \
great_circle.hpp \
trnsfdtls.hpp

##
##  Source files (distributed).
## ----------------------------
##
dist_libgeodesy_la_SOURCES = \
top2daz.cpp \
top2car.cpp \
vincenty_direct.cpp \
vincenty_inverse.cpp \
karney_direct.cpp \
karney_inverse.cpp \
trnsfdtls.cpp \
meridian_arc_length.cpp \
great_circle.cpp \
wgs84_to_pz90.cpp
