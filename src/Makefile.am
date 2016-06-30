##
##  The library will be called: (lib)"ggeodesy".
## ----------------------------------------------
##
lib_LTLIBRARIES = libggeodesy.la

##
##  Only use the following flags for debuging purposes.
## ----------------------------------------------------
##
libggdatetime_la_CXXFLAGS = \
	-std=c++14 \
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
        ellipsoid.hpp \
        geodesy.hpp \
        car2ell.hpp \
        car2top.hpp \
        ell2car.hpp \
	vincenty.hpp

##
##  Source files (distributed).
## ----------------------------
##
dist_libggeodesy_la_SOURCES = \
	top2daz.cpp
