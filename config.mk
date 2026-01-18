#-------------------------------------------------------------------------------
#	$Id: config.mk.in 227 2015-08-12 00:25:27Z pwessel $
#
#	Standard Makefile Macro Setup for GMTSAR
#
# Edit this file only, NOT the makefile itself.
#-------------------------------------------------------------------------------
# The purpose of this section is to contain common make macros
# that should be processed by every execution of that utility.
#-------------------------------------------------------------------------------

# POSIX shell.  On some platforms it is not /bin/sh.
SHELL		= /bin/bash

# Installation Directories:
PACKAGE_TARNAME = gmtsar
GMTSARHOME	= /home/ysdong/Software/GMTSAR
ORBITS_DIR	= /usr/local/orbits
TIFF_INC	= /usr/local/include
TIFF_LIB	= /usr/local/lib
LALIBS		= -llapack -lblas
prefix          = /home/ysdong/Software/GMTSAR
exec_prefix     = ${prefix}
bindir		= $(DESTDIR)${exec_prefix}/bin
sharedir	= $(DESTDIR)/home/ysdong/Software/GMTSAR/share/gmtsar

# GMTSAR version number
GMTSAR_VERSION = 6.2.0

# Compilers, if $CC not set in environment
CC		= gcc

# Preprocessing:
CPP		= gcc -E
CPPFLAGS	= $(INCLUDES) $(DEFINES) 
HDF5_CPPFLAGS	= -I/usr/include/hdf5/serial -I/usr/include
HDF5_LDFLAGS	= -L/usr/lib/x86_64-linux-gnu/hdf5/serial 

#-------------------------------------------------------------------------------
#	Math library specification 
#	(Will most of the time be -lm, on Solaris -lsunmath -lm)
#-------------------------------------------------------------------------------
#
LIBS		= -lm 
HDF5_LIBS	= -lm -ldl -lz -lsz -lpthread -lcurl -lcrypto  -lhdf5 -lhdf5_hl -lhdf5_cpp
#
#-------------------------------------------------------------------------------
#	Miscellaneous Standard Utilities
#-------------------------------------------------------------------------------
#
INSTALL		= /usr/bin/install -c
AWK		= mawk
AR		= ar
RANLIB		= ranlib
CSH		= csh
LD		= ld
LN_S		= ln -s
GNUTAR		= :
#
#-------------------------------------------------------------------------------
#	Required directives for GMT includes and library
#-------------------------------------------------------------------------------
GMT_INC		= -I/usr/include/gmt
GMT_LIB		= -L/usr/lib/x86_64-linux-gnu -lgmt
#
#-------------------------------------------------------------------------------
#	Required directives for GMTSAR library
#-------------------------------------------------------------------------------
GMTSAR		= -L/home/ysdong/Software/GMTSAR/gmtsar -lgmtsar
#
#-------------------------------------------------------------------------------
#	Compiler switches and linker flags
#-------------------------------------------------------------------------------
#
CFLAGS		=  -O3 -Wall -fPIC -fno-strict-aliasing -std=gnu99 -march=native -ffast-math  -fopenmp
LDFLAGS		=  -s -Wl,--allow-multiple-definition -Wl,-rpath=/usr/lib/x86_64-linux-gnu
#
#-------------------------------------------------------------------------------
#	Shared library file extension. Examples:-march=native -ffast-math 
#	Solaris, Linux, IRIX: so	HP: sl		Mac OS X: dylib
#	SL_VERSION is extension for versioned shared libraries. Examples:
#	Linux: so.<version>		Mac OS X: <version>.dylib
#	Set SL_VERSION = $(SL) to prevent the creation of versioned shared libraries
#-------------------------------------------------------------------------------
#
SL		= so
SL_VERSION	= $(SL)
#
#-------------------------------------------------------------------------------
#	Set  LIBEXT = a      to create static libraries (default)
#	Set  LIBEXT = $(SL)  to create shared libraries (and set SL below)
#-------------------------------------------------------------------------------
#
LIBEXT		= a
#-------------------------------------------------------------------------------
#	Linker switch for building shared libraries
#	SunOS: -G  HP-UX: -b  Linux,IRIX: -shared
#-------------------------------------------------------------------------------
#
LD_OPT		= 
#-------------------------------------------------------------------------------
