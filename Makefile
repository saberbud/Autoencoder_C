################################################################################
# Makefile for building: A simple neural network for pattern recognition       #
################################################################################

MAKEFILE      = Makefile

####### Compiler, tools and options
CC            = gcc
CXX           = g++
FORTRAN       = gfortran
FFLAGS        = -fPIC
#CXXFLAGS      = -shared -std=c++11 -O2 -g -pipe -Wall -m64 -fPIC
#CXXFLAGS      = -shared -std=c++11 -O2 -g -pipe -Wall -m64 -fPIC -frtti \
                   -DLINUXVERS -I/usr/local/Cellar/root/5.34.36_3/include/root -O
CXXFLAGS      = -shared -std=c++11 -O2 -g -pipe -Wall -m64 -fPIC -frtti \
                   -DLINUXVERS -I$(ROOTSYS)/include -O
LIBFLAGS      = -shared -std=c++11 -Wall -O2 -g -pipe -m64 -fPIC
#TESTFLAGS     = -std=c++11 -Wall -frtti -fno-exceptions -fPIC \
                   -DLINUXVERS -I/usr/local/Cellar/root/5.34.36_3/include/root -O
TESTFLAGS     = -std=c++11 -Wall -frtti -fno-exceptions -fPIC \
                   -DLINUXVERS -I$(ROOTSYS)/include -O
INCPATH       = -Iinclude
DEL_FILE      = rm -f
CHK_DIR_EXISTS= test -d
MKDIR         = mkdir -p
COPY          = cp -f
COPY_FILE     = cp -f
COPY_DIR      = cp -f -R
INSTALL_FILE  = install -m 644 -p
INSTALL_PROGRAM = install -m 755 -p
INSTALL_DIR   = cp -f -R
DEL_FILE      = rm -f
SYMLINK       = ln -f -s
DEL_DIR       = rmdir
MOVE          = mv -f
TAR           = tar -cf
COMPRESS      = gzip -9f
LINK          = g++
LFLAGS        = -shared
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
#ROOTLIBS      = -L/usr/local/Cellar/root/5.34.36_3/lib/root -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread -lm -ldl
#ROOTGLIBS     = -L/usr/local/Cellar/root/5.34.36_3/lib/root -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lpthread -lm -ldl
RLIBS          = $(ROOTLIBS)
GLIBS         = $(ROOTGLIBS)
LIBS          = $(GLIBS) $(RLIBS)
AR            = ar cqs
RANLIB        = 
SED           = sed
STRIP         = 

####### Targets, add new objects here
#TARGET_LIB    = libcneural.so
TARGET_test    = test
LIB_OBJ_DIR   = obj
LIB_CLASSES   = xy_utils \
				xyNeuron \
                                xyrbmunit \
				xy_W_bag \
				xy_MLP \
				xydata_mini_batch \
				xy_Autoencoder

LIB_OBJECTS   = $(addprefix $(LIB_OBJ_DIR)/, $(LIB_CLASSES:=.o))


####### Build rules
first: all

all: lib

lib: Makefile $(TARGET_test)

$(LIB_OBJ_DIR)/%.o: src/%.cpp
	$(CXX) -c $(LIBFLAGS) $(INCPATH) -o $@ $<


$(TARGET_test): $(TARGET_test).o  $(LIB_OBJECTS) $(TARGET_test).cpp
	$(LINK) $(TESTFLAGS) -o $@ $(LIB_OBJECTS) $(TARGET_test).o $(LIBS)

####### Clean
clean: cleanobj cleantest

cleanobj:
	$(DEL_FILE) $(LIB_OBJ_DIR)/*.o

cleanlib:
	$(DEL_FILE)

cleantest:
	$(DEL_FILE) $(TARGET_test)
	$(DEL_FILE) $(TARGET_test).o

####### Install

install:   FORCE

uninstall:   FORCE

FORCE:

