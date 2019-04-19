#PATHS
CURRENT_DIR 	= $(shell pwd)/
#CONSTANTS		= /Users/isafa/Programs/Physics_Constants/
GSL             = /usr/local/Cellar/gsl/2.5/
LHAPDF          = /usr/local/Cellar/lhapdf/6.2.1/

SOURCES		    = $(wildcard src/*.cpp)
OBJECTS 	    = $(SOURCES:.cpp=.o)

INCLUDE_PATH	=  -I/usr/local/include -I./inc
#INCLUDE_PATH   +=  -I$(CONSTANTS)inc
INCLUDE_PATH   +=  -I$(GSL)include
INCLUDE_PATH   += -I$(LHAPDF)include

#Libraries
LIBS 		   += -lgsl -lgslcblas #-lPhysicsConstants
LIB_PATH	    = -L/usr/local/lib/ -L/usr/lib
LIB_PATH	   += -L$(CONSTANTS)lib
LIB_PATH	   += -L./lib
LIB_PATH       +=  -Wl,-rpath -Wl,$(GSL)lib -L$(GSL)lib 
LIB_PATH    +=  -Wl,-rpath -Wl,$(CURRENT_DIR)/lib -L$(CURRENT_DIR)/lib/

#Compiler
CC 		        =  gcc
CXX		        =  g++
AR              =  ar

NAME		    = NewDipole
STAT_PRODUCT  	= lib$(NAME).a
DYN_PRODUCT	  	= lib$(NAME)$(DYN_SUFFIX)

#Dynamic Library
OS_NAME=$(shell uname -s)
ifeq ($(OS_NAME),Linux)
DYN_SUFFIX=.so
DYN_OPT=-shared -Wl,-soname,lib$(NAME).so
endif
ifeq ($(OS_NAME),Darwin)
DYN_SUFFIX=.dylib
DYN_OPT=-dynamiclib -install_name $(CURRENT_DIR)/lib/$(DYN_PRODUCT)
endif

#Flags
CXXFLAGS 	    =  $(INCLUDE_PATH) -O3 -fPIC -std=c++11
CXXFLAGS       += -I$(GSL)include
CFLAGS 	    =  $(INCLUDE_PATH) -O3 -fPIC -std=c++11
CFLAGS       += -I$(GSL)include

LD		        =  g++
#LDFLAGS		    =  -L/usr/local/lib/  -L/usr/lib -Wl,-rpath -Wl,$(CURRENT_DIR)/lib -L$(CURRENT_DIR)/lib/
LDFLAGS		    =  -L/usr/local/lib/  -L/usr/lib -L/usr/local/Cellar/lib
LDFLAGS		    =  -L$(LHAPDF)lib
LDFLAGS        +=  -Wl,-rpath -Wl,$(GSL)lib -L$(GSL)lib 
LDFLAGS        +=  -lgsl -lgslcblas
LDFLAGS        +=  -lLHAPDF

MAIN_SRC	    = $(wildcard mains/*.cpp)
MAIN_OBJ	    = $(MAIN_SRC:.cpp=.o)
MAIN_EXE	    = $(MAIN_SRC:.cpp=.exe)

.PHONY: all clean

all:  $(STAT_PRODUCT) $(DYN_PRODUCT) $(MAIN_EXE)

%.exe : %.cpp
	$(CXX) $(CFLAGS) $(LIB_PATH) $< -o $@ $(LIBS) -l$(NAME)
	mv $@ bin/

$(DYN_PRODUCT) : $(OBJECTS)
	@echo Linking $(DYN_PRODUCT)
	@$(LD) $(DYN_OPT)  $(LDFLAGS) -o $(DYN_PRODUCT) $(OBJECTS)
	mv $(DYN_PRODUCT) lib/$(DYN_PRODUCT)

$(STAT_PRODUCT) : $(OBJECTS)
	@echo Linking $(STAT_PRODUCT)
	@$(AR) -rcs $(STAT_PRODUCT) $(OBJECTS)
	mv $(STAT_PRODUCT) lib/$(STAT_PRODUCT)

%.o:%.cpp
	$(CXX) $(CXXFLAGS) -c $^ -o $@

clean:
	rm -f src/*.o lib/* bin/*.exe mains/*.o

