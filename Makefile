# Build libmonoHZZ4L
# Created 13 Jun 2015 HBP - Les Houches
# Updated 01 Jan 2017 HBP - for RGS paper
# ----------------------------------------------------------------------------
ifndef ROOTSYS
$(error * Define ROOTSYS)
endif

ifndef DELPHES
$(error * define DELPHES. Edit setup.(c)sh then source setup.(c)sh)
endif
# ----------------------------------------------------------------------------
NAME	:= monoHZZ4L
incdir	:= include
srcdir	:= src
libdir	:= lib
bindir	:= bin

$(shell mkdir -p lib)

# get lists of sources
APPS	:=	analyzeHZZ4L
PROGRAMS:=	$(foreach x,$(APPS),$(bindir)/$x)
APPSRCS	:=	$(foreach x,$(APPS),$(srcdir)/$x.cc)

SRCS	:= 	$(srcdir)/monoHZZ4L.cc \
		$(srcdir)/nic.cc \
		$(srcdir)/LeptonEfficiency.cc

CINTSRCS:= $(wildcard $(srcdir)/*_dict.cc)

OTHERSRCS:= $(filter-out $(CINTSRCS) $(APPSRCS) $(SRCS),$(wildcard $(srcdir)/*.cc))

# list of dictionaries to be created
DICTIONARIES	:= $(SRCS:.cc=_dict.cc)

# get list of objects
OBJECTS		:= $(OTHERSRCS:.cc=.o) $(SRCS:.cc=.o) $(DICTIONARIES:.cc=.o)

APPOBJS	:= $(APPSRCS:.cc=.o)

ALLOBJECTS	:= $(OBJECTS) $(APPOBJS)



#say := $(shell echo "DICTIONARIES:     $(DICTIONARIES)" >& 2)
#say := $(shell echo "" >& 2)
#say := $(shell echo "APPS: $(APPS)" >& 2)
#say := $(shell echo "APPOBJS: $(APPOBJS)" >& 2)
#say := $(shell echo "ALLOBJECTS: $(ALLOBJECTS)" >& 2)
#$(error bye)
# ----------------------------------------------------------------------------
ROOTCINT	:= rootcint

# check for clang++, otherwise use g++
COMPILER	:= $(shell which clang++ >& $(HOME)/.cxx; tail $(HOME)/.cxx)
COMPILER	:= $(shell basename "$(COMPILER)")
ifeq ($(COMPILER),clang++)
CXX		:= clang++
LD		:= clang++
else
CXX		:= g++
LD		:= g++
endif
CPPFLAGS	:= -I. -I$(DELPHES) -I$(DELPHES)/external -I$(incdir)
CXXFLAGS	:= -O2 -Wall -fPIC -ansi -Wshadow -Wextra \
$(shell root-config --cflags)
LDFLAGS		:= -g
# ----------------------------------------------------------------------------
# which operating system?
OS := $(shell uname -s)
ifeq ($(OS),Darwin)
	LDFLAGS += -dynamiclib
	LDEXT	:= .dylib
else
	LDFLAGS	+= -shared
	LDEXT	:= .so
endif	
LDFLAGS += $(shell root-config --ldflags) -L$(DELPHES)
LIBS 	:= -lDelphes -lPyROOT $(shell root-config --libs --nonew)
LIBRARY	:= $(libdir)/lib$(NAME)$(LDEXT)
# ----------------------------------------------------------------------------
all: $(PROGRAMS)

lib: $(LIBRARY)

$(PROGRAMS)	: 	$(bindir)/%	:	$(srcdir)/%.o 	$(LIBRARY)
	@echo ""
	@echo "=> Linking program $@"
	$(LD) $(LDFLAGS) $^ $(LIBS) -L$(libdir) -l$(NAME) -o $@

$(LIBRARY)	: $(OBJECTS)
	@echo ""
	@echo "=> Linking shared library $@"
	$(LD) $(LDFLAGS) $^ $(LIBS)  -o $@

$(ALLOBJECTS)	: %.o	: 	%.cc
	@echo ""
	@echo "=> Compiling $<"
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

$(DICTIONARIES)	: $(srcdir)/%_dict.cc	: $(incdir)/%.h
	@echo ""
	@echo "=> Building dictionary $@"
	$(ROOTCINT) -f $@ -c -I$(DELPHES) $^
	find $(srcdir) -name "*.pcm" -exec mv {} $(libdir) \;

tidy:
	rm -rf $(srcdir)/*_dict*.* $(srcdir)/*.o 

clean:
	rm -rf $(libdir)/* $(srcdir)/*_dict*.* $(srcdir)/*.o $(PROGRAMS)
