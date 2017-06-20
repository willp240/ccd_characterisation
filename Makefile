LIBNAME  = DmtpcSkim

CXX=g++
LDFLAGS=`root-config --libs` `dmtpc-config --lib-core`
CXXFLAGS= `root-config --cflags` -fPIC -O3 -g   `dmtpc-config --inc-core`
#LDFLAGS+= -Wl,-z,defs 
LDFLAGS+= -rdynamic

# add JanEvent location
#CXXFLAGS+= -I../DmtpcJanEve/include 
#LDFLAGS+=-L/home/balewski/dmtpc-software//DmtpcJanEve/lib -lM3Event

VPATH=src:include:build
LIBDIR=lib
BUILDDIR=build
INCLUDEDIR=include
BINDIR=bin

#SRCS_J     = WeventIO.cc CcdPed_init.cc CcdPedMaker.cc CcdCalibMaker.cc CcdCluster_init.cc CcdClusterMaker.cc  AnaM3EveBlue.cc
SRCS_J     = CcdPed_init.cc CcdPedMaker.cc 
OBJS_J   = $(notdir $(patsubst %.cc,%.o,$(SRCS_J)))

OBJS := $(addprefix $(BUILDDIR)/, $(LIBNAME)Cint.o  $(OBJS_J) )
INCLUDES := $(addprefix $(INCLUDEDIR)/, $(shell ls $(INCLUDEDIR)))


all: shared

shared: $(LIBDIR)/lib$(LIBNAME).so

# UGLY solution - why libDmtpcCore.so must be loaded here  and not in .C

$(LIBDIR)/lib$(LIBNAME).so: $(OBJS) | $(LIBDIR)
	@echo Building shared library $@
	$(CXX) $(LDFLAGS) $(OBJS)  -g  -shared ../DmtpcCore/lib/libDmtpcCore.so -o $@

$(OBJS): | $(BUILDDIR)

$(LIBDIR): 
	mkdir -p $(LIBDIR)

$(BUILDDIR): 
	mkdir -p $(BUILDDIR)

$(BUILDDIR)/$(LIBNAME)Cint.cc: $(INCLUDES) LinkDef.h | $(BUILDDIR)
	@echo Running rootcint
	rootcint -f $(BUILDDIR)/$(LIBNAME)Cint.cc -c $(INCLUDES) LinkDef.h

$(BUILDDIR)/%.o: %.cc $(INCLUDES) Makefile | $(BUILDDIR) 
	@echo Compiling  $< 
	$(CXX) $(CXXFLAGS) -I./include -I./ -o $@ -c $< 

clean:
	rm -rf build
	rm -rf lib



