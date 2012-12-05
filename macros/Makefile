LIBDIR=lib
INCDIR1=../src/Extra_code
INCDIR2=../src/MadGraphSrc
INCDIR3=../src
INCDIR4=../interface
PROC_DIR=../src/MadGraphSrc
MODELLIB=model_HEF_UFO
CXXFLAGS= -fPIC -O -I$(INCDIR1) -I$(INCDIR2) -I$(INCDIR3) -I$(INCDIR4) -I.
CFLAGS=-fPIC
LIBFLAGS= -L$(LIBDIR) -l$(MODELLIB)
UNAME := $(shell uname)

########## ROOT PART ##########
ROOTLIB=$(ROOTSYS)/lib
ROOTINCLUDES=$(ROOTSYS)/include
###############################

objects=$(PROC_DIR)/BKG_DN_OF.o $(PROC_DIR)/BKG_UP_OF.o $(PROC_DIR)/BKG_DN_SF.o $(PROC_DIR)/BKG_UP_SF.o $(PROC_DIR)/Spin0_OF.o $(PROC_DIR)/Spin0_SF.o $(PROC_DIR)/Spin2_OF.o $(PROC_DIR)/Spin2_SF.o 

#../src/MEKD_MG.o

targetlib=$(LIBDIR)/libMEKD_MG.so
target=runKD_MAD


all: $(target)


ifeq ($(UNAME), Linux)
$(targetlib): $(objects) $(LIBDIR)/lib$(MODELLIB).a
	$(CXX) -fPIC -shared -Wl,-soname,libMEKD_MG.so -o $(targetlib) $(objects) $(LIBFLAGS)

ifeq ($(MEKD_COMPILE_WITH_ROOT), Yes)
$(target): $(targetlib)
	$(CXX) -D MEKD_with_ROOT `$(ROOTSYS)/bin/root-config --ldflags` -I $(ROOTINCLUDES) $(target).cpp -o $(target) -L. $(LIBFLAGS) -L $(ROOTLIB) `$(ROOTSYS)/bin/root-config --libs --glibs` $(targetlib)
else
$(target): $(targetlib)
	$(CXX) $(target).cpp -o $(target) -L. $(LIBFLAGS) $(targetlib)
endif

endif

ifeq ($(UNAME), Darwin)
$(targetlib): $(objects) $(LIBDIR)/lib$(MODELLIB).a
	$(CXX) -fPIC -shared -Wl,-install_name,libMEKD_MG.so -o $(targetlib) $(objects) $(LIBFLAGS)

ifeq ($(MEKD_COMPILE_WITH_ROOT), Yes)
$(target): $(targetlib)
	$(CXX) -D MEKD_with_ROOT `$(ROOTSYS)/bin/root-config --ldflags` -I$(INCDIR3) -I$(INCDIR4) -I $(ROOTINCLUDES) $(target).cpp -o $(target) -L. $(LIBFLAGS) -L $(ROOTLIB) `$(ROOTSYS)/bin/root-config --libs --glibs` $(targetlib)
else
$(target): $(targetlib)
	$(CXX) $(target).cpp -o $(target) -L. $(LIBFLAGS) $(targetlib)
endif

endif



$(LIBDIR)/lib$(MODELLIB).a:
	cd ../src ; \
	make

.PHONY: clean

clean:
	rm -f $(targetlib); \
	rm -f $(target); \
	rm -f $(objects); \
	cd ../src; \
	make clean
