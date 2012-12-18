LIBDIR=lib
INCDIR1=../src
INCDIR2=../src/Extra_code
INCDIR3=../src/higgs_properties
INCDIR4=../src/MadGraphSrc
INCDIR5=../interface
PDF_DIR=../src/PDFTables
CXXFLAGS= -pedantic -fPIC -O -I$(INCDIR1) -I$(INCDIR2) -I$(INCDIR3) -I$(INCDIR4) -I$(INCDIR5) -I.
CFLAGS=-pedantic -fPIC
LIBFLAGS= -L$(LIBDIR)
UNAME := $(shell uname)

########## ROOT PART ##########
ROOTLIB=$(ROOTSYS)/lib
ROOTINCLUDES=$(ROOTSYS)/include
###############################

objects1=$(INCDIR4)/BKG_DN_OF.o $(INCDIR4)/BKG_UP_OF.o $(INCDIR4)/BKG_DN_SF.o $(INCDIR4)/BKG_UP_SF.o $(INCDIR4)/Spin0_OF.o $(INCDIR4)/Spin0_SF.o $(INCDIR4)/Spin1_OF.o $(INCDIR4)/Spin1_SF.o $(INCDIR4)/Spin2_OF.o $(INCDIR4)/Spin2_SF.o $(INCDIR4)/HelAmps_HEF_UFO_spin0.o $(INCDIR4)/HelAmps_HZZ_Unitary_spin1.o $(INCDIR4)/HelAmps_HEF_UFO_spin2.o $(INCDIR4)/HelAmps_HEF_UFO_bkg.o $(INCDIR4)/Parameters_HEF_UFO.o $(INCDIR4)/read_slha.o

objects2=$(INCDIR2)/MEKD_CalcHEP_PDF.o $(INCDIR2)/MEKD_CalcHEP_Extra_functions.o $(INCDIR2)/MEKD_MG_Boosts.o $(INCDIR3)/hggeffective.o $(PDF_DIR)/pdt.o

objects3=$(INCDIR1)/MEKD_MG.o

objects=$(objects1) $(objects2) $(objects3)


targetlib=$(LIBDIR)/libMEKD_MG.so
target=runKD_MAD



all: $(target)




ifeq ($(UNAME), Linux)
# do something Linux-y
ifeq ($(MEKD_STANDALONE), Yes)
CXXFLAGS+= -D MEKD_STANDALONE
endif
$(targetlib): $(objects)
	$(CXX) $(CXXFLAGS) -fPIC -shared -Wl,-soname,libMEKD_MG.so -o $(targetlib) $(objects) $(LIBFLAGS)


ifeq ($(MEKD_COMPILE_WITH_ROOT), Yes)
CXXFLAGS+= -D MEKD_with_ROOT -std=c++0x
$(target): $(targetlib)
	$(CXX) $(CXXFLAGS) `$(ROOTSYS)/bin/root-config --ldflags` -I $(ROOTINCLUDES) $(target).cpp -o $(target) -L. $(LIBFLAGS) -L $(ROOTLIB) `$(ROOTSYS)/bin/root-config --libs --glibs` $(targetlib)
else
$(target): $(targetlib)
	$(CXX) $(CXXFLAGS) $(target).cpp -o $(target) -L. $(LIBFLAGS) $(targetlib)
endif

endif



ifeq ($(UNAME), $(filter $(UNAME),OSX Darwin))
# do something OSX-y
ifeq ($(MEKD_STANDALONE), Yes)
CXXFLAGS+= -D MEKD_STANDALONE
endif
$(targetlib): $(objects)
	$(CXX) -fPIC -shared -Wl,-install_name,libMEKD_MG.so -o $(targetlib) $(objects) $(LIBFLAGS)


ifeq ($(MEKD_COMPILE_WITH_ROOT), Yes)
CXXFLAGS+= -D MEKD_with_ROOT -std=c++0x
$(target): $(targetlib)
	$(CXX) $(CXXFLAGS) `$(ROOTSYS)/bin/root-config --ldflags` -I $(ROOTINCLUDES) $(target).cpp -o $(target) -L. $(LIBFLAGS) -L $(ROOTLIB) `$(ROOTSYS)/bin/root-config --libs --glibs` $(targetlib)
else
$(target): $(targetlib)
	$(CXX) $(CXXFLAGS) $(target).cpp -o $(target) -L. $(LIBFLAGS) $(targetlib)
endif

endif



.PHONY: clean

clean:
	rm -f $(targetlib); \
	rm -f $(target); \
	rm -f $(objects)
