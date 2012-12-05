EXTRA_CODE=Extra_code
EXTRA_HIGGS=higgs_properties
PDF_DIR=PDFTables
MG_CODE=MadGraphSrc
LIBDIR=../macros/lib
CXXFLAGS= -fPIC -O -I.
CFLAGS=-fPIC

target=$(LIBDIR)/libmodel_HEF_UFO.a

all: $(target)

objects=$(EXTRA_CODE)/MEKD_CalcHEP_PDF.o $(EXTRA_CODE)/MEKD_MG_Boosts.o $(EXTRA_HIGGS)/hggeffective.o $(MG_CODE)/HelAmps_HEF_UFO_spin0.o $(MG_CODE)/HelAmps_HEF_UFO_spin2.o $(MG_CODE)/HelAmps_HEF_UFO_bkg.o $(MG_CODE)/Parameters_HEF_UFO.o $(MG_CODE)/read_slha.o $(PDF_DIR)/pdt.o
$(target): $(objects)
	$(AR) cru $@ $(objects)
	ranlib $(target)

.PHONY: clean

clean:
	rm -f $(target)
	rm -f $(objects)

