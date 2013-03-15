FFLAGS = -O2
FC = f77
LDIR = /ai43/jcd/47tuc/evolpack_071206/gensr

LOFILES = \
$(LDIR)/blstio.o $(LDIR)/ofiles.o $(LDIR)/length.o $(LDIR)/skpcom.o \
$(LDIR)/vinta.d.o $(LDIR)/derive.d.o $(LDIR)/derivk.new.d.o \
$(LDIR)/leq.d.o $(LDIR)/vintk.d.o $(LDIR)/nrk.d.o \
$(LDIR)/zero.d.o $(LDIR)/store.d.o $(LDIR)/lir.d.o 

OFILES1 = \
kerexact.o calci.o rdaa5.d.o stop1.o

default: kerexact.x

kerexact.x: $(OFILES1) lofiles
	$(FC) $(FFLAGS) -o kerexact.x $(OFILES1) $(LOFILES)

ofiles: $(OFILES)

lofiles:
	cd $(LDIR); make $(LOFILES)
