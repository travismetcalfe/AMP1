FFLAGS = -O2
FC = f77
JCD = /local/i/jcd/evolpack_10/evolpack_ref
TJCD = $(JCD)
EDIR = $(TJCD)/prg
MDIR = $(EDIR)/mhd
LDIR = $(EDIR)/livermore
GDIR = $(JCD)/gensr
ADI = $(JCD)/adi
ADISR = $(JCD)/auxprg/adisr

GFILES = \
$(GDIR)/blstio.o $(GDIR)/length.o $(GDIR)/lir.d.o $(GDIR)/ofiles.o \
$(GDIR)/ofiles.c.o \
$(GDIR)/rnmean.d.o $(GDIR)/skpcom.o $(GDIR)/spline.d.o $(GDIR)/vinta.d.o

ASR = intmdl.d.o \
	rdamdl.n.d.o \
	rdfrqe.d.o \
	setssm.d.o \
	srdist.c.d.o \
	rseta4.c.d.o \
	wrfreq.d.o \
	rdfreq.d.o \
	rdfrqm.d.o 

FORM-GONG.N.D.O = form-gong.n.d.o $(EDIR)/setcns-g.d.o $(GDIR)/zero.d.o \
	$(GDIR)/skpcom.o $(EDIR)/wrtbc.d.o $(EDIR)/wdtgng.d.o \
	$(EDIR)/wdtmod.d.o $(EDIR)/rdgong.n.d.o

DIFF-GONG.N.D.O = diff-gong.n.d.o $(GDIR)/lir.d.o $(EDIR)/rdgong.n.d.o 

DIFF-FGONG.D.O = diff-fgong.d.o $(GDIR)/lir.d.o

all: gfiles bscan-emdl.n.d.x scan-emdl.n.d.x sel-emdl.n.d.x \
	form-gong.n.d.x diff-fgong.d.x diff-gong.n.d.x \
	redistrb.d.x redistrb.n.d.x redistrb.c.d.x scan-agsm.d.x scan-amdl.d.x \
	scan-amde.d.x set-obs.d.x

bscan-emdl.n.d.x: bscan-emdl.n.d.o
	$(FC) $(FFLAGS) -o bscan-emdl.n.d.x bscan-emdl.n.d.o \
	$(EDIR)/rdemdl.n.d.o $(GDIR)/blstio.o

scan-emdl.n.d.x: scan-emdl.n.d.o
	$(FC) $(FFLAGS) -o scan-emdl.n.d.x scan-emdl.n.d.o \
	$(EDIR)/rdemdl.n.d.o $(GDIR)/blstio.o

sel-emdl.n.d.x: sel-emdl.n.d.o
	$(FC) $(FFLAGS) -o sel-emdl.n.d.x sel-emdl.n.d.o \
	$(EDIR)/rdemdl.n.d.o $(EDIR)/wremdl.n.d.o  \
	$(EDIR)/wdtmod.d.o $(GDIR)/blstio.o

form-gong.n.d.x: $(FORM-GONG.N.D.O)
	$(FC) $(FFLAGS) -o form-gong.n.d.x $(FORM-GONG.N.D.O)

diff-fgong.d.x: $(DIFF-FGONG.D.O)
	$(FC) $(FFLAGS) -o diff-fgong.d.x $(DIFF-FGONG.D.O)

diff-gong.n.d.x: $(DIFF-GONG.N.D.O)
	$(FC) $(FFLAGS) -o diff-gong.n.d.x $(DIFF-GONG.N.D.O)

redistrb.d.x: redistrb.d.o asr
	$(FC) $(FFLAGS) -o redistrb.d.x redistrb.d.o \
	$(ADISR)/rdamdl.n.d.o \
	$(GDIR)/ofiles.o $(GDIR)/length.o \
	$(GDIR)/lir.d.o $(GDIR)/vinta.d.o \
	$(GDIR)/rnmean.d.o $(GDIR)/spline.d.o $(GDIR)/blstio.o

redistrb.n.d.x: redistrb.n.d.o
	$(FC) $(FFLAGS) -o redistrb.n.d.x redistrb.n.d.o \
	$(GDIR)/ofiles.o $(GDIR)/length.o \
	$(GDIR)/lir.d.o $(GDIR)/vinta.d.o \
	$(GDIR)/rnmean.d.o $(GDIR)/blstio.o

redistrb.c.d.x: redistrb.c.d.o asr
	$(FC) $(FFLAGS) -o redistrb.c.d.x redistrb.c.d.o \
	$(ADISR)/srdist.c.d.o $(ADISR)/rseta4.c.d.o \
	$(GDIR)/ofiles.c.o $(GDIR)/length.o \
	$(GDIR)/lir.d.o $(GDIR)/vinta.d.o \
	$(GDIR)/rnmean.d.o $(GDIR)/blstio.o

scan-agsm.d.x: scan-agsm.d.o
	$(FC) $(FFLAGS) -o scan-agsm.d.x scan-agsm.d.o $(GDIR)/blstio.o

scan-amdl.d.x: scan-amdl.d.o
	$(FC) $(FFLAGS) -o scan-amdl.d.x scan-amdl.d.o $(GDIR)/blstio.o

scan-amde.d.x: scan-amde.d.o
	$(FC) $(FFLAGS) -o scan-amde.d.x scan-amde.d.o $(GDIR)/blstio.o

set-obs.d.x: set-obs.d.o asr
	$(FC) $(FFLAGS) -o set-obs.d.x set-obs.d.o \
	$(ADISR)/rdfreq.d.o $(GDIR)/skpcom.o $(GDIR)/blstio.o

asr:
	cd $(ADISR); make $(ASR)

gfiles:
	cd $(GDIR); make $(GFILES)

#  Files depending on include file
