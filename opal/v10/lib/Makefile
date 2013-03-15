#
# Makefile for the lib-dir of the 
# opacity-interpolation package
#
# G. Houdek
# 16.9.1992
#
#----------------------------------------
LIB  =  libopint.a

set_prec = ./set_single
#dble#set_prec = ./set_double
#----------------------------------------
$(LIB): $(LIB_OBJ)
	( \
         ar ruv $(LIB) $? ;\
         $(RANLIB) $(LIB) ;\
)

dopinit.f: dopinit.rf
	$(set_prec) dopinit.rf dopinit.f

dopintc.f: dopintc.rf
	$(set_prec) dopintc.rf dopintc.f

dopintf.f: dopintf.rf
	$(set_prec) dopintf.rf dopintf.f

dopints.f: dopints.rf
	$(set_prec) dopints.rf dopints.f

rat2d.f: rat2d.rf
	$(set_prec) rat2d.rf rat2d.f

rbivpd.f: rbivpd.rf
	$(set_prec) rbivpd.rf rbivpd.f
