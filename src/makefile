#
f77=gfortran
#

NEW_OBJ=icamax.o cgefa.o scacum.o csscal.o cgesl.o\
        cgeco.o cdotc.o caxpy.o cscal.o \
        ran0.o ran1.o ran2.o ran3.o izamax.o dcabs1.o

NEW_OBJ_make=icamax cgefa scacum csscal cgesl \
        cgeco cdotc caxpy cscal \
        ran0 ran1 ran2 ran3 izamax dcabs1


all:	$(NEW_OBJ) vgffmco vgffmc vgfphz vgfin vgfkv vgffmctest  

vgffmco: vgffmc.f
	f77 -c -o vgffmc.o vgffmc.f

vgffmc:	vgffmc.o $(NEW_OBJ) 
	f77 -o vgffmc vgffmc.o $(NEW_OBJ)
 
vgffmctest:	vgffmctest.o $(NEW_OBJ) 
	f77 -o vgffmctest vgffmctest.o $(NEW_OBJ)
 
vgfin:	vgfin.f  
	f77 -o vgfin vgfin.f 
 
vgfphz:	vgfphz.f  
	f77  -o vgfphz vgfphz.f 
 
vgfkv:	vgfkv.f 
	f77  -o vgfkv vgfkv.f 
 
icamax: icamax.f
	f77 -c -o icamax.o icamax.f

cgefa: cgefa.f
	f77 -c -o cgefa.o cgefa.f

scacum: scacum.f
	f77 -c -o scacum.o scacum.f

csscal: csscal.f
	f77 -c -o csscal.o csscal.f

cgesl: cgesl.f
	f77 -c -o cgesl.o cgesl.f

cgeco: cgeco.f
	f77 -c -o cgeco.o cgeco.f

cdotc: cdotc.f
	f77 -c -o cdotc.o cdotc.f

caxpy: caxpy.f
	f77 -c -o caxpy.o caxpy.f

cscal: cscal.f
	f77 -c -o cscal.o cscal.f

ran0: ran0.f
	f77 -c -o ran0.o ran0.f

ran1: ran1.f
	f77 -c -o ran1.o ran1.f

ran2: ran2.f
	f77 -c -o ran2.o ran2.f

ran3: ran3.f
	f77 -c -o ran3.o ran3.f

izamax: izamax.f
	f77 -c -o izamax.o izamax.f

dcabs1: dcabs1.f
	f77 -c -o dcabs1.o dcabs1.f

cleano: 
	rm $(NEW_OBJ) vgffmc.o vgfphz.o vgfin.o  

clean: 
	rm $(NEW_OBJ) vgffmc.o vgfphz.o vgfin.o vgffmc vgfphz vgfin vgfkv  
