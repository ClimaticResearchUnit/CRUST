# Uses gfortran windows/linux and dislin graphics
# TYP=WIN
TYP=LIN
# TYP=MAC

ifeq ($(TYP),WIN)
  LIB=C:\dislin\disgf_d.a \
  C:\MinGW32\lib\libcomdlg32.a \
  C:\MinGW32\lib\libgdi32.a \
  C:\MinGW32\lib\libopengl32.a 
  FO = gfortran -static -o 
endif
ifeq ($(TYP),LIN) 
  LIB=/usr/local/dislin/lib/dislin_d-10.0.a \
    /usr/lib/libXm.so.2
  FO = gfortran -o 
endif
ifeq ($(TYP),MAC) 
  LIB=/usr/local/dislin/lib/dislin_d-10.0.a \
    /usr/lib/libXm.so.4
  FO = gfortran -o 
endif
FC = gfortran -c -O2

all: CRUST
CRUST: crust.o crustprocs.o rfigs.o crustutil.o stand.o detdata.o
	$(FO) CRUST crust.o crustprocs.o rfigs.o crustutil.o stand.o detdata.o $(LIB)
ifeq ($(TYP),WIN)
	copy CRUST.exe ..
	del *.o *.mod CRUST.exe
else
	mv CRUST ../
	rm *.o *.mod  
endif

crust.o: crust.f90 crustprocs.o  
	$(FC) crust.f90 
crustprocs.o: crustprocs.f90 rfigs.o
	$(FC) crustprocs.f90  
rfigs.o: rfigs.f90 crustutil.o 
	$(FC) rfigs.f90 
crustutil.o: crustutil.f90 stand.o 
	$(FC) crustutil.f90 
stand.o: stand.f90 detdata.o dislin.o
	$(FC) stand.f90 
detdata.o: detdata.f90 
	$(FC) detdata.f90
dislin.o: dislin.f90
	$(FC) dislin.f90
# Cleaning the project
clean:
ifeq ($(TYP),WIN)
	del ..\CRUST.exe
else
	rm ../CRUST   
endif
# End of the makefile
