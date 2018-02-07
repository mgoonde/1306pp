gcc -c ffnautyex1.c ffnautyex1_bis.c ffnautyex1_tris.c -I../../nauty/nauty26r11/
gfortran -c interface.f90 fnautyex1.f90
gfortran interface.o ffnautyex1.o ffnautyex1_bis.o ffnautyex1_tris.o -o fnautyex1.x fnautyex1.o ../../nauty/nauty26r11/nauty.a 
