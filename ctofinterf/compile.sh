gcc -c ffnautyex1.c ffnautyex1_bis.c ffnautyex1_tris.c ffnautyex1_quadris.c -I../nauty/nauty26r11/
gfortran -c interface.f90 sort_index.f90 fnautyex2.f90 
gfortran interface.o sort_index.o ffnautyex1.o ffnautyex1_bis.o ffnautyex1_tris.o ffnautyex1_quadris.o -o fnautyex2.x fnautyex2.o ../nauty/nauty26r11/nauty.a 
