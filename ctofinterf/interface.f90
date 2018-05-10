module f90nautyinterf

interface

subroutine c_ffnautyex1(n) BIND(C, name="ffnautyex1")
use ISO_C_BINDING
integer(kind=C_INT) :: n
end subroutine

subroutine c_ffnautyex1_bis(n,matrix) BIND(C, name="ffnautyex1_bis")
use ISO_C_BINDING
integer(kind=C_INT) :: n
integer(kind=C_INT), dimension(n,n) :: matrix
end subroutine

subroutine c_ffnautyex1_tris(n,matrix,vector1,vector2) BIND(C, name="ffnautyex1_tris")
use ISO_C_BINDING
integer(kind=C_INT) :: n
integer(kind=C_INT), dimension(n,n) :: matrix
integer(kind=C_INT), dimension(n) :: vector1, vector2
end subroutine

subroutine c_ffnautyex1_quadris(n,matrix,vector1,vector2) BIND(C, name="ffnautyex1_quadris")
use ISO_C_BINDING
integer(kind=C_INT) :: n
integer(kind=C_INT), dimension(n,n) :: matrix
integer(kind=C_INT), dimension(n) :: vector1, vector2
end subroutine

subroutine c_ffnautyex1_petka(n,matrix,vector1,vector2,hv1,hv2,hv3) BIND(C, name="ffnautyex1_petka")
use ISO_C_BINDING
integer(kind=C_INT) :: n
integer(kind=C_INT) :: hv1,hv2,hv3
integer(kind=C_INT), dimension(n,n) :: matrix
integer(kind=C_INT), dimension(n) :: vector1, vector2
end subroutine

subroutine c_ffnautyex1_sestic(n,matrix,vector1,vector2,vector3,hv1,hv2,hv3) BIND(C, name="ffnautyex1_sestic")
use ISO_C_BINDING
integer(kind=C_INT) :: n
integer(kind=C_INT) :: hv1,hv2,hv3
integer(kind=C_INT), dimension(n,n) :: matrix
integer(kind=C_INT), dimension(n) :: vector1, vector2, vector3
end subroutine

end interface

end module
