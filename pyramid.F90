!**********************************************************************
program pyramid 
!**********************************************************************
implicit none
integer, parameter :: wp = 8
integer, parameter :: ndim = 3
real(wp) :: x(ndim)
real(wp) :: base
real(wp) :: height
real(wp) :: v1
real(wp) :: v2
real(wp) :: alpha
real(wp) :: ans_pot
real(wp) :: ans_gauss



!*************************BEGIN INPUT PARAMETERS***********************



! Coordinates:

x(1) = 0.0
x(2) = 0.0
x(3) = 0.0

! Length of the base:
base = 1.0

! Height of the pyramid:
height = 10.0

! Potential inside the pyramid:
v1 = 0.0d0

! Potential outside the pyramid:
v2 = 1.0d6

! Alpha value for the gaussian:
alpha = 0.00250e0_wp


!*************************END INPUT PARAMETERS*************************


call calc_pyramid_pot(ndim,base,height,v1,v2,x,ans_pot)
call calc_gaussian(ndim,base,height,x,alpha,ans_gauss)

write(*,*) ans_pot, ans_gauss

end



!**********************************************************************
subroutine calc_pyramid_pot(ndim,base,h,v1,v2,r,ans)
!**********************************************************************

!This subroutine calculates the potential for a given point
implicit none
integer, parameter :: wp = 8
integer, intent(in) :: ndim
real(wp), intent(in) :: base
real(wp), intent(in) :: h
real(wp), intent(in) :: v1
real(wp), intent(in) :: v2
real(wp), intent(in) :: r(ndim)
real(wp), intent(out) :: ans



real(wp) :: xmin
real(wp) :: xmax
real(wp) :: ymin
real(wp) :: ymax
real(wp) :: a
real(wp) :: b
real(wp) :: x
real(wp) :: y
real(wp) :: z
real(wp) :: theta



x = r(1)
y = r(2)
z = r(3)


a = base/2.0e0_wp
b = base/2.0e0_wp



call min_max_calculation(z,a,h,xmin,xmax)
call min_max_calculation(z,b,h,ymin,ymax)



if(x .ge. xmin .and. x .le. xmax .and. y .ge. ymin &
&.and. y .le. ymax .and. z .ge. 0. .and. z .le. h) then
  theta = 1.0e0_wp
else 
  theta = 0.0e0_wp
end if


ans = (theta*v1) + ((1.0e0_wp-theta)*v2)


end subroutine calc_pyramid_pot


!**********************************************************************
subroutine calc_gaussian(ndim,base,h,r,alpha,ans)
!**********************************************************************

!This subroutine calculates the gaussian for a given point
implicit none
integer, parameter :: wp = 8
integer, intent(in) :: ndim
real(wp), intent(in) :: base
real(wp), intent(in) :: h
real(wp), intent(in) :: r(ndim)
real(wp), intent(in) :: alpha
real(wp), intent(out) :: ans

real(wp) :: xmin
real(wp) :: xmax
real(wp) :: ymin
real(wp) :: ymax
real(wp) :: a
real(wp) :: b
real(wp) :: x
real(wp) :: y
real(wp) :: z
real(wp) :: theta



x = r(1)
y = r(2)
z = r(3)

a = base/2.0e0_wp
b = base/2.0e0_wp

call min_max_calculation(z,a,h,xmin,xmax)
call min_max_calculation(z,b,h,ymin,ymax)

ans = (x-xmin)*(xmax-x)*(y-ymin)*(ymax-y)*Exp(-alpha*(x**2+y**2))

write(*,*) 'ans',ans,alpha 

end subroutine calc_gaussian



!**********************************************************************
subroutine min_max_calculation(z,a,h,xmin,xmax)
!**********************************************************************

! This subroutine calculates xmin and xmax 
implicit none
integer, parameter :: wp = 8
real(wp), intent(in) :: z
real(wp), intent(in) :: a, h
real(wp), intent(out) :: xmin, xmax

 
xmin = -a*(1.0e0_wp-(z/h))
xmax =  a*(1.0e0_wp-(z/h))


end subroutine min_max_calculation
