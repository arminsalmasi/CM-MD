!****************************************************************************
!
! PROGRAM: gaussianrand
!
! PURPOSE: generate gaussian random numbers with given sigma
! Uses FUNCTION RAND(0) to generate uniform random numbers
! between 0 and 1, and ROUTINE SEED to set the seed
! Uses the Box-Muller algorithm to generate numbers with
! a Gaussian distribution of probability. 
!
!****************************************************************************
program gaussianrand
implicit none

integer,parameter:: nrand=1000000
REAL*8 UNIF1(nrand), UNIF2(nrand)
!
REAL*8 Y1,Y2,SIGMA,binsize, binhalf
INTEGER I,nbins, bins(201),ix

! supply seed to get the same pseudorandom sequence each time
! program is run. To get a different sequence, change the seed
! by calling SEED with a different integer value.
!call seed(1)
! set the sigma value you want
SIGMA=2
!
! get the random numbers. RGAUSS returns two gaussian random numbers
! from the same distribution for each call. 
!
DO I=1,NRAND
CALL RGAUSS(SIGMA,Y1,Y2)
UNIF1(I)=Y1
UNIF2(I)=Y2
END DO
!
! After this loop, there are NRAND values in each array UNIF1 and UNIF2, each with
! a Gaussian distribution with standard deviation SIGMA
!
! To display the resultsfor UNIF1, bin the results usingbin size=0.05*sigma
nbins=201
binsize=0.05*sigma
binhalf=0.5*binsize

do i=1,nrand
! compute which bineach number goes into
ix=101+floor((UNIF1(i)+binhalf)/binsize)
if(ix.gt.0.and.ix.lt.201) bins(ix)=bins(ix)+1
end do

! save the numbers in each bin to display later
open(unit=10,file='rgauss.dat',form='formatted',status='unknown')
write(10,1000) (bins(i),i=1,201)
1000 format(10i8)
close(10)

end program gaussianrand

subroutine rgauss(sigma, y1,y2)
real*8 x1, x2, w, y1, y2, sigma

do while ( (w .ge. 1.0).or.(w.eq.0) )
x1 = 2.0 * rand(0) - 1.0
x2 = 2.0 * rand(0) - 1.0
w = x1 * x1 + x2 * x2
end do

w = sigma*sqrt( (-2.0 * log( w ) ) / w )
y1 = x1 * w
y2 = x2 * w

end subroutine rgauss




!program test
!
!  implicit none
!  real, dimension(:), allocatable :: tp
!   integer :: i
!
!  allocate(tp(0:1000))
!  
!  do i = 0,size(tp)
!    tp(i)= 20+i
!  end do
!
!  print *, tp
!
!end program test
