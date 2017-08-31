program array_passing_test
use testmod

INTEGER, ALLOCATABLE :: a(:)

	call pass_array(a)
	
	write(*,*) a

end program array_passing_test

