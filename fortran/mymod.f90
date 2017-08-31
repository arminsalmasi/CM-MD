Module testmod
contains

subroutine pass_array(b)
	implicit none
	inTEGER, ALLOCATABLE :: b(:)

	allocate(b(10))

	b(1:5) = 150
	b(6:10) = 200

end subroutine pass_array
end module testmod