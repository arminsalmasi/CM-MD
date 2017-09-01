program array_passing_test
!use testmod

INTEGER :: a(5,5) , b(1:5*5)
integer ix, jx, kx

	!call pass_array(a)
	
	!write(*,*) a
	
	a(1,:) = (/1,2,3,4,5/)
	a(2,:) = (/6,7,8,9,10/)
	a(3,:)=(/11,12,13,14,15/)
	a(4,:)=(/16,17,18,19,20/)
	a(5,:)=(/21,22,23,24,25/)
	
	
    kx=1
	do ix = 1 , 5
	  !do jx = 1 , 5
		b(:) = a(ix,:)
		kx=kx+1
	  !end do
	end do
!	print*, a
	print*, b

end program array_passing_test

