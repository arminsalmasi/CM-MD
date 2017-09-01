program array_passing_test
!use testmod

real :: a(2,3) , b(2,3), res1(2,3), r(2), sumr
integer ix, jx, kx

	!call pass_array(a)
	
	!write(*,*) a
	
	a(1,:) = (/1,2,3/)
	
	a(2,:) = (/3,4,5/)
	b(1,:) = (/6,7,8/)
	b(2,:) = (/9,10,11/)
	!a(3,:)=(/11,12,13,14,15/)
	!a(4,:)=(/16,17,18,19,20/)
	!a(5,:)=(/21,22,23,24,25/)
	do ix = 1,2
		!res1(ix,:) = sqrt( sum( (a(ix,:) - b(ix,:) )**2))
	    !res1(ix,:) = res1(ix,:) **2 
		sumr = sqrt(sum( (a(ix,:) - b(ix,:) )**2))
		    print*, sumr
		!r(ix) = sqrt(sumr)
		!print*, r(ix)
	end do
    
	!print*, res1(1,:)
	!print*, res1(2,:)

	 

end program array_passing_test

