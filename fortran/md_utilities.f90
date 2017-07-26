module utilities
  
  use datastructure
  
  contains  
!##########################################################
    
    subroutine do_allocate_dtyp(tmstp)
    
    ! allocate tmstp datatype vector
    ! allocate all vector fileds of all cells of tmstp vector
    
      type(tmstp_holder), dimension(:), allocatable :: tmstp
      allocate(tmstp(0:floor(t/dt)))
     
    ! loops over all timesteps to allocate variables in 
    ! the tmstp elements(datatpe)
      do i = 0,N_tstps 
        allocate(tmstp(i)%xyz(N_atms,3)) ! N*3 matrix
        allocate(tmstp(i)%vel(N_atms,3)) ! N*3 matrix
        allocate(tmstp(i)%acc(N_atms,3)) ! N*3 matrix
      end do
      
    end subroutine do_allocate_dtyp
   
   !##########################################################
   
   subroutine do_loop_tstps(tmstp)
   
  
 
   ! loops over all timesteps(N_tstp=t/dt)
   !   in tstp=0:
   !   in 0<tstp<N_tstp
 

     type(tmstp_holder), dimension(:), allocatable :: tmstp
     integer :: i ,j 
     real :: r(3)
 
     do i= 0 , N_tstps
       if (i == 0) then
         ! TODO : intialize tstp 0 
         ! call do_rand_xyz()
         do j = 1, N_atms
           
           ! randomize xyz of atoms in cell 0 
           call do_genRand(r,3,j)
           tmstp(i)%xyz(j,:) = r(:) * ( box_vol**(1.0/3.0))
           
           ! randomize velocities of atoms in cell zero
           !vstd = sqrt( tmp * kb / atm_masses(j) )
           !call do_genRand(r,3,j)
           !tmstp(i)%vel(j,:) = 0 + r * vstd


         end do 
       else
         ! TODO : do md on each timestep
         ! do_md()
!         print *, tmstp(i)%xyz(j,:)            
         ! save_tstp() 
       end if
     end do
     
   end subroutine do_loop_tstps
   
   !##########################################################
   
   subroutine do_genRand(r,m,j)
   
   !
   !
   !
     integer :: i, j, n, clock
     integer, allocatable :: seed(:)
     real :: r(m)
     
     call random_seed(size=n)
     allocate(seed(n))
     call system_clock(count=clock)
     clock = clock * j 
     seed = clock + 37 * (/ (i - 1, i = 1, n) /)
     call random_seed(put=seed)
     call random_number(r)
     deallocate(seed)
   end subroutine do_genRand
   
  !########################################################## 


end module utilities








































































