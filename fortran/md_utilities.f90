module utilities
  
  use datastructure
  
  contains  
!##########################################################
    subroutine do_allocate_dtyp()

    ! allocate tmstp datatype vector
    ! allocate all vector fileds of all cells of tmstp vector
    
      allocate(tmstp(N_tstps))
      
    ! loops over all timesteps to allocate variables in 
    ! the tmstp elements(datatpe)
      do i = 1, N_tstps 
        allocate(tmstp(i)%xyz(N_atms,3)) ! N*3 matrix
        allocate(tmstp(i)%vel(N_atms,3)) ! N*3 matrix
        allocate(tmstp(i)%acc(N_atms,3)) ! N*3 matrix
      end do
      
   end subroutine do_allocate_dtyp
   
!##########################################################
   
   subroutine do_loop_tstps()

   ! loops over all timesteps(N_tstp=t/dt)
   !   in tstp=0:
   !   in 0<tstp<N_tstp

   integer :: j,seed
     
     do i= 0 , N_tstps
       if (i == 0) then
         ! TODO : intialize tstp 0 
         ! call do_rand_xyz()
         seed = 86456
         do j = 1, N_atms
           call srand(seed) 
           print* , j, rand(seed)*((vol_box)**(1/3))
         end do 
         ! call do_rand_vels()
       else
         ! TODO : do md on each timestep
         ! do_md()
         ! save_tstp() 
       end if
     end do
     
   end subroutine do_loop_tstps
  
  !##########################################################

end module utilities
