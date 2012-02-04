  !----------------------------------------------------------
  ! Particle-particle_velocity_gradient_tensor_fluid_fluid_ip
  !
  ! Velocity gradient tensor calculation of particle ip.
  !
  ! Here we handle currently only number density 
  ! as in Espanol&Revenga Phys. Rev. E 2003 or 
  !       Hu&Adams J. Comput.Phys. 2006 or
  !       Hu&Adams Phys.Fluids 2006  formulations.
  !----------------------------------------------------------
  
  DO b = 1, num_dim  ! ---, row direction
     
     DO a =1, num_dim ! |,  column direction
        
        this%vgt(a+num_dim*(b-1),ip) = &
             this%vgt(a+num_dim*(b-1),ip) + &
             (v_jp(a) - v_ip(a)) * &
             (this%x(b,ip) - this%x(b,jp)) * &
             gradw / rho_ip / dij

     END DO
     
  END DO
