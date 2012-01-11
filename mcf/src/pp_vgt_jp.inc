  
  !----------------------------------------------------------
  ! Particle-particle_velocity_gradient_tensor_fluid_fluid_jp
  !
  ! Velocity gradient tensor calculation of particle jp.
  !  
  ! Here we handle currently only number density 
  ! as in Espanol&Revenga 2003 or Hu&Admas 2006 formulation.
  !----------------------------------------------------------

  
  DO b = 1, num_dim  ! ---, row direction
     DO a =1, num_dim ! |,  column direction
        
        this%vgt(a+num_dim*(b-1),jp) = &
             this%vgt(a+num_dim*(b-1),jp) + &
             (v_ip(a) - v_jp(a)) * &
             (this%x(b,jp) - this%x(b,ip)) * &
             gradw / rho_jp / dij
        
     END DO
  END DO




