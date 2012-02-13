  !----------------------------------------------------------
  ! particle-particle_interaction_colloid_colloid
  !
  ! ip is colloid boundary particle, 
  ! jp is colloid boundary particle. 
  !
  ! Note that ip and jp are not from the same colloid.
  !
  ! ip and jp's pressure are from equation of state.
  !
  ! If Non-Newtonian fluids, 
  ! not implemented.
  !
  !
  ! c_c interaction takes place the same way
  ! as f_f interaction.
  !----------------------------------------------------------

  
  v_ip(1:num_dim) = this%v(1:num_dim,ip)
  v_jp(1:num_dim) = this%v(1:num_dim,jp)
  
  
  IF ( Newtonian ) THEN

     CALL rhs_force_ff(this%rhs,&
          this%x(1:num_dim,ip),&
          this%x(1:num_dim,jp),dij,&
          v_ip(1:num_dim),&
          v_jp(1:num_dim),&
          this%rho(ip),this%rho(jp),&
          this%p(ip),this%p(jp),&
          this%m(ip),this%m(jp),&
          w,gradw,fip,fjp,stat_info=stat_info_sub)

#if 0     
     ! trying lubrication on individual boundary particles.
     rij(1:num_dim) = rij(1:num_dim) / dij
     uij(1:num_dim) = v_ip(1:num_dim) - v_jp(1:num_dim)
     h = dij - 2.0_MK* dx(1)
     IF ( h < dx(1)/1000.0 ) THEN
        h = dx(1)/1000.0_MK
     END IF
     
     fip(1:num_dim) = -3.0_MK * mcf_pi * eta * &
          DOT_PRODUCT(uij(1:num_dim), rij(1:num_dim)) * &
          dx(1)**2/2.0_MK/h * rij(1:num_dim)
     fjp(1:num_dim) = -fip(1:num_dim)
     
#endif

  ELSE ! non-Newtonian
     
     rho_ip = this%rho(ip)
     rho_jp = this%rho(jp)
     
#include "pp_vgt_ip.inc"
     
     IF ( symmetry ) THEN
        
#include "pp_vgt_jp.inc"
        
     END IF
     
     this%pt(1:num_dim,1:num_dim,jp) = &
          this%pt(1:num_dim,1:num_dim,ip)
     
     CALL rhs_force_ff(this%rhs, &
          this%x(1:num_dim,ip), &
          this%x(1:num_dim,jp),dij, &
          v_ip(1:num_dim), &
          v_jp(1:num_dim), &
          rho_ip,rho_jp, &
          this%pt(1:num_dim,1:num_dim,ip), &
          this%pt(1:num_dim,1:num_dim,jp), &
          this%m(ip),this%m(jp), &
          w,gradw,fip,fjp,stat_info=stat_info_sub)
     
  END IF ! Newtonian
  
  IF( stat_info_sub /=0 ) THEN
     PRINT *, "pp_inteaction_cc : ",& 
          "Pairwise force has problem !"
     stat_info = -1
     GOTO 9999
  END IF

 !----------------------------------------------------------
  ! Calculate stress tensor on each particle.
  ! For mi=mj, rij=-rji,
  ! sip(:) = sjp(:)
  !----------------------------------------------------------
  

  IF ( stress_tensor ) THEN

     DO sa = 1, num_dim
        DO sb = 1, num_dim
           
           sip((sb-1)*num_dim+sa) = rij(sa) * fip(sb)
           !sjp((sb-1)*num_dim+sa) = -rij(sa) * fjp(sb)
           sjp((sb-1)*num_dim+sa) = sip((sb-1)*num_dim+sa)
           
        END DO
     END DO
     
#ifdef __PARTICLES_STRESS_SEPARATE
     DO sa = 1, num_dim
        DO sb = 1, num_dim
           
           spip((sb-1)*num_dim+sa) = rij(sa) * fpip(sb)
           !spjp((sb-1)*num_dim+sa) = -rij(sa) * fpjp(sb)
           spjp((sb-1)*num_dim+sa) = spip((sb-1)*num_dim+sa)
           
           svip((sb-1)*num_dim+sa) = rij(sa) * fvip(sb)
           !svjp((sb-1)*num_dim+sa) = -rij(sa) * fvjp(sb)
           svjp((sb-1)*num_dim+sa) = svip((sb-1)*num_dim+sa)
           
           IF ( Brownian ) THEN
              
              srip((sb-1)*num_dim+sa) = rij(sa) * frip(sb)
              !srjp((sb-1)*num_dim+sa) = -rij(sa) * frjp(sb)
              srjp((sb-1)*num_dim+sa) = srip((sb-1)*num_dim+sa)
              
           END IF
           
        END DO
     END DO
#endif

  END IF
