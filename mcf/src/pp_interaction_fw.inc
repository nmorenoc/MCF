  !----------------------------------------------------------
  ! particle_particle_interaction_fluid_wall
  !
  ! ip is fluid particle, 
  ! jp is wall boundary particle. 
  !
  ! In case of solid wall boundary,
  ! substitute jp's velocity with an artificial 
  ! velocity for no slip condition,
  ! pressure of jp with Neumann condition.
  !
  ! If non-Newtonian fluids,
  ! calculate velocity gradient tensor.
  ! pressure tensor of jp with Neumann condition.
  !
  !
  ! f_w interaction(force) takes place the same way
  ! as f_f interaction.
  !----------------------------------------------------------
  

  v_ip(1:num_dim) = this%v(1:num_dim,ip)
  v_jp(1:num_dim) = this%v(1:num_dim,jp)
  
  
  IF ( num_wall_solid > 0 ) THEN
     
     !-------------------------------------------------------
     ! wall boundary using solid particles with
     ! Morris no slip boundary condition.
     !-------------------------------------------------------
     
     IF ( wall_noslip == 2 ) THEN
        
        CALL boundary_noslip(tboundary, &
             this%x(1:num_dim,ip),this%x(1:num_dim,jp), &
             this%v(1:num_dim,ip),v_jp(1:num_dim), &
             this%id(this%sid_idx,jp), &
             stat_info_sub)
        
        IF( stat_info_sub /=0 ) THEN
           PRINT *, "pp_interaction_fw : ", &
                "Wall no slip of jp has problem !"
           stat_info = -1
           GOTO 9999
        END IF
        
     END IF ! wall_noslip
     
  END IF ! num_wall_solid > 0
  
  
  !----------------------------------------------------------
  ! For Newtonian fluid, pressure using Neumann conditon.
  !----------------------------------------------------------

  IF ( Newtonian ) THEN

     !this%rho(jp) = this%rho(ip)     
     this%p(jp)   = this%p(ip)

     
#ifdef __PARTICLES_FORCE_SEPARATE

     CALL rhs_force_ff(this%rhs, &
          this%x(1:num_dim,ip), &
          this%x(1:num_dim,jp),dij, &
          v_ip(1:num_dim),v_jp(1:num_dim), &
          this%rho(ip),this%rho(jp), &
          this%p(ip),this%p(jp), &
          this%m(ip),this%m(jp), &
          w,gradw,fip,fjp,&
          fpip,fpjp,fvip,fvjp,frip,frjp,&
          stat_info=stat_info_sub)
#else

     CALL rhs_force_ff(this%rhs, &
          this%x(1:num_dim,ip), &
          this%x(1:num_dim,jp),dij, &
          v_ip(1:num_dim),v_jp(1:num_dim), &
          this%rho(ip),this%rho(jp), &
          this%p(ip),this%p(jp), &
          this%m(ip),this%m(jp), &
          w,gradw,fip,fjp,stat_info=stat_info_sub)

#endif
     
  ELSE ! non-Newtonian
     
     !-------------------------------------------------------
     ! For non-Newtonian fluid, density and pressure 
     ! using Neumann conditon.
     !-------------------------------------------------------
     
     rho_ip = this%rho(ip)
     rho_jp = this%rho(ip)
     
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
     PRINT *, "pp_inteaction_fw : ",& 
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
           
           sip((sb-1)*num_dim+sa) = rij(sa) * fip(sb) / 2.0_MK
           !sjp((sb-1)*num_dim+sa) = -rij(sa) * fjp(sb) / 2.0_MK
           sjp((sb-1)*num_dim+sa) = sip((sb-1)*num_dim+sa)
           
        END DO
     END DO
     
#ifdef __PARTICLES_STRESS_SEPARATE
     DO sa = 1, num_dim
        DO sb = 1, num_dim
           
           spip((sb-1)*num_dim+sa) = rij(sa) * fpip(sb) / 2.0_MK
           !spjp((sb-1)*num_dim+sa) = -rij(sa) * fpjp(sb) / 2.0_MK
           spjp((sb-1)*num_dim+sa) = spip((sb-1)*num_dim+sa)
           
           svip((sb-1)*num_dim+sa) = rij(sa) * fvip(sb) / 2.0_MK
           !svjp((sb-1)*num_dim+sa) = -rij(sa) * fvjp(sb)  / 2.0_MK
           svjp((sb-1)*num_dim+sa) = svip((sb-1)*num_dim+sa)
           
           IF ( Brownian ) THEN
              
              srip((sb-1)*num_dim+sa) = rij(sa) * frip(sb) / 2.0_MK
              !srjp((sb-1)*num_dim+sa) = -rij(sa) * frjp(sb) / 2.0_MK
              srjp((sb-1)*num_dim+sa) = srip((sb-1)*num_dim+sa)
              
           END IF
           
        END DO
     END DO
#endif

  END IF
