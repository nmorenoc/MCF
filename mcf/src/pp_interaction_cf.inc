  !----------------------------------------------------------
  ! particle-particle_interaction_colloid_fluid
  !
  ! ip is colloid boundary particle. 
  ! jp is fluid particle.
  !
  ! Substitute ip's velocity with an artificial 
  ! velocity for no slip condition,
  ! pressure jp with Neumann condition.
  !
  ! If Non-Newtonian fluids, 
  ! calculat velocity gradient tensor.
  ! pressure tensor jp with Neumann condition.
  !
  ! c_f interaction takes place the same way
  ! as f_f interaction.
  !----------------------------------------------------------
  

  v_ip(1:num_dim) = this%v(1:num_dim,ip)
  v_jp(1:num_dim) = this%v(1:num_dim,jp)
  

  !----------------------------------------------------------
  ! Morris no slip conditin for colloid.
  !----------------------------------------------------------

  IF ( coll_noslip == 2 ) THEN

     CALL colloid_noslip(colloids,&
          this%x(1:num_dim,jp),this%x(1:num_dim,ip),&
          this%v(1:num_dim,jp),v_ip(1:num_dim),&
          this%id(this%sid_idx,ip),&
          stat_info_sub)
     
     IF( stat_info_sub /=0 ) THEN
        PRINT *, "pp_interaction_cf : ",& 
             "Colloid no slip of ip has problem !"
        stat_info = -1
        GOTO 9999
     END IF
     
  END IF  ! coll_noslip
  
  !----------------------------------------------------------
  ! For Newtonian fluid, pressure using Neumann conditon.
  !----------------------------------------------------------

  IF ( Newtonian ) THEN
     
     this%p(ip) = this%p(jp)
     
     CALL rhs_force_ff(this%rhs,&
          this%x(1:num_dim,ip),&
          this%x(1:num_dim,jp),dij,&
          v_ip(1:num_dim),&
          v_jp(1:num_dim),&
          this%rho(ip),this%rho(jp),&
          this%p(ip),this%p(jp),&
          this%m(ip),this%m(jp),&
          w,gradw,fip,fjp,stat_info=stat_info_sub)
     
  ELSE ! non-Newnotian
     
     !-------------------------------------------------------
     ! For non-Newtonian fluid, density and pressure 
     ! using Neumann conditon.
     !-------------------------------------------------------
     
     rho_ip = this%rho(jp)
     rho_jp = this%rho(jp)
     
#include "pp_vgt_ip.inc"
     
     IF ( symmetry ) THEN
        
#include "pp_vgt_jp.inc"
        
     END IF
     
     this%pt(1:num_dim,1:num_dim,ip) = &
          this%pt(1:num_dim,1:num_dim,jp)
     
     CALL rhs_force_ff(this%rhs,&
          this%x(1:num_dim,ip),&
          this%x(1:num_dim,jp),dij,&
          v_ip(1:num_dim),&
          v_jp(1:num_dim),&
          rho_ip,rho_jp,&
          this%pt(1:num_dim,1:num_dim,ip), &
          this%pt(1:num_dim,1:num_dim,jp), &
          this%m(ip),this%m(jp),&
          w,gradw,fip,fjp,stat_info=stat_info_sub)
     
  END IF ! Newtonian
  
  IF( stat_info_sub /=0 ) THEN
     PRINT *, "pp_interactoin_cf : ",& 
          "Pairwise force has problem !"
     stat_info = -1
     GOTO 9999
  END IF
