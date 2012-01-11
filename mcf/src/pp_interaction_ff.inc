  !------------------------------------------------
  ! particle_particle_interaction_fluid_fluid
  !
  ! Remark:
  ! 
  ! In case of non-Newtonian fluid,
  ! we calculate velocity gradient tensor;
  ! For force calculation, we use pressure tensor.
  !------------------------------------------------
  
 
  v_ip(1:num_dim) = this%v(1:num_dim,ip)
  v_jp(1:num_dim) = this%v(1:num_dim,jp)
  
  IF ( Newtonian ) THEN
     
     CALL rhs_force_ff(this%rhs, &
          this%x(1:num_dim,ip), &
          this%x(1:num_dim,jp),dij, &
          v_ip(1:num_dim), v_jp(1:num_dim), &
          this%rho(ip),this%rho(jp), &
          this%p(ip),this%p(jp), &
          this%m(ip),this%m(jp), &
          w,gradw,fip,fjp,stat_info=stat_info_sub)
     
  ELSE ! non-Newtonian
     
     rho_ip = this%rho(ip)
     rho_jp = this%rho(jp)
     
#include "pp_vgt_ip.inc"
     
     IF ( symmetry ) THEN
        
#include "pp_vgt_jp.inc"
        
     END IF

     CALL rhs_force_ff(this%rhs,&
          this%x(1:num_dim,ip),&
          this%x(1:num_dim,jp),dij,&
          v_ip(1:num_dim), &
          v_jp(1:num_dim), &
          rho_ip,rho_jp,&
          this%pt(1:num_dim,1:num_dim,ip), &
          this%pt(1:num_dim,1:num_dim,jp),&
          this%m(ip),this%m(jp),&
          w,gradw,fip,fjp,stat_info=stat_info_sub)
     
  END IF
  
  IF(stat_info_sub /=0) THEN
     PRINT *, "pp_interaction_ff : ", & 
          "Pairwise force has problem !"
     stat_info = -1
     GOTO 9999
  END IF
  
