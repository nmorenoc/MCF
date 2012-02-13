  !----------------------------------------------------------
  ! This file is used/included as "subroutine"/"repeated part"
  ! particle-particle-interaction :
  !
  !
  ! Remark 1:
  !
  ! This .inc file is included by 
  ! particles_compute_interaction() routine,
  !
  ! Remark 2: 
  !
  ! f : fluid particle.
  ! c : colloidal boundary particle.
  ! w : wall boundary particle.
  !
  ! Posibility of interaction pairs :
  !
  ! f - f  Yes
  ! f - c  Yes
  ! f - w  Yes
  !
  ! c - f  Yes
  ! c - c  depending on pp_interact_cc
  ! c - w  depending on pp_interact_cw
  !
  ! w - f  Yes
  ! w - c  depending on pp_interact_cw
  ! w - w  no
  !
  !----------------------------------------------------------
  
  !----------------------------------------------------------
  ! Distance of particles ip and jp.
  !----------------------------------------------------------
  
  rij(1:num_dim) = &
       this%x(1:num_dim,ip) - this%x(1:num_dim,jp)
  
  dij = DOT_PRODUCT(rij(1:num_dim), rij(1:num_dim))
  
  !----------------------------------------------------------
  ! Skip 2 particles beyond cuf off.
  !----------------------------------------------------------
  
  IF ( dij >= cut_off2 ) THEN
     CYCLE
  ELSE
     dij = SQRT(dij)
  END IF
  
  !----------------------------------------------------------
  ! Particle ip interacts with jp.
  ! Get kernel value and its derivative first.
  !----------------------------------------------------------
  
  CALL kernel_kernel(this%kern,dij,w,gradw,stat_info_sub)
  
  !----------------------------------------------------------
  ! Check type of particle ip.
  !----------------------------------------------------------
  
  SELECT CASE( this%id(this%sid_idx,ip) ) 
     
     !-------------------------------------------------------
     ! ip is fluid particle. (ip-jp interaction)
     !-------------------------------------------------------
     
  CASE (mcf_particle_type_fluid)
     
     !-------------------------------------------------------
     ! Check type of particle jp. 
     !-------------------------------------------------------
     
     SELECT CASE( this%id(this%sid_idx,jp) )
        
        !----------------------------------------------------
        ! jp is fluid particle. f-f interaction.
        !----------------------------------------------------
        
     CASE (mcf_particle_type_fluid)
        
#include "pp_interaction_ff.inc"
        
        !----------------------------------------------------
        ! jp is colloidal boundary particle. f-c interaction.
        !----------------------------------------------------
        
     CASE (mcf_particle_type_colloid:)
        
#include "pp_interaction_fc.inc"
        
        !----------------------------------------------------
        ! jp is wall boundary particle. f-w interaction. 
        !----------------------------------------------------
        
     CASE(:mcf_particle_type_wall)
        
#include "pp_interaction_fw.inc"
        
     CASE DEFAULT
        
        PRINT *, "pp_interaction : ", &
             "interaction pair does not exist !"
        stat_info = -1
        GOTO 9999
        
     END SELECT ! this%id(this%sid_idx,jp)
     
     !-------------------------------------------------------
     ! ip is colloidal boundary particle. (ip-jp interaction)
     !-------------------------------------------------------
     
  CASE (mcf_particle_type_colloid:) 
     
     !-------------------------------------------------------
     ! Check type of particle jp.
     !-------------------------------------------------------
     
     SELECT CASE (this%id(this%sid_idx,jp))
        
        !----------------------------------------------------
        ! jp is fluid particle. c-f interaction.
        !----------------------------------------------------
        
     CASE (mcf_particle_type_fluid)
        
#include "pp_interaction_cf.inc"                             
        
        !----------------------------------------------------
        ! jp is colloidal boundary particle, c-c interaction.
        ! ip and jp must be from different colloid.        
        !----------------------------------------------------
        
     CASE (mcf_particle_type_colloid:)
        
        IF ( this%id(this%sid_idx,jp) == &
             this%id(this%sid_idx,ip) ) THEN
           
           CYCLE
           
        ELSE IF ( this%pp_interact_cc ) THEN
           
#include "pp_interaction_cc.inc"
           
        END IF
        
        !----------------------------------------------------
        ! jp is wall boundary particle. c-w interaction. 
        !----------------------------------------------------
        
     CASE(:mcf_particle_type_wall)
        
        IF ( this%pp_interact_cw ) THEN
           
#include "pp_interaction_cw.inc"
           
        END IF
        
     CASE DEFAULT
        
        PRINT *, "pp_interaction : ", &
             "interaction pair does not exist !"
        stat_info = -1
        GOTO 9999
        
     END SELECT ! this%id(this%sid_idx,jp)
     
     !-------------------------------------------------------
     ! ip is wall boundary particle. (ip-jp interaction)
     !-------------------------------------------------------
     
  CASE (:mcf_particle_type_wall)
     
     SELECT CASE (this%id(this%sid_idx,jp))
        
        !----------------------------------------------------
        ! jp is fluid particle. w-f interaction.
        !----------------------------------------------------
        
     CASE (mcf_particle_type_fluid)
        
#include "pp_interaction_wf.inc"
        
        !----------------------------------------------------
        ! jp is colloid boundary particle.  w-c interaction. 
        !----------------------------------------------------
        
     CASE (mcf_particle_type_colloid:)
        
        IF ( this%pp_interact_cw ) THEN
           
#include "pp_interaction_wc.inc"
        
        END IF
        
        !----------------------------------------------------
        ! jp is wall boundary particle. w-w interaction, 
        ! i.e., no interaction at all.
        !----------------------------------------------------
        
     CASE (:mcf_particle_type_wall)
        
        CYCLE
        
     CASE DEFAULT
        
        PRINT *, "pp_interaction : ", &
             "interaction pair does not exist !"
        stat_info = -1
        GOTO 9999
        
     END SELECT ! this%id(this%sid_idx,jp)
     
     
  CASE DEFAULT
     
     PRINT *, "pp_interaction : ", &
          "interaction pair does not exist !"
     stat_info = -1
     GOTO 9999
     
  END SELECT ! this%id(this%sid_idx,ip)
