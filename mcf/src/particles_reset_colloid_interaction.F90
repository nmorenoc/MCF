      SUBROUTINE particles_reset_colloid_interaction(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_reset_colloid_interaction
        !----------------------------------------------------
        !
        ! Purpose     : Reset colloidal boundary particles
        !               interaction, e.g., acceleration
        !               to its instaneous one according
        !               to its host colloid. 
        !
        ! Refernece   :
        !
        ! Remark      :
        !
        ! Revision    : V0.3 02.12.2009, including 
        !               Lees-Edwards boundary.
        !
        !               V0.2 08.07.2009, Check work flow
        !               is correct and supply with comments.
        !          
        !               V0.1 24.06.2009, original version.
        !
        !----------------------------------------------------
        ! Author      : Xin Bian
        ! Contact     : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------


        !----------------------------------------------------
        ! Arguments
        !
        ! this       : an object of Particles Class.
        ! stat_info  : return flag of status.
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)          :: this
        INTEGER, INTENT(OUT)                    :: stat_info

        
	!----------------------------------------------------
        ! Local variables
	!----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: dim, num_colloid
        TYPE(Colloid), POINTER                  :: colloids
        REAL(MK), DIMENSION(:,:,:), POINTER     :: coll_f
        REAL(MK), DIMENSION(:,:,:), POINTER     :: coll_alpha

        REAL(MK), DIMENSION(3)                  :: rx, rf
        INTEGER                                 :: i,ip,sid
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        dim         = &
             physics_get_num_dim(this%phys,stat_info_sub)
        num_colloid = &
             physics_get_num_colloid(this%phys,stat_info_sub)
        
        NULLIFY(colloids)
        NULLIFY(coll_f)
        NULLIFY(coll_alpha)
        
        IF ( num_colloid > 0 ) THEN
           
           CALL physics_get_colloid(this%phys,colloids,stat_info_sub)
           
           CALL colloid_get_f(colloids,coll_f,stat_info_sub)
           CALL colloid_get_alpha(colloids,coll_alpha,stat_info_sub)
        
           DO i =1, this%num_part_colloid
              
              !----------------------------------------------
              ! Get index of this boundary particle
              ! and its species ID.
              !----------------------------------------------
           
              ip  = this%part_colloid_list(1,i)
              sid = this%part_colloid_list(2,i)
              
              this%f(1:dim,ip) = coll_f(1:dim,sid,1)
              
              CALL colloid_in_relative_position(colloids,&
                   this%x(1:dim,ip),sid,rx(1:dim),stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *, "particles_reset_colloid_interaction: ", &
                      "Calculating colloid particle position failed !"
                 stat_info = -1
                 GOTO 9999
              END IF
              
              CALL tool_cross_product(this%tool,&
                   coll_alpha(1:3,sid,1),rx(1:3), &
                   rf(1:3),stat_info_sub)
              
              this%f(1:dim,ip) = this%f(1:dim,ip) + &
                   rf(1:dim)
              
           END DO ! i = 1, num_part_colloid
        
        END IF ! num_colloid > 0

        
9999    CONTINUE
        
        IF(ASSOCIATED(coll_f)) THEN
           DEALLOCATE(coll_f)
        END IF
        
        IF(ASSOCIATED(coll_alpha)) THEN
           DEALLOCATE(coll_alpha)
        END IF
        
        RETURN
        
      END SUBROUTINE particles_reset_colloid_interaction
      
      
