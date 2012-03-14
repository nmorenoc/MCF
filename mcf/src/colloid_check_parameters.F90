      LOGICAL FUNCTION colloid_check_parameters(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_check_parameters
        !----------------------------------------------------
        !
        ! Purpose     : Check if colloid parameters are
        !               given resonably.
        !
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revisions   :  V0.1 21.10.2009, original version.
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
        
        TYPE(colloid), INTENT(IN)       :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        
        LOGICAL                         :: lcheck
        INTEGER                         :: num_dim


        INTEGER                         :: i,j
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info = 0
        
        lcheck    = .TRUE.
        num_dim   = this%num_dim
        
        !----------------------------------------------------
        ! Check lubrication cut off or on for colloid-colloid.
        !----------------------------------------------------

        IF ( this%cc_lub_type /= 0 .AND. &
             this%cc_lub_cut_off < mcf_machine_zero) THEN
           PRINT *, "colloid_check_parameters : ",&
                "cc_lub_cut_off should be positive ! "
           lcheck = .FALSE.
           GOTO 9999
        END IF
        
        IF ( this%cc_lub_type /= 0 .AND. &
             this%cc_lub_cut_on < mcf_machine_zero) THEN
           PRINT *, "colloid_check_parameters : ",&
                "cc_lub_cut_on should be positive ! "
           lcheck = .FALSE.
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Check repulsive force cut off for colloid-colloid.
        !----------------------------------------------------

        IF ( this%cc_repul_type /= 0 .AND. &
             this%cc_repul_cut_off < mcf_machine_zero) THEN
           PRINT *, "colloid_check_parameters : ",&
                "cc_repul_cut_off should be positive ! "
           lcheck = .FALSE.
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Check repulsive force cut on for colloid-colloid.
        !----------------------------------------------------

        IF ( this%cc_repul_type /= 0 .AND. &
             this%cc_repul_cut_on < mcf_machine_zero) THEN
           PRINT *, "colloid_check_parameters : ",&
                "cc_repul_cut_on should be positive ! "
           lcheck = .FALSE.
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Check maximum repulsive force for collod-colloid
        !----------------------------------------------------

        IF ( this%cc_repul_type /= 0 .AND. &
             this%cc_repul_F0 < mcf_machine_zero) THEN
           PRINT *, "colloid_check_parameters : ",&
                "cc_repul_F0 should be positive ! "
           lcheck = .FALSE.
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Check lubrication cut off or on for colloid-wall.
        !----------------------------------------------------

        IF ( this%cw_lub_type /= 0 .AND. &
             this%cw_lub_cut_off < mcf_machine_zero) THEN
           PRINT *, "colloid_check_parameters : ",&
                "cw_lub_cut_off should be positive ! "
           lcheck = .FALSE.
           GOTO 9999
        END IF
        
        IF ( this%cw_lub_type /= 0 .AND. &
             this%cw_lub_cut_on < mcf_machine_zero) THEN
           PRINT *, "colloid_check_parameters : ",&
                "cw_lub_cut_on should be positive ! "
           lcheck = .FALSE.
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Check repulsive force cut off between colloid-wall.
        !----------------------------------------------------
        
        IF ( this%cw_repul_type /= 0 .AND. &
             this%cw_repul_cut_off < mcf_machine_zero) THEN
           PRINT *, "colloid_check_parameters : ",&
                "cw_repul_cut_off should be positive ! "
           lcheck = .FALSE.
           GOTO 9999
        END IF

        !----------------------------------------------------
        ! Check repulsive force cut on between colloid-wall.
        !----------------------------------------------------
        
        IF ( this%cw_repul_type /= 0 .AND. &
             this%cw_repul_cut_on < mcf_machine_zero) THEN
           PRINT *, "colloid_check_parameters : ",&
                "cw_repul_cut_on should be positive ! "
           lcheck = .FALSE.
           GOTO 9999
        END IF

      
        !----------------------------------------------------
        ! Check maximum repulsive force between colloid-wall.
        !----------------------------------------------------
        
        IF ( this%cw_repul_type /= 0 .AND. &
             this%cw_repul_F0 < mcf_machine_zero) THEN
           PRINT *, "colloid_check_parameters : ",&
                "cw_repul_F0 should be positive ! "
           lcheck = .FALSE.
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Loop over all colloids:
        !----------------------------------------------------
        
        DO j = 1, this%num_colloid
           
           !-------------------------------------------------
           ! normally colloids' center must be given 
           ! inside the domain,
           ! except it is symmetry boundary on that side.
           !-------------------------------------------------
           
           DO i = 1, this%num_dim
              
              IF ( this%x(i,j) < this%min_phys(i) ) THEN
                 
                 PRINT *,&
                      "colloid_check_parameters : ", &
                      "center should be inside domain for colloid ", j
                 lcheck = .FALSE.
                 GOTO 9999
                 
              END IF
              
              IF ( this%x(i,j) > this%max_phys(i) ) THEN
                 
                 PRINT *,&
                      "colloid_check_parameters : ", &
                      "center should be inside domain for colloid ", j
                 
                 lcheck = .FALSE.
                 GOTO 9999
                 
              END IF ! 
              
           END DO ! i = 1, num_dim
           
           
           !-------------------------------------------------
           ! The distance from surface to the center should
           ! be bigger than din (approximately cut off, 
           ! a little bit bigger for caution).
           !-------------------------------------------------
           
           SELECT CASE( this%shape(j) ) 
              
              !----------------------------------------------
              ! cylinder, disk/sphere
              !----------------------------------------------
              
           CASE (mcf_colloid_shape_cylinder)
              
              DO i = 1, num_dim
                 
                 IF ( this%radius(i,j) < mcf_machine_zero ) THEN
                    PRINT *, "colloid_check_parameters: ", &
                         "radius must be bigger than 0 !"
                    stat_info = -1
                    GOTO 9999
                 END IF
                 
              END DO

              IF ( num_dim == 3 .AND. &
                   this%radius(1,j) < this%radius(2,j) ) THEN
                 
                 PRINT *, "colloid_check_parameters: ", &
                      "bigger radius must be in front !"
                 stat_info = -1
                 GOTO 9999
                 
              END IF
              
           CASE (mcf_colloid_shape_sphere)
              
              IF ( this%radius(1,j) < mcf_machine_zero ) THEN
                 PRINT *, "colloid_check_parameters: ", &
                      "radius must be bigger than 0 !"
                 stat_info = -1
                 GOTO 9999
              END IF
              
              
           CASE (mcf_colloid_shape_ellipsoid)
              
              DO i = 1, num_dim
                 
                 IF ( this%radius(i,j) < mcf_machine_zero ) THEN
                    PRINT *, "colloid_check_parameters: ", &
                         "radius must be bigger than 0 !"
                    stat_info = -1
                    GOTO 9999
                 END IF
                 
              END DO

              DO i = 1, num_dim - 1
                 
                 IF ( this%radius(i,j) < this%radius(i+1,j) ) THEN
                    PRINT *, "colloid_check_parameters: ", &
                         "bigger radius must be in front !"
                    stat_info = -1
                    GOTO 9999
                 END IF
                 
              END DO

           CASE (mcf_colloid_shape_dicolloid)
              
              DO i = 1, num_dim
                 IF ( this%radius(i,j) < mcf_machine_zero ) THEN
                    PRINT *, "colloid_check_parameters: ", &
                         "radius must be bigger than 0 !"
                    stat_info = -1
                    GOTO 9999
                 END IF
              END DO
              
              DO i = 1, num_dim - 1
                 
                 IF ( this%radius(i,j) < this%radius(i+1,j) ) THEN
                    PRINT *, "colloid_check_parameters: ", &
                         "bigger radius must be in front !"
                    stat_info = -1
                    GOTO 9999
                 END IF
                 
              END DO
              
           CASE (mcf_colloid_shape_star)
              
              DO i = 1, num_dim
                 IF ( this%radius(i,j) < mcf_machine_zero ) THEN
                    PRINT *, "colloid_check_parameters: ", &
                         "radius must be bigger than 0 !"
                    stat_info = -1
                    GOTO 9999
                 END IF
              END DO
              
              IF ( this%radius(2,j) < 0.0_MK ) THEN
                 PRINT *, "colloid_check_parameters : Radius (" , &
                      this%radius(2,j), ") of colloid ", j, &
                      "should be bigger than 0 "
                 lcheck = .FALSE.
                 GOTO 9999
              END IF
              
              IF ( this%freq(j) <= 0 ) THEN
                 
                 PRINT *, "colloid_check_parameters : Freq (" , &
                      this%freq(j), ") of colloid ", j, &
                      "should be bigger than 0 "
                 lcheck = .FALSE.
                 GOTO 9999
              END IF
              
           CASE DEFAULT
              
              PRINT *, "colloid_check_parameters: ", &
                   "No shape available !"
              stat_info = -1
              GOTO 9999
              
           END SELECT ! shape(j)
           
           
        END DO ! j = 1, num_colloid


9999    CONTINUE
        
        colloid_check_parameters = lcheck
        
        RETURN
        
      END FUNCTION colloid_check_parameters
      
