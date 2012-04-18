      SUBROUTINE particles_apply_body_force(this,num,stat_info)
        !----------------------------------------------------
        !  Subroutine :  particles_apply_body_force
        !----------------------------------------------------
        !
        !  Purpose    : Apply external / body force to first
        !               num particles.
        !
        !  Reference  :
        !
        !  Remark     :
        !
        !
        !  Revisions  : V0.3 30.07.2009, adding one more type 
        !               of body force. F*Sin(ky) for Kolmogorov
        !
        !               V0.3 27.07.2009, if there is a body force,
        !               it should apply on colloid also,
        !               instead of applying only on fluid particles.
        !
        !               V0.2 08.07.2009, 
        !               check again the work flow is correct and
        !               supply with more comments for code.
        !
        !               V0.1 01.03.2009, original version.
        !
        !----------------------------------------------------
        ! Author       : Xin Bian
        ! Contact      : xin.bian@aer.mw.tum.de
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
        ! num        : first num particles needed body force.
        ! stat_info  : return flag of status.
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)  :: this
        INTEGER, INTENT(IN)             :: num
        INTEGER, INTENT(OUT)            :: stat_info
        
        !------------------------------------------
      	! Local variables
      	!------------------------------------------
      
        INTEGER                         :: stat_info_sub
        INTEGER                         :: body_force_type
        
        
#ifdef __DEBUG
        INTEGER                         :: rank
        INTEGER                         :: debug_flag
        REAL(MK)                        :: time_routine_start
#endif
        
        !------------------------------------------
        ! Initialization of variables.
        !------------------------------------------
        
        stat_info      = 0
        stat_info_sub  = 0        
        
        
#ifdef __DEBUG
                
        debug_flag  = debug_get_flag(global_debug,stat_info_sub)
        rank  = technique_get_rank(this%tech,stat_info_sub)           
        IF( debug_flag == 2 ) THEN
           CALL debug_substart(global_debug,rank,&
                "particles_apply_body_force",&
                time_routine_start,stat_info_sub)
        END IF
#endif
        
        body_force_type = &
             physics_get_body_force_type(this%phys,stat_info_sub)
        
        !------------------------------------------
        ! For normal body force, first type
        ! should be used;
        !
        ! For poiseuille in periodic boundaries,
        ! 2 direction of body force can be used.
        !
        ! For oscillating shear with sin function
        ! then 3rd should be used.
        !------------------------------------------
        
        SELECT CASE(body_force_type)
           
        CASE (0)
        CASE (1)
           CALL particles_apply_body_force_1direction(this,num,stat_info_sub)
        CASE (2)
           CALL particles_apply_body_force_2direction(this,num,stat_info_sub)
        CASE (3)
           CALL particles_apply_body_force_sin(this,num,stat_info_sub)
	CASE (4)
	   CALL particles_apply_body_force_sm(this,num,stat_info_sub)
           
        END SELECT
        
        
#ifdef __DEBUG
        IF( debug_flag == 2 ) THEN
           CALL debug_substop(global_debug,rank,&
                "particles_apply_body_force",&
                time_routine_start,stat_info_sub)
        END IF
#endif
        
        RETURN
        
      END SUBROUTINE particles_apply_body_force
      
      
      
      SUBROUTINE particles_apply_body_force_1direction(this,num,stat_info)
        !----------------------------------------------------
        !  Subroutine :  particles_apply_body_force_1direction
        !----------------------------------------------------
        !
        !  Purpose    :  The subroutines for applying body force
        !                 only for fluid particles.
        !	 	      	 
        !                
        !  Remarks    :
        !
        !  References :
        !
        !  Revisions  : 0.1 08.02.2009, original version.
        !----------------------------------------------------
        !  Author       : Xin Bian
        !  Contact      : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)          :: this
        INTEGER, INTENT(IN)                     :: num
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: num_dim
        REAL(MK), DIMENSION(:), POINTER         :: body_force
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: i,j
        
        
        stat_info     = 0
        stat_info_sub = 0
        NULLIFY(body_force)
        
        num_dim  = physics_get_num_dim(this%phys,stat_info_sub)
        
        CALL physics_get_body_force(this%phys,body_force,stat_info_sub)
        
        DO i = 1, num_dim
           
           IF( ABS(body_force(i)) > mcf_machine_zero ) THEN
              
              DO j = 1, num
                 
                 IF (this%id(this%sid_idx,j) == mcf_particle_type_fluid ) THEN
                    
                    this%f(i,j) = &
                         this%f(i,j) + body_force(i)
                    
                 END IF
                 
              END DO
              
           END IF
           
        END Do
        
        
        IF( ASSOCIATED(body_force) ) THEN
           DEALLOCATE(body_force)
        END IF
        
        RETURN
        
      END SUBROUTINE particles_apply_body_force_1direction
      
        
      SUBROUTINE particles_apply_body_force_2direction(this, num,stat_info)
        !--------------------------------------------------
        ! Subroutine  :  particls_apply_body_force_sin
        !--------------------------------------------------
        !
        ! Purpose     :  Applying body force in two direction,
        !                x+, x-.
        !                
        !	 	      	 
        !                
        !  Remarks      :
        !
        !  References   :
        !
        !  Revisions    : 0.1 08.02.2009 
        !-------------------------------------------------------------
        !  Author       : Xin Bian
        !  Contact      :  xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !-------------------------------------------------------------
          
        TYPE(Particles), INTENT(INOUT)          :: this
        INTEGER, INTENT(IN)                     :: num
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: num_dim
        REAL(MK), DIMENSION(:), POINTER         :: max_phys
        REAL(MK), DIMENSION(:), POINTER         :: body_force
        
        INTEGER                                 :: stat_info_sub
        INTEGER	 		                :: i
        
        stat_info     = 0
        stat_info_sub = 0
        NULLIFY(max_phys)
        NULLIFY(body_force)
        
        num_dim  = physics_get_num_dim(this%phys,stat_info_sub)
        CALL physics_get_max_phys(this%phys,max_phys,stat_info_sub)        
        CALL physics_get_body_force(this%phys,body_force,stat_info_sub)
        
        DO i =1, num
           
           IF( ABS(body_force(1)) > mcf_machine_zero) THEN
              
              IF ( this%x(num_dim, i) < max_phys(num_dim) * 0.5_MK .AND. &
                   this%id(this%sid_idx, i) == mcf_particle_type_fluid ) THEN
                 
                 this%f(1,i) = this%f(1,i) +&
                      body_force(1)
                 
              ELSE IF (this%x(num_dim,i) > max_phys(num_dim) * 0.5_MK .AND. &
                   this%id(this%sid_idx, i) == mcf_particle_type_fluid ) THEN
                 
                 this%f(1,i) = this%f(1,i) - &
                      body_force(1)
                 
              END IF
              
           END IF
           
        END DO
        
        IF ( ASSOCIATED(max_phys)) THEN
           DEALLOCATE(max_phys)
        END IF
        
        IF( ASSOCIATED(body_force) ) THEN
           DEALLOCATE(body_force)
        END IF
        
        RETURN
        
      END SUBROUTINE particles_apply_body_force_2direction
        
        
      SUBROUTINE particles_apply_body_force_sin(this,num,stat_info)
        !----------------------------------------------------
        !  Subroutine	:  particles_apply_body_force_sin
        !----------------------------------------------------
        !
        !  Purpose      :  Applying body force  as sin() function 
        !                  along last dimension in x-direction.
        !	 	      	 
        !                
        !  Remarks      :
        !
        !  References   :
        !
        !  Revisions    : 0.1 30.07.2009 
        !----------------------------------------------------
        !  Author       : Xin Bian
        !  Contact      : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)          :: this
        INTEGER, INTENT(IN)                     :: num
        INTEGER, INTENT(OUT)                    :: stat_info
        
        REAL(MK), DIMENSION(:), POINTER         :: body_force
        REAL(MK), DIMENSION(:), POINTER         :: min_phys
        REAL(MK), DIMENSION(:), POINTER         :: max_phys
        
        
        REAL(MK)                                :: i, k
        INTEGER                                 :: stat_info_sub
        
        stat_info     = 0
        stat_info_sub = 0
        NULLIFY(min_phys)
        NULLIFY(max_phys)          
        NULLIFY(body_force)
        
        CALL physics_get_min_phys(this%phys,min_phys,stat_info_sub)
        CALL physics_get_max_phys(this%phys,max_phys,stat_info_sub)
        CALL physics_get_body_force(this%phys,body_force,stat_info_sub)
        
        
        k = 2.0_MK* mcf_pi / &
             ( max_phys(2) - min_phys(2) )
        
        IF( ABS(body_force(1)) > mcf_machine_zero) THEN
           
           DO i =1, num
              
              IF (this%id(this%sid_idx,i) == mcf_particle_type_fluid ) THEN
                 this%f(1,i) =  this%f(1,i) + &
                      body_force(1) * sin( k * this%x(2,i))
              END IF
              
           END DO
           
        END IF
        
        IF( ASSOCIATED(min_phys) ) THEN
           DEALLOCATE(min_phys)
        END IF
        
        IF( ASSOCIATED(max_phys) ) THEN
           DEALLOCATE(max_phys)
        END IF
        
        IF( ASSOCIATED(body_force) ) THEN
           DEALLOCATE(body_force)
        END IF
        
        RETURN
        
      END SUBROUTINE particles_apply_body_force_sin


     SUBROUTINE particles_apply_body_force_sm(this,num,stat_info)

        !----------------------------------------------------
        !  Subroutine	:  particles_apply_body_force_sm
        !----------------------------------------------------
        !
        !  Purpose      :  first test case of stochastic model 
        !                  external force for simple Langevin equation
        !	 	      	 
        !                
        !  Remarks      :
        !
        !  References   :
        !
        !  Revisions    : 02.11.2010 
        !----------------------------------------------------

        
        TYPE(Particles), INTENT(INOUT)          :: this
        INTEGER, INTENT(IN)                     :: num
        INTEGER, INTENT(OUT)                    :: stat_info
        
        INTEGER                                 :: num_dim
	INTEGER                                 :: i
        INTEGER                                 :: j         
        
        REAL(MK)                                :: k_sm
	REAL(MK)                                :: tau_sm
	REAL(MK)                                :: dt
	REAL(MK)                                :: gauss_num
        INTEGER                                 :: stat_info_sub
        
        stat_info     = 0
        stat_info_sub = 0
        
        num_dim  = physics_get_num_dim(this%phys,stat_info_sub)  
        tau_sm = physics_get_tau_sm(this%phys,stat_info_sub)
	k_sm = physics_get_k_sm(this%phys,stat_info_sub)
	dt = physics_get_dt(this%phys,stat_info_sub)

	DO i = 1,num_dim
	   DO j = 1,num
	      gauss_num = random_random(this%random,stat_info_sub)
!	      PRINT *, "gauss_num : ", gauss_num
	      this%f(i,j) =  this%f(i,j) - &
                	  1.0_MK / tau_sm * this%v(i,j) + sqrt(4.0_MK * k_sm / 3.0_MK / tau_sm / dt) * gauss_num
	   END DO
	END DO

      END SUBROUTINE particles_apply_body_force_sm
