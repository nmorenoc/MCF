      SUBROUTINE particles_compute_act(this,num,stat_info)
        !----------------------------------------------------
        ! Subroutine :  particles_compute_act
        !----------------------------------------------------
        !
        ! Purpose   :  Compute the accleration of conformation
        !              tensor of particles.
        !
        ! Reference : Vazquez-Quesada, Ellero and Espanol
        !             Phyical Review E 79. 056707, 2009.
        !
        ! Remark    :
        !
        ! Revision  :  V0.1  31.07.2009, original version.
        !
        !              15.04.2010, Adolfo Vázquez-Quesada
        !                 I changed the coefficients of the 
        !                 velocity gradient tensor in act
        !----------------------------------------------------
        ! Author    : Xin Bian
        ! Contact   : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        !----------------------------------------------------
        ! Arguments
        !
        ! this           : an object of Particles Class.
        ! num            : number of particles updated,
        !                  i.e. first num particles in this%x 
        !                  are operated.
        ! lambda         : coefficient required.
        ! stat_info      : return flag of status.
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)          :: this
        INTEGER, INTENT(IN)                     :: num
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables
        !----------------------------------------------------
        
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: dim
        REAL(MK)                                :: tau
        INTEGER                                 :: i,j,k,p
        REAL(MK), DIMENSION(3,3)                :: t_c
        REAL(MK), DIMENSION(3,3)                :: t_vgt
        REAL(MK), DIMENSION(3,3)                :: Im
        REAL(MK)                                :: t_act

        !-------------------------------
        ! Initialization of variables.
        !-------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        !--------------------------------
        ! Calculation only for real particles.
        !--------------------------------
        
        IF( num > this%num_part_real) THEN
           PRINT *, "particles_compute_act : ", &
                "num > num_part_real, wrong !"
           stat_info = -1
           GOTO 9999      
        END IF
        
        dim = physics_get_num_dim(this%phys,stat_info_sub)
        tau = physics_get_tau(this%phys,stat_info_sub)
        
        IF (ABS(tau) < mcf_machine_zero) THEN
           PRINT *, "particles_compute_act : ", &
                "tau should not be zero !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !-------------------------------
        ! Allocate memory for act.
        !-------------------------------
        
        IF(ASSOCIATED(this%act)) THEN
           DEALLOCATE(this%act)
        ENd IF
        
        ALLOCATE(this%act(dim**2,num))
        
        this%act(:,:) = 0.0_MK
        
        !----------------------
        ! Identity matrix
        !----------------------
        
        Im(1:3,1:3) = 0.0_MK
        Im(1,1) = 1.0_MK
        Im(2,2) = 1.0_MK
        Im(3,3) = 1.0_MK
        
        !------------------------------------------
        ! Compute accleration of conformation 
        ! tensor using velocity gradient tensor.
        !------------------------------------------
        
        DO p = 1, num
           
           !---------------------------------------
           ! Convert array notation to matrix
           ! for clarity.
           !---------------------------------------
           
           DO j = 1, dim              
              DO i = 1, dim
                 t_c(i,j)   = this%ct(i+dim*(j-1),p)
                 t_vgt(i,j) = this%vgt(i+dim*(j-1),p)
              END DO
           END DO
                      
           DO j =1, dim  !--- row direction
              
              Do i=1, dim ! | column direction
                 
                 t_act = 0.0_MK
                 
                 DO k= 1, dim
                    
                    !------------------------------
                    ! c <dot> k + k^t <dot> c
                    !------------------------------
                    
!!$                    t_act = t_act + &
!!$                         t_c(i,k) * t_vgt(k,j) + &
!!$                         t_vgt(k,i) * t_c(k,j)

                    t_act = t_act + &
                         t_c(i,k) * t_vgt(j,k) + &
                         t_vgt(i,k) * t_c(k,j)
                    
                 END DO ! k
                 
                 
                 this%act(i+dim*(j-1),p) = t_act - &
                      ( t_c(i,j) - Im(i,j) ) /tau                    
                 
              END DO ! i
              
           END DO ! j
           
        END DO ! p
        
9999    CONTINUE      
        
        RETURN
        
      END SUBROUTINE particles_compute_act
      
      
      
