      SUBROUTINE particles_compute_evgt(this,num,stat_info)
        !----------------------------------------------------
        ! Subroutine :  particles_compute_evgt
        !----------------------------------------------------
        !
        ! Purpose    :  Compute egenvector dynamics matrix 
        !              elemetens from velocity gradient tensor.
        !
        ! Reference  : Vazquez-Quesada, Ellero, Espanol
        !              Phyical Review E 79. 056707, 2009.
        !
        ! Remark     :
        !
        ! Revision   :  V0.1  20.08.2010, an innecesary line
        !               to assign the 0 value to this%evgt 
        !               was deleted. (Adolfo)
        !
        !               V0.1  04.08.2009, original version.
        !
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
        ! num            : number of particles needed to be 
        !                  updated, i.e. first num particles 
        !                  in this%x  are operated.
        ! lambda         : coefficient required.
        ! stat_info      : return flag of status.
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)          :: this
        INTEGER, INTENT(IN)                     :: num
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------
        ! Local variables
	!----------------------

        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: dim
        INTEGER                                 :: i,j,p
        REAL(MK), DIMENSION(3,3)                :: t_vgt
        REAL(MK), DIMENSION(3,3)                :: t_evec

        !-------------------------------
        ! Initialization of variables.
        !-------------------------------
        
        stat_info     = 0
        stat_info_sub = 0

        !--------------------------------
        ! Calculation can be for ghosts
        ! particles also.
        !--------------------------------
        
        IF( num > this%num_part_all) THEN
           PRINT *, "particles_compute_evgt : ", &
                "num > num_part_all, wrong !"
           stat_info = -1
           GOTO 9999      
        END IF
        
        dim = &
             physics_get_num_dim(this%phys,stat_info_sub)
        
        
        IF(ASSOCIATED(this%evgt)) THEN
           DEALLOCATE(this%evgt)
        ENd IF
        
        ALLOCATE(this%evgt(dim**2,num))
!!$        this%evgt(:,:) = 0.0_MK
        
        !------------------------------------------
        ! Compute martrix elements of egenvector
        ! dynamics velocity gradient tensor.
        !------------------------------------------
        
        
        DO p = 1, num ! particle index
           
           !---------------------------------------
           ! Convert array notation to matrix
           ! notation for clarity.
           !---------------------------------------
           
           DO j = 1, dim
              DO i = 1, dim
                 t_vgt(i,j)  = this%vgt(i+dim*(j-1),p)
                 t_evec(i,j) = this%evec(i+dim*(j-1),p)
              END DO
           END DO
           
           DO j =1, dim ! --- row direction
              
              Do i=1, dim ! | column direction
                 
                 this%evgt(i+dim*(j-1),p) = &
                      DOT_PRODUCT(MATMUL(t_evec(1:dim,i),&
                      t_vgt(1:dim,1:dim)), &
                      t_evec(1:dim,j))
                 
              END DO ! i
              
           END DO ! j
           
        END DO ! p
        
        
9999    CONTINUE      
        
        RETURN
        
      END SUBROUTINE particles_compute_evgt
      
      
      
