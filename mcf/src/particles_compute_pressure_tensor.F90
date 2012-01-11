      SUBROUTINE particles_compute_pressure_tensor(this,num,stat_info)
        !----------------------------------------------------
        ! Subroutine :  particles_compute_pressure_tensor
        !----------------------------------------------------
        !
        ! Purpose    :  Compute the pressure tensor
        !               of Non-Newtonian
        !               oldroyd-B model.
        !               
        !
        ! Reference  :  Vazquez-Quesada, Ellero, Espanol
        !               Phyical Review E 79. 056707, 2009.
        !
        ! Remark     :
        !
        ! Revision   :  V0.1  25.08.2010
        !               The osmotic pressure is added in the
        !               computation of the pressure tensor.
        !               Also was deleted the dim2 variable
        !               because it was unnecesary.
        !               (Adolfo)
        !               
        !               V0.1  16.04.2010, original version.
        !               Originally was done in 
        !               particle-particle interaction,
        !               which saved memory for a matrix of
        !               each particle, but
        !               noted by Adolf Vazquez-Quesada that
        !               it wasn't efficient, since pressure
        !               tensure of particle i is not related
        !               to particle j at all.
        !               Considering memory and speed of the
        !               code, we decided to put it seperatedly.
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
        ! num            : number of particles needed updated,
        !                  i.e. first num particles in this%x 
        !                  are operated.
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
        INTEGER                                 :: i,j,k
        REAL(MK)                                :: n_p, kt_p, G

        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        !----------------------------------------------------
        ! Calculation only for real particles.
        !----------------------------------------------------
        
        IF( num > this%num_part_all) THEN
           PRINT *, "particles_compute_pt : ", &
                "num > num_part_all, wrong !"
           stat_info = -1
           GOTO 9999      
        END IF
        

        dim = &
             physics_get_num_dim(this%phys,stat_info_sub)
!        dim2 = dim**2
        
        n_p   = physics_get_n_p(this%phys,stat_info_sub)
        kt_p  = physics_get_kt_p(this%phys,stat_info_sub)
        G     = n_p * kt_p
        
        !----------------------------------------------------
        ! Allocate memory for pt.
        !----------------------------------------------------
        
        IF(ASSOCIATED(this%pt)) THEN
           DEALLOCATE(this%pt)
        ENd IF
        
!!$        ALLOCATE(this%pt(dim2,dim2,num))
        ALLOCATE(this%pt(dim,dim,num))
        
        this%pt(1:dim,1:dim,1:num) = 0.0_MK
        
        !----------------------------------------------------
        ! Compute pressure tensor.
        !----------------------------------------------------
        
        
        DO k = 1, num ! particle index
           
           DO j =1, dim ! --- row direction
              
              Do i=1, dim ! | column direction
                 
                 this%pt(i,j,k) = -G*this%ct(i+dim*(j-1),k)
                 
                 IF ( i==j )  THEN
                    
                    !-- The osmotic pressure G is also added --
                    this%pt(i,j,k) = this%pt(i,j,k) + this%p(k) + G

                 END IF
                 
              END DO
              
           END DO
           
        END DO
        
        
9999    CONTINUE      
        
        RETURN
        
      END SUBROUTINE particles_compute_pressure_tensor
      
      
