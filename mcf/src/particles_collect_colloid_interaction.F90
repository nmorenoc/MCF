      SUBROUTINE particles_collect_colloid_interaction(this,&
           drag,torque,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_collect_colloid_interaction
        !----------------------------------------------------
        !
        ! Purpose     : Summerizing up force/torque contributed
        !               from current process on colloids.
        !               
        ! Routines    :   
        !
        ! References  :  
        !
        ! Remarks     : Since the force exerted on each single
        !               numerical particle of a colloid object
        !               has already been calculated as pair-wise
        !               interaction, here we accumulate 
        !               the total force/torque on parts 
        !               of colloids on each process.
      
        !
        ! Revisions   : V0.1 23.10 2009, original version.
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
        ! Arguments:
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(IN)             :: this
        REAL(MK),DIMENSION(:,:),INTENT(INOUT)   :: drag
        REAL(MK),DIMENSION(:,:),INTENT(INOUT)   :: torque
        INTEGER, INTENT(OUT)                    :: stat_info

        !----------------------------------------------------
        ! Local parameters: 
        !----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        INTEGER                         :: dim, j
        INTEGER                         :: num
        TYPE(Colloid), POINTER          :: colloids
        INTEGER                         :: ip,sid
        REAL(MK), DIMENSION(3)          :: F
        REAL(MK), DIMENSION(3)          :: rx
        REAL(MK), DIMENSION(3)          :: T

        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info = 0
        
        NULLIFY(colloids)
        
        !----------------------------------------------------
        ! Physics parameters, including colloid parameters
        !----------------------------------------------------
        
        dim = this%num_dim
        
        num = &
             physics_get_num_colloid(this%phys,stat_info_sub)
        
        IF ( num > 0 ) THEN
           
           CALL physics_get_colloid(this%phys,colloids,stat_info_sub)
           
        END IF
        
        
        !----------------------------------------------------
        ! Reset drag and torque.
        !----------------------------------------------------
        
        drag(1:dim,1:num) = 0.0_MK
        torque(1:3,1:num) = 0.0_MK
        
        
        !----------------------------------------------------
        ! Sum up each boundary particle's contribution.
        !----------------------------------------------------
        
        DO j = 1, this%num_part_colloid
           
           ip  = this%part_colloid_list(1,j)
           sid = this%part_colloid_list(2,j)
           
           !-------------------------------------------------
           ! Sum up force of colloidal boundary particles.
           !-------------------------------------------------
           
           F(1:3) = 0.0_MK
           
           F(1:dim) = this%f(1:dim,ip) * this%m(ip)
           
           drag(1:dim,sid) = drag(1:dim,sid) + F(1:dim)
           
           !-------------------------------------------------
           ! Return relative position of boundary particle
           ! to its colloid center.
           !-------------------------------------------------
           
           rx(1:3) = 0.0_MK
           CALL colloid_in_relative_position(colloids,&
                this%x(1:dim,ip),sid,rx(1:dim),stat_info_sub)
           
           IF ( stat_info_sub /= 0 ) THEN
              PRINT *, __FILE__, ":", __LINE__
              stat_info = -1
              GOTO 9999
           END IF
           !-------------------------------------------------
           ! Sum up torque of colloidal boundary particles
           !-------------------------------------------------
           
           T(1:3) = 0.0_MK
           
           CALL tool_cross_product(this%tool,&
                rx(1:3),F(1:3),&
                T(1:3),stat_info_sub)
           
           torque(1:3,sid) = &
                torque(1:3,sid) + T(1:3)
           
           
        END DO ! j = 1, num_part_colloid
        
        
9999    CONTINUE
        
        
        RETURN
        
      END SUBROUTINE particles_collect_colloid_interaction
        
      
