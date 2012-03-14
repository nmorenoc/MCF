      SUBROUTINE  particles_set_colloid_on_lattice(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_set_colloid_on_lattice
        !----------------------------------------------------
        !
        ! Purpose     : Put colloid on lattice, to have its
        !               shape symmetrically represented by
        !               boundary particle.
        !
        ! Routines    :
        !
        ! Remarks     : 1) No images of colloid due to 
        !               different boundaries are considered yet.
        !               2) Only spherical colloid.
        !               3) Only for lattice type 1.
        !
        ! References  :
        !
        ! Revisions   : V0.1 9.12.2010, original version.
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
        ! Arguments :
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)          :: this
        INTEGER,INTENT(OUT)	                :: stat_info
        
        !----------------------------------------------------
    	! Local variables start here :
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        
        INTEGER                                 :: dim, i, j
        INTEGER                                 :: lattice_type
        REAL(MK), DIMENSION(:), POINTER         :: min_phys
        REAL(MK), DIMENSION(:), POINTER         :: max_phys
        
        INTEGER                                 :: num_colloid
        TYPE(Colloid), POINTER                  :: colloids
        REAL(MK), DIMENSION(:,:), POINTER       :: coll_x
        REAL(MK), DIMENSION(:,:), POINTER       :: p_x
        REAL(MK), DIMENSION(:), POINTER         :: dist_x
        REAL(MK), DIMENSION(3)                  :: ax,bx
        REAL(MK)                                :: dx
        
        !----------------------------------------------------
    	! Initialization of variables.
    	!----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        NULLIFY(min_phys)
        NULLIFY(max_phys)
        NULLIFY(colloids)
        NULLIFY(coll_x)
        NULLIFY(p_x)
        NULLIFY(dist_x)

        
        !----------------------------------------------------
        ! Get physics including colloid parameters.
        !----------------------------------------------------
        dim = this%num_dim
        lattice_type = &
             physics_get_lattice_type(this%phys,stat_info_sub)        
        CALL physics_get_min_phys(this%phys, &
             min_phys,stat_info_sub)
        CALL physics_get_max_phys(this%phys, &
             max_phys,stat_info_sub)
        num_colloid = &
             physics_get_num_colloid(this%phys,stat_info_sub)
        
        IF ( lattice_type ==1 .AND. &
             num_colloid > 0 ) THEN
           
           CALL physics_get_colloid(this%phys, &
                colloids,stat_info_sub)
           
           CALL colloid_get_x(colloids,coll_x,stat_info_sub)
           
           ALLOCATE(p_x(dim,num_colloid))
           ALLOCATE(dist_x(num_colloid))
           
           dist_x(:) = (max_phys(1) - min_phys(1))**2
           
           !-------------------------------------------------
           ! Loop over all particles and find each nearest
           ! one to each colloid center.
           !-------------------------------------------------
           
           DO j = 1, this%num_part_real
              
              ax(1:dim) = this%x(1:dim,j)
              
              DO i = 1, num_colloid
                 
                 bx(1:dim) = ax(1:dim) - coll_x(1:dim,i)
                 
                 dx= DOT_PRODUCT(bx(1:dim), bx(1:dim))
                 
                 IF ( dx < dist_x(i) ) THEN
                    dist_x(i)    = dx
                    p_x(1:dim,i) = ax(1:dim)
                 END IF
                 
              END DO ! i
              
           END DO ! j
           
           CALL colloid_set_x(colloids,p_x(1:dim,&
                1:num_colloid),stat_info_sub)
           
        END IF ! num_colloid > 0
        
        !----------------------------------------------------
        ! Return.
        !----------------------------------------------------
        
9999	CONTINUE	
        
        IF(ASSOCIATED(min_phys)) THEN
           DEALLOCATE(min_phys) 
        END IF
        
        IF(ASSOCIATED(max_phys)) THEN
           DEALLOCATE(max_phys) 
        END IF
        
        IF(ASSOCIATED(coll_x)) THEN
           DEALLOCATE(coll_x)
        END IF

        IF(ASSOCIATED(p_x)) THEN
           DEALLOCATE(p_x)
        END IF

        IF(ASSOCIATED(dist_x)) THEN
           DEALLOCATE(dist_x)
        END IF
        
        RETURN
        
      END SUBROUTINE particles_set_colloid_on_lattice
      
      
