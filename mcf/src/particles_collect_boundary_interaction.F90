      SUBROUTINE particles_collect_boundary_interaction(this,&
           drag,drag_p,drag_v,drag_r,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_collect_boundary_interaction
        !----------------------------------------------------
        !
        ! Purpose     :  Summing up interaction(force, etc.)
        !                exterted on parts of walls on this
        !                process.
        !                  
        !  References  :                   
        !
        !  Remarks     : Wall boundary is also assumed to be 
        !                a shear boundary.
        !
        !  Revisions   : V0.2 07.12 2009, including general
        !                boundary conditions.
        !
        !                V0.1 23.10 2009, original version.
        !----------------------------------------------------
        !  Author      :  Xin Bian 
        !  Contact     :  xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------

        !----------------------------------------------------
        ! Arguments :
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(IN)     :: this
        REAL(MK),DIMENSION(:,:)         :: drag
        REAL(MK),DIMENSION(:,:),OPTIONAL:: drag_p
        REAL(MK),DIMENSION(:,:),OPTIONAL:: drag_v
        REAL(MK),DIMENSION(:,:),OPTIONAL:: drag_r
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Local variables start here :
        !----------------------------------------------------     
        
        INTEGER                         :: stat_info_sub
        LOGICAL                         :: symmetry
        LOGICAL                         :: Brownian
        INTEGER                         :: num_dim, j
        TYPE(Boundary), POINTER         :: tboundary
        INTEGER                         :: num_wall_sym
        INTEGER                         :: num_wall_solid
        INTEGER                         :: num_le
        INTEGER                         :: ip,sid
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        
        !----------------------------------------------------
        ! Control parameters :
        !----------------------------------------------------

        symmetry = &
             control_get_symmetry(this%ctrl,stat_info_sub)
        Brownian = &
             control_get_Brownian(this%ctrl,stat_info_sub)
      
        !----------------------------------------------------
        ! Physics parameters :
        !----------------------------------------------------
        
        num_dim = this%num_dim
        
        IF ( SIZE(drag,1) /= num_dim ) THEN
           PRINT *, "particles_collect_boundary_interaction : ", &
                "drag dimension does not match ! "
           stat_info = -1
           GOTO 9999
        END IF
        
        IF ( SIZE(drag,2) /= num_dim*2 ) THEN
           PRINT *, "particles_collect_boundary_interaction : ", &
                "drag number does not match ! "
           stat_info = -1
           GOTO 9999
        END IF
        
        drag(:,:) = 0.0_MK
        
#ifdef __WALL_FORCE_SEPARATE
        drag_p(:,:) = 0.0_MK
        drag_v(:,:) = 0.0_MK
        drag_r(:,:) = 0.0_MK        
#endif        
        
        NULLIFY(tboundary)
        
        CALL physics_get_boundary(this%phys,tboundary,&
             stat_info_sub)
        
        num_wall_sym = &
             boundary_get_num_wall_sym(tboundary,stat_info_sub)
        num_wall_solid = &
             boundary_get_num_wall_solid(tboundary,stat_info_sub)
        num_le         = &
             boundary_get_num_le(tboundary,stat_info_sub)
        
        
        !----------------------------------------------------
        ! Sum up each shear boundary particle's contribution.
        !----------------------------------------------------
        
        IF ( .NOT. symmetry ) THEN
           
           IF ( num_wall_sym  > 0 ) THEN
              
              DO j = 1, this%num_part_wall_sym
                 
                 ip  = this%part_wall_sym_list(1,j)
                 sid =  ABS(this%part_wall_sym_list(2,j))
                 
                 drag(1:num_dim,sid) = &
                      drag(1:num_dim,sid) + &
                      this%f(1:num_dim,ip) * this%m(ip)
                 
#ifdef __WALL_FORCE_SEPARATE
                 drag_p(1:num_dim,sid) = &
                      drag_p(1:num_dim,sid) + &
                      this%fp(1:num_dim,ip) * this%m(ip)
                 drag_v(1:num_dim,sid) = &
                      drag_v(1:num_dim,sid) + &
                      this%fv(1:num_dim,ip) * this%m(ip)
                 IF ( Brownian ) THEN
                    drag_r(1:num_dim,sid) = &
                         drag_r(1:num_dim,sid) + &
                         this%fr(1:num_dim,ip) * this%m(ip)
                 END IF
#endif
                 
              END DO ! j = 1, num_part_wall_sym
              
           END IF ! num_wall_sym

        END IF
        
        
        IF ( num_wall_solid > 0 ) THEN
           
           DO j = 1, this%num_part_wall_solid_real
              
              ip  = this%part_wall_solid_real_list(1,j)
              sid = ABS(this%part_wall_solid_real_list(2,j))
              
              drag(1:num_dim,sid) = &
                   drag(1:num_dim,sid) + &
                   this%f(1:num_dim,ip) * this%m(ip)
              
#ifdef __WALL_FORCE_SEPARATE
              drag_p(1:num_dim,sid) = &
                   drag_p(1:num_dim,sid) + &
                   this%fp(1:num_dim,ip) * this%m(ip)
              drag_v(1:num_dim,sid) = &
                   drag_v(1:num_dim,sid) + &
                   this%fv(1:num_dim,ip) * this%m(ip)
              IF ( Brownian ) THEN
                 drag_r(1:num_dim,sid) = &
                      drag_r(1:num_dim,sid) + &
                      this%fr(1:num_dim,ip) * this%m(ip)
              END IF
#endif

           END DO ! j = 1, num_part_wall_solid_real
           
        END IF ! num_wall_solid
        
        
9999    CONTINUE
        
        
        RETURN
        
      END SUBROUTINE particles_collect_boundary_interaction
      
