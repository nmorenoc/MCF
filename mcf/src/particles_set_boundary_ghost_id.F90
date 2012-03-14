      SUBROUTINE particles_set_boundary_ghost_id(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_set_boundary_ghost_id
        !----------------------------------------------------
        !
        ! Purpose     : Setting IDs for ghost particles which
        !               are boundary particles,
        !               according to different boundary
        !               conditions.
        !               Record boundary ghosts particles
        !               index and Species ID.
        !                 
        !
        ! References  :
        !
        ! Remarks     : 
        !               For symmetry boundary, species ID
        !               is changed according to different
        !               colloid index, e.g., 
        !               two colloids approaching each other
        !               along x-axis, only one is inside
        !               domain with symmetry boundary,
        !               when the symmetry particle is outside
        !               of x_max, change species ID to 2.
        !
        !               For wall using symmetry,
        !               check each dimension at boundary and
        !               give negative species ID to ghost 
        !               particle which is outside of physical
        !               boundary for wall symmetry boundary,
        !               i.e., 
        !               x_min,x_max,y_min,y_max(,z_min,z_max)
        !               with -1,-2,-3,-4(,-5,-6), 
        !               respectively.
        !               
        !
        ! Revisions   : V0.3 21.05.2010, change colloid species
        !               ID for symmetry boundary condition.
        !
        !               V0.2 10.05.2010, negative species ID
        !               only make sense to wall_symmetry,
        !               remove the rest.
        !
        !               V0.1 07.12 2009, original version.
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
        ! Argumetns :
      	!----------------------------------------------------
        
        TYPE(Particles),INTENT(INOUT)   :: this
        INTEGER,INTENT(OUT)             :: stat_info
        
        !----------------------------------------------------
        ! Physics parameters :
     	!----------------------------------------------------
        
        INTEGER                         :: num_dim
        REAL(MK), DIMENSION(:), POINTER :: min_phys
        REAL(MK), DIMENSION(:), POINTER :: max_phys
        
        !----------------------------------------------------
        ! Boundary parameters :
     	!----------------------------------------------------
        
        INTEGER, DIMENSION(:), POINTER  :: bcdef
        TYPE(Boundary), POINTER         :: tboundary
        INTEGER                         :: num_sym
        INTEGER                         :: num_wall_sym
        INTEGER                         :: num_le

        !----------------------------------------------------
        ! Colloid paramters :
        !----------------------------------------------------
        
        INTEGER                         :: num_colloid
        
        !----------------------------------------------------
        ! Auxiliary parameters :
     	!----------------------------------------------------
        
        INTEGER                         :: i, j
        INTEGER                         :: num_ghost
        INTEGER                         :: stat_info_sub
        
        !----------------------------------------------------
        ! Initialization of variables.
      	!----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        NULLIFY(min_phys)
        NULLIFY(max_phys)

        NULLIFY(bcdef)
        NULLIFY(tboundary)        
        
        !----------------------------------------------------
        ! Physics parameters :
        !----------------------------------------------------
        
        num_dim = this%num_dim
        
        CALL physics_get_min_phys(this%phys,min_phys,stat_info_sub)
        CALL physics_get_max_phys(this%phys,max_phys,stat_info_sub)
        
        !----------------------------------------------------
        ! Boundary parameters :
        !----------------------------------------------------
       
        CALL physics_get_bcdef(this%phys,bcdef,stat_info_sub)
        CALL physics_get_boundary(this%phys,tboundary,stat_info_sub)
        
        num_sym        = &
             boundary_get_num_sym(tboundary,stat_info_sub)
        num_wall_sym   = &
             boundary_get_num_wall_sym(tboundary,stat_info_sub)
        num_le         = &
             boundary_get_num_le(tboundary,stat_info_sub)
        
        !----------------------------------------------------
        ! Colloid paramters :
        !----------------------------------------------------
        
        num_colloid    = & 
             physics_get_num_colloid(this%phys,stat_info_sub)
         
        !----------------------------------------------------
        ! Reset indices of 
        ! num_part_sym, num_part_wall_sym, 
        ! num_part_wall_solid_ghost, num_part_le to zero.
        !
        ! Prepare to count and record symmetry, wall using
        ! symmetry/mirror, solid wall and Lees-Edwards 
        ! boundary particles which ghosts of current 
        ! process(or processor).
        !
        ! Allocate memory to record them.
        ! num_ghost is the maximum it can be.
        !
        ! Note that these must be all ghosts particles.
        !----------------------------------------------------
        
        this%num_part_sym               = 0
        this%num_part_wall_sym          = 0
        this%num_part_le                = 0
        
        num_ghost = this%num_part_ghost
        
        IF ( num_sym > 0 ) THEN
           
           IF (ASSOCIATED(this%part_sym_list))THEN
              DEALLOCATE(this%part_sym_list)
           END IF
           
           ALLOCATE(this%part_sym_list(2,num_ghost))
           
        END IF
        
        IF ( num_wall_sym > 0 ) THEN
           
           IF (ASSOCIATED(this%part_wall_sym_list))THEN
              DEALLOCATE(this%part_wall_sym_list)
           END IF
           
           ALLOCATE(this%part_wall_sym_list(2,num_ghost))
           
        END IF
        
        IF ( num_le > 0 ) THEN
           
           IF (ASSOCIATED(this%part_le_list))THEN
              DEALLOCATE(this%part_le_list)
           END IF
           
           ALLOCATE(this%part_le_list(2,num_ghost))
           
        END IF
        
        
        !----------------------------------------------------
        ! Loop over each ghost particle.
        !----------------------------------------------------
        
        IF ( num_sym > 0 .OR. num_wall_sym > 0 .OR. &
             num_le > 0 ) THEN
           
           DO j = this%num_part_real + 1, this%num_part_all
              
              !----------------------------------------------
              ! Check each dimension at boundary and
              ! give species ID to ghost particle which is
              ! outside of boundary accordingly.
              !
              ! Negative species ID is given to each
              ! ghost layer,
              ! i.e., 
              ! x_min, x_max, y_min,y_max, z_min, z_max 
              ! with -1,-2,-3,-4,-5,-6, respectively.
              !----------------------------------------------
              
              DO i = 1, num_dim
                 
                 !-------------------------------------------
                 ! The one outside of lower boundary.
                 !-------------------------------------------
                 
                 IF( this%x(i,j) < min_phys(i) ) THEN
                    
                    SELECT CASE ( bcdef(2*i-1) ) 
                       
                    CASE ( ppm_param_bcdef_symmetry )
                       
                       this%num_part_sym = &
                            this%num_part_sym + 1
                       
                       this%part_sym_list(1, &
                            this%num_part_sym) = j
                       this%part_sym_list(2, &
                            this%num_part_sym) = &
                            this%id(this%sid_idx,j)
                       
                    CASE ( ppm_param_bcdef_wall_sym )
                       
                       this%id(this%sid_idx,j) = 1-2*i
                       
                       this%num_part_wall_sym = &
                            this%num_part_wall_sym + 1
                       
                       this%part_wall_sym_list(1, &
                            this%num_part_wall_sym) = j
                       this%part_wall_sym_list(2, &
                            this%num_part_wall_sym) = &
                            this%id(this%sid_idx,j)
                       
                    CASE ( ppm_param_bcdef_LE )
                       
                       !this%id(this%sid_idx,j) = 1-2*i
                       
                       this%num_part_le = &
                            this%num_part_le + 1
                       
                       this%part_le_list(1, &
                            this%num_part_le) = j
                       this%part_le_list(2, &
                            this%num_part_le) = &
                            this%id(this%sid_idx,j)
                       
                    END SELECT
                    
                    EXIT ! Found one dimension, quit
                    
                    !----------------------------------------
                    ! The one outside of upper boundary.
                    !----------------------------------------
                    
                 ELSE IF( this%x(i,j) >= max_phys(i) ) THEN
                    
                    SELECT CASE ( bcdef(2*i) ) 
                       
                    CASE ( ppm_param_bcdef_symmetry )
                       
                       this%num_part_sym = &
                            this%num_part_sym + 1
                       
                       this%part_sym_list(1, &
                            this%num_part_sym) = j
                       this%part_sym_list(2, &
                            this%num_part_sym) = &
                            this%id(this%sid_idx,j)
                       
                    CASE ( ppm_param_bcdef_wall_sym )
                       
                       this%id(this%sid_idx,j) = -2*i
                       
                       this%num_part_wall_sym = &
                            this%num_part_wall_sym + 1
                       
                       this%part_wall_sym_list(1, &
                            this%num_part_wall_sym) = j
                       this%part_wall_sym_list(2, &
                            this%num_part_wall_sym) = &
                            this%id(this%sid_idx,j)
                       
                    CASE ( ppm_param_bcdef_LE )
                       
                       !this%id(this%sid_idx,j) = -2*i
                       
                       this%num_part_le = &
                            this%num_part_le + 1
                       
                       this%part_le_list(1, &
                            this%num_part_le) = j
                       this%part_le_list(2, &
                            this%num_part_le) = &
                            this%id(this%sid_idx,j)
                       
                    END SELECT ! bcdef()
                    
                    EXIT ! Found one dimension, quit
                    
                 END IF ! min_phys, max_phys
                 
              END DO ! i = 1, num_dim
              
              !----------------------------------------------
              ! Assign  ghost particles with negative
              ! particle id to identify in the futural use.
              !
              ! Simple '-' would not do so,
              ! since one ghost particle could be a ghost 
              ! particle in last step also, '-' would 
              ! change its negative id to positive again.
              ! Therefore, use -ABS() function instead,
              ! which make sure that ghost id is always 
              ! negative.
              !----------------------------------------------
              
              !this%id(this%pid_idx,j) = &
              !     - ABS(this%id(this%pid_idx,j))
              
           END DO ! j = real+1, all
           
        END IF
        
9999    CONTINUE
        
        !----------------------------------------------------
        ! Release all dynamics memories.
        !----------------------------------------------------
        
        IF(ASSOCIATED(min_phys)) THEN
           DEALLOCATE(min_phys)
        END IF
        
        IF(ASSOCIATED(max_phys)) THEN
           DEALLOCATE(max_phys)
        END IF
        
        IF(ASSOCIATED(bcdef)) THEN
           DEALLOCATE(bcdef)
        END IF
        
        RETURN
        
      END SUBROUTINE particles_set_boundary_ghost_id
      
