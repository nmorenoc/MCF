      SUBROUTINE particles_init_global_assign_id(this,d_rank,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_init_global_assign_id
        !----------------------------------------------------
        !
        ! Purpose     : Create particles on certain lattice 
        !               globally and internally,
        !               x, y,(z), vx,vy,(vz), p_id,s_id.
        !
        ! 
        ! Routines    : 
        !
        !
        ! Remarks     : 
        !               We generate particles using different
        !               lattice.
        !               However, for colloid, a little bit 
        !               more attention has to be paid:
        !
        !               colloid boundary particles can be 
        !               on lattice, also can be on 
        !               the layers parallel to the surface,
        !               which further has different situations:
        !               1) each layer's particles have
        !                  equal distance, 
        !               2) or each layer has a fixed 
        !                  number of particles
        !                  (Zhu et al. 1999)
        !
        !               3) 4) 5)
        !
        ! References  :
        !
        ! Revisions   : 
        !               V0.1  08.03 2012, original version.
        !               (separated from particles_init_global_inter.F90)
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
        ! Modules :
        !----------------------------------------------------
        
        USE ppm_module_find_duplicates
        
        !----------------------------------------------------
        ! Arguments :
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)          :: this
        INTEGER, INTENT(IN)                     :: d_rank
        INTEGER, INTENT(OUT)	                :: stat_info
        
        !----------------------------------------------------
        ! Local variables start here :
        !
        ! num_dim      : number of dimension.
        ! dx           : initial distance between two particles;        
        ! tboundary    : pointer to a boundary object.
        ! num_wall_solid : number of solid wall boundaries.
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: num_dim,j
        REAL(MK), DIMENSION(:), POINTER         :: dx
        TYPE(Boundary), POINTER                 :: tboundary
        INTEGER                                 :: num_wall_solid
        
        !----------------------------------------------------
        ! Colloids parameters.
        !----------------------------------------------------
        
        INTEGER                                 :: num_colloid
        TYPE(Colloid), POINTER                  :: colloids
        INTEGER, DIMENSION(:),  ALLOCATABLE     :: coll_num_numerical_part
        INTEGER                                 :: coll_place
        REAL(MK), DIMENSION(:,:), POINTER       :: c_x
        INTEGER, DIMENSION(:), POINTER          :: c_sid
        
        !----------------------------------------------------
        ! Counters & Flags :
        !
        ! num_total : number of total particles.
        ! num       : counter for particles.
        ! num_extra :
        !             number of colloid boundary particles
        !             which are paralle layers to surface.
        ! num_more  : counter of more particles.
        !
        ! l_w   : particle is a solid wall boundary particle.
        !
        ! l_sur : particle is inside a colloid surface.
        ! l_out : particle is inside outter ring of a colloid,
        !         i.e., outside but close to the surface.
        ! l_in  : particle is inside inner ring of a colloid,
        !         i.e., inside but far from the surface which
        !         is computationally useless.
        !
        ! p_v   : particle's velocity.
        ! p_sid : particle's species ID.
        !----------------------------------------------------
        
        INTEGER                                 :: num_total
        INTEGER                                 :: num
        INTEGER                                 :: num_extra
        INTEGER                                 :: num_more
        LOGICAL                                 :: counted_wall_solid
        LOGICAL                                 :: counted_colloid
        LOGICAL                                 :: counted_ignore
        
        LOGICAL                                 :: l_w
        LOGICAL                                 :: l_sur
        LOGICAL                                 :: l_out
        LOGICAL                                 :: l_in
        INTEGER                                 :: p_sid
     
        !----------------------------------------------------
        ! Intermediate variables.
        !----------------------------------------------------
        
        REAL(MK), DIMENSION(:,:), ALLOCATABLE   :: tx
        INTEGER, DIMENSION(:,:),  ALLOCATABLE   :: tid 
        
        REAL(MK), DIMENSION(:,:), ALLOCATABLE   :: tc_x
        INTEGER, DIMENSION(:,:),  ALLOCATABLE   :: tc_id 
    
        !----------------------------------------------------
        ! Check parameters.
        !----------------------------------------------------        
        
        INTEGER , DIMENSION(:,:), POINTER       :: ide
        INTEGER                                 :: nid
        
        
#ifdef __DEBUG
        REAL(MK), DIMENSION(:,:), ALLOCATABLE   :: ppp
#endif 
        !----------------------------------------------------
        ! Check if this routine is called by root process.
        !----------------------------------------------------

        IF ( d_rank /= 0 ) THEN
           PRINT *,&
                "particles_init_global_assign_id : ", &
                "can only be called by root process !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
    	! Initialization of variables.
    	!----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        num_dim = this%num_dim
        
        NULLIFY(dx)
        NULLIFY(tboundary)
        
        !----------------------------------------------------
        ! Initialization of variables of colloids.
        !----------------------------------------------------
        
        NULLIFY(colloids)
        
        NULLIFY(c_x)
        NULLIFY(c_sid)
        
        
        NULLIFY(ide)
        
        !----------------------------------------------------
        ! Get physics variables,
        ! dimension and which lattice.
        ! number of solid walls;
        ! number of colloids.
        !----------------------------------------------------
        
        num_dim        = this%num_dim
        CALL physics_get_dx(this%phys,dx,stat_info_sub)
        CALL physics_get_boundary(this%phys,tboundary,stat_info_sub)
        num_wall_solid = &
             boundary_get_num_wall_solid(tboundary,stat_info_sub)
        num_colloid    = &
             physics_get_num_colloid(this%phys,stat_info_sub)
        
#ifdef __DEBUG

#ifdef __DEBUG_INIT
        !----------------------------------------------------
        ! For debug purpose, may be written into files.
        !----------------------------------------------------
        
        num = this%num_part_real
        
        ALLOCATE(ppp(num_dim+1,num))
        
        ppp(1:num_dim,1:num) = this%x(1:num_dim,1:num)
        ppp(num_dim+1,1:num) = this%id(2,1:num)
        
        CALL debug_write_output(global_debug,d_rank,&
             "particles_init_global_inter ", &
             "ppp_real",0,ppp,1,this%num_part_real,stat_info_sub)
        
#endif
        
#endif
        !----------------------------------------------------
        ! If there is colloid
        !----------------------------------------------------
        
        IF ( num_colloid > 0 ) THEN

#ifdef __COLLOID_ON_LATTICE
        
           !-------------------------------------------------
           ! Adjust colloid position to make its center
           ! on lattice, therefore the generated particles
           ! can be symmetrical according to the center.
           !-------------------------------------------------
           
           CALL particles_set_colloid_on_lattice(this,stat_info_sub)
           
           IF ( stat_info_sub /=0 ) THEN
              
              PRINT * , "particles_init_global_assign_id : ", &
                   "calling particles_set_colloid_on_lattice failed !"
              stat_info = -1
              GOTO 9999
           END IF
           
#endif
           
           !-------------------------------------------------
           ! Get colloid parameters
           ! 1: how the boundary particles are placed
           ! 2: compute locations of colloid images 
           !    according to boundary conditions, 
           !    such as periodic and Lees-Edwards.
           !-------------------------------------------------
           
           CALL physics_get_colloid(this%phys, &
                colloids,stat_info_sub)
           coll_place= &
                colloid_get_place(colloids,stat_info_sub)
           
           CALL colloid_compute_image(colloids,stat_info_sub)
           
        END IF
        
        !----------------------------------------------------
    	! Allocate temporary memory (large enough) of
        ! fluid particles to keep particles.
        !
	! tx   : position
       	! tid  : particle ID,  species ID
        ! this%num_part_real: 
        !       total number of particles generated.
    	!----------------------------------------------------
        
        num_total = this%num_part_real
        
        ALLOCATE(tx(num_dim,num_total), STAT=stat_info_sub)
        ALLOCATE(tid(this%num_id,num_total),STAT=stat_info_sub)
        
        IF( stat_info_sub /= 0 ) THEN
           PRINT *, &
                "particles_init_global_assign_id : ", &
                "Allocating memory for variables has problem !"
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Reset all counters to zero,
        ! prepare for counting 
        ! number of fluid particles;
        ! number of wall solid boundary particles;
        ! number of colloid boundary particles.
        !----------------------------------------------------
        
        this%num_part_fluid = 0
        this%num_part_wall_solid = 0
        
        IF ( num_colloid > 0 ) THEN
           
           ALLOCATE(coll_num_numerical_part(1:num_colloid))
           coll_num_numerical_part(1:num_colloid) = 0
           
        END IF
        
        !----------------------------------------------------
        ! Loop over all generated particles to decide 
        ! species ID, i.e., which type of particle it is.
        !----------------------------------------------------
        
        num = 0

        DO j = 1, num_total
           
           !-------------------------------------------------
           ! Before checking, this particle is assumed solvent
           ! and its velocity is zero.
           !-------------------------------------------------

           p_sid              = 0
           counted_wall_solid = .FALSE.
           counted_colloid    = .FALSE.
           counted_ignore     = .FALSE.           
           
           !-------------------------------------------------
           ! Check if x is a wall(solid) boundary particle.
           !-------------------------------------------------
           
           IF ( num_wall_solid > 0 ) THEN
              
              l_w = .FALSE.
              
              CALL boundary_check_wall_solid_particle(tboundary, &
                   this%x(1:num_dim,j),l_w,p_sid,stat_info_sub)
              
              IF ( stat_info_sub /= 0 ) THEN
                 PRINT *, "paricles_init_global_inter : ", &
                      "Checking wall solid particle failed !"
                 stat_info = -1
                 GOTO 9999
              END IF
              
              IF ( l_w ) THEN
                 
                 counted_wall_solid = .TRUE.
                 
              END IF
              
           END IF ! num_wall_solid           
           
           !-------------------------------------------------
           ! If sx is not solid wall boundary particle,
           ! check if it is inside any colloid geometry.
           !-------------------------------------------------
           
           IF ( .NOT. counted_wall_solid .AND. &
                num_colloid > 0 ) THEN
              
              !----------------------------------------------
              ! loop over all colloids 
              !----------------------------------------------
              
              l_sur = .FALSE.
              l_out = .FALSE.
              l_in  = .FALSE.
              
              CALL colloid_check_boundary_particle(colloids, &
                   this%x(1:num_dim,j),l_sur,l_out,l_in, p_sid,&
                   stat_info_sub)
              
              IF ( l_sur .AND. &
                   coll_place == mcf_colloid_place_lattice ) THEN
                 
                 !-------------------------------------------
                 ! particle j is inside surface of a colloid
                 ! and the placement of it is on lattice.
                 !-------------------------------------------
                 
                 IF( .NOT. l_in ) THEN
                    
                    counted_colloid = .TRUE.
                 
                 ELSE
                    
                    counted_ignore = .TRUE.
                    
                 END IF ! .NOT. l_in
                 
              ELSE IF( l_out .AND. &
                   coll_place /= mcf_colloid_place_lattice ) THEN
                 
                 !-------------------------------------------
                 ! If j is outside of a colloid, but close
                 ! enough to the surface, it should not be 
                 ! created as fluid particle when colloid
                 ! boundary particles are required to be 
                 ! placed on paralle surfaces, since sx would
                 ! be too close to the boundary particles
                 ! on the surface.
                 !-------------------------------------------
                 
                 counted_ignore = .TRUE.
                 
              END IF ! l_sur
              
           END IF ! NOT counted_wall_solid AND num_colloid > 0
           
#if __PARTICLES_NO_SOLVENT
           !-------------------------------------------------
           ! ignore solvent particles, for special usage only.
           !-------------------------------------------------
           
           IF ( (.NOT. counted_wall_solid ) .AND. &
                (.NOT. counted_colloid ) ) THEN
              
              counted_ignore = .TRUE.
              
           END IF
#endif
           
           !-------------------------------------------------
           ! If the particle is not ignored,
           ! save the particle's information.
           !-------------------------------------------------
           
           IF ( .NOT. counted_ignore ) THEN
              
              num = num + 1
              tx(1:num_dim,num)     = this%x(1:num_dim,j)
              tid(this%pid_idx,num) = num
              tid(this%sid_idx,num) = p_sid
              
              !----------------------------------------------
              ! Increase counters of different species.
              !----------------------------------------------
              
              IF ( counted_wall_solid )  THEN
                 
                 this%num_part_wall_solid = &
                      this%num_part_wall_solid + 1
                 
              ELSE IF ( counted_colloid )  THEN
                 
                 coll_num_numerical_part(p_sid) = &
                      coll_num_numerical_part(p_sid) + 1
                 
              ELSE
                 
                 this%num_part_fluid = &
                      this%num_part_fluid + 1
                 
              END IF ! counted
              
           END IF ! .NOT. counted_ignore
           
        END DO ! j = 1, num_total
        
        !----------------------------------------------------
        ! Record number of solid wall boundary particles.
        !----------------------------------------------------
        
        CALL boundary_set_num_part_wall_solid(tboundary, &
             this%num_part_wall_solid,stat_info_sub)
        
        !----------------------------------------------------
        ! Check if we need to create colloid boundary 
        ! particles in other special ways,
        ! e.g., on parallel surfaces of colloid.
        !
        ! If needed, we call colloid_create_boundary_particle*
        ! to create them.
        !----------------------------------------------------
        
        num_extra = 0
        
        IF ( num_colloid > 0 .AND. &
             coll_place /=  mcf_colloid_place_lattice ) THEN
           
           IF ( num_dim == 2 ) THEN
              
              CALL colloid_create_boundary_particle_2D(colloids,&
                   dx(1:num_dim),c_x,c_sid,stat_info_sub)
              
           ELSE IF ( num_dim == 3 ) THEN
              
              CALL colloid_create_boundary_particle_3D(colloids,&
                   dx(1:num_dim),c_x,c_sid,stat_info_sub)
              
           END IF

           IF ( stat_info_sub /= 0 ) THEN

              PRINT *, "particles_init_global_assign_id: ", &
                   "colloid_create_boundary_particle has problem !"
              stat_info = -1
              GOTO 9999
              
           END IF
           
           num_extra = SIZE(c_sid,1)
           
        END IF ! num_colloid > 0 AND coll_place /= mcf_colloid_place_lattice

        !----------------------------------------------------
        ! If there is extra boundary particles created
        ! in other special ways,
        ! loop over all generated particles to decide 
        ! which boundary particle is really needed and
        ! which one is computationally useless.
        !----------------------------------------------------
        
        num_more = 0
        
        IF ( num_extra > 0 ) THEN
           
           !----------------------------------------------------
           ! Allocate temporary memory (large enough) of
           ! to keep boundary particles.
           !
           ! tc_x   : position
           ! tc_id  : particle ID,  species ID
           !----------------------------------------------------
           
           ALLOCATE(tc_x(num_dim,num_extra), STAT=stat_info_sub)
           ALLOCATE(tc_id(this%num_id,num_extra),STAT=stat_info_sub)
           
           IF( stat_info_sub /= 0 ) THEN
              PRINT *, &
                   "particles_init_global_assig_id: ", &
                   "Allocating memory for variables has problem !"
              stat_info = -1
              GOTO 9999
           END IF
           
           DO j = 1, num_extra
              
              counted_ignore    = .FALSE.
              
              !----------------------------------------------
              ! Check if it is inside any colloid geometry.
              ! loop over all colloids 
              !----------------------------------------------
#if 0
              l_in  = .FALSE.
              
              CALL colloid_check_boundary_particle(colloids, &
                   c_x(1:num_dim,j),l_sur,l_out,l_in, &
                   c_sid(j),stat_info_sub)
           
              
              IF ( l_in ) THEN
                 
                 counted_ignore = .TRUE.
                 
              END IF
#endif 
              !----------------------------------------------
              ! If it is compuationally usefull,
              ! save the boundary particle's information.
              !----------------------------------------------
              
              IF ( .NOT. counted_ignore ) THEN
                 
                 num_more = num_more + 1
                 tc_x(1:num_dim,num_more)      = c_x(1:num_dim,j)
                 tc_id(this%pid_idx,num_more)  = num + num_more
                 tc_id(this%sid_idx,num_more)  = c_sid(j)
                 
                 coll_num_numerical_part(c_sid(j)) = &
                      coll_num_numerical_part(c_sid(j)) + 1
                 
              END IF ! .NOT. counted_ignore
              
           END DO ! j = 1, num_total

        END IF ! num_extra > 0
        
        
        !----------------------------------------------------
        ! Copy the data from temporary varaiable 
        ! to particles members.
	!----------------------------------------------------
        
        num_total = num + num_more
        
        IF (ASSOCIATED(this%x)) THEN
           DEALLOCATE(this%x,STAT=stat_info_sub)
        END IF
        IF (ASSOCIATED(this%v)) THEN
           DEALLOCATE(this%v,STAT=stat_info_sub)
        END IF
        IF (ASSOCIATED(this%id)) THEN
           DEALLOCATE(this%id,STAT=stat_info_sub)
        END IF
        
        ALLOCATE(this%x(num_dim,num_total), STAT=stat_info_sub)
        ALLOCATE(this%v(num_dim,num_total), STAT=stat_info_sub)
        ALLOCATE(this%id(this%num_id,num_total),STAT=stat_info_sub)
        
        this%x(1:num_dim,1:num) = tx(1:num_dim,1:num)
        this%v(1:num_dim,1:num) = 0.0_MK
        this%id(1:this%num_id,1:num) = &
             tid(1:this%num_id,1:num)
        
        !----------------------------------------------------
        ! If there is extra particles, i.e.,
        ! created on layers parallel to colloid surface.
        !----------------------------------------------------
        
        IF ( num_more > 0 ) THEN
           
           this%x(1:num_dim,num+1:num_total) = &
                tc_x(1:num_dim,1:num_more)
           this%v(1:num_dim,num+1:num_total) = 0.0_MK
           this%id(1:this%num_id,num+1:num_total) = &
                tc_id(1:this%num_id, 1:num_more)
           
        END IF ! num_more > 0
        
        !----------------------------------------------------
        ! Count the number of numerical particles which 
        ! consititute colloids.
        !----------------------------------------------------
        
        this%num_part_colloid = 0
        
        IF ( num_colloid > 0 ) THEN
           
           DO j = 1, num_colloid
              this%num_part_colloid = &
                   this%num_part_colloid + coll_num_numerical_part(j)
           END DO
           
           CALL colloid_set_num_numerical_part(colloids,&
                coll_num_numerical_part(:),stat_info_sub)
        
        END IF
        
        !----------------------------------------------------
        ! Set number of real, all and ghost particles.
        !----------------------------------------------------
        
        this%num_part_real  = num_total
        this%num_part_all   = this%num_part_real
        this%num_part_ghost = 0
        
#ifdef __DEBUG

#ifdef __DEBUG_INIT
        !----------------------------------------------------
        ! For debug purpose, could be written into files.
        !----------------------------------------------------
        
        num = this%num_part_real

        IF(ALLOCATED(ppp)) THEN
           DEALLOCATE(ppp)
        END IF
        
        ALLOCATE(ppp(num_dim+1,num))
        
        ppp(1:num_dim,1:num) = this%x(1:num_dim,1:num)
        ppp(num_dim+1,1:num) = this%id(2,1:num)
        
        CALL debug_write_output(global_debug,d_rank,&
             "particles_init_global_assign_id: ", &
             "ppp_real",1,ppp,1,this%num_part_real,stat_info_sub)

#endif
        
#endif
        
        !----------------------------------------------------
    	! Check if each particles is unique.
    	!----------------------------------------------------
        
        nid = 0
        
        CALL ppm_find_duplicates(this%x,num_dim,&
             this%num_part_real,nid,ide,stat_info_sub)
        
        IF ( stat_info_sub /= 0 ) THEN
           PRINT *, &
                "particles_init_global_assign_id : ", &
                "Error by checking duiplicates !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF ( nid /= 0) THEN
           PRINT *, &
                "particles_init_global_assign_id : ", &
                "Found collocating particles !"
           stat_info = -1
           GOTO 9999
        END IF
        
        
9999    CONTINUE
        
        IF(ASSOCIATED(dx)) THEN
           DEALLOCATE(dx)
        END IF
          
        IF(ASSOCIATED(c_x)) THEN
           DEALLOCATE(c_x)
        END IF
        
        IF(ASSOCIATED(c_sid)) THEN
           DEALLOCATE(c_sid)
        END IF
        
        IF ( ASSOCIATED(ide) ) THEN
           DEALLOCATE(ide)
        END IF
        
        RETURN
        
      END SUBROUTINE particles_init_global_assign_id
  
