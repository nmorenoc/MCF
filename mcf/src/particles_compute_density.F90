      SUBROUTINE particles_compute_density(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particls_compute_density
        !----------------------------------------------------
        !
        ! Purpose     : Computing 2D/3D density using 
        !               symmetric or non-symmetric
        !               inter-process communication.
        !	 	      	 
        !                
        ! Reference   :
        !
        ! Remark      :
        !               
        !
        ! Revisions   : V0.7 07.12 2009, remove count
        !               of boundary particles to
        !               to particles_set_boudary_ID().
        !
        !               V0.6 09.10 2009, merge 2D,3D
        !               symmetry, non-symmetry together.
        !
        !               V0.5 29.09 2009, boundary particles'
        !               (including colloid and wall boundary
        !               particles) density are changed to be
        !               constant.
        !
        !               V0.4 23.07 2009, boundary partices' 
        !               density also evlove as fluid particles,
        !               same as Morris et al. 1997 and
        !               Zhu et al. 1999.
        ! 
        !               V0.3 13.07 2009, kinetic energy is not
        !               constant for simple pressure force.
        !               check again the work flow is correct 
        !               and supply with more comments.
        !
        !               V0.2 08.07 2009, 
        !               check again the work flow is correct 
        !               and supply with more comments.
        !
        !               V0.1 06.02.2009, original version.
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
        
        USE ppm_module_neighlist
        
	!----------------------------------------------------
      	! Arguments :
      	!----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)  :: this
        INTEGER, INTENT(OUT)		:: stat_info	
        
      	!----------------------------------------------------
      	! Local variables start here :
        !----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        
        !----------------------------------------------------
        ! Control parameters :
        !
        ! rhs_density_type : formulation of density.
        ! symmetry         : indicate if we use symmetric
        !                    inter-communication or not.
        !----------------------------------------------------
        
        INTEGER                         :: rhs_density_type
        LOGICAL                         :: symmetry
        
        !----------------------------------------------------
        ! Physics parameters :
        !
        ! num_dim      : physics dimension.
        ! cut_off      : compact support domain length.
        ! cut_off2     : cut_off * cut_off.
        ! init_density : initial density.
        !----------------------------------------------------
        
        INTEGER                         :: num_dim
        REAL(MK)                        :: cut_off
        REAL(MK)                        :: cut_off2
        REAL(MK)                        :: init_density
        
        !----------------------------------------------------
        ! Boundary parameters :
        !
        ! num_wall_solid : number of solid wall boundaries.
        !----------------------------------------------------
        
        INTEGER, DIMENSION(:), POINTER  :: bcdef
        TYPE(Boundary), POINTER         :: tboundary
        INTEGER                         :: wall_rho_type
        INTEGER                         :: num_wall_solid
        INTEGER                         :: num_part_wall_solid
        
        !----------------------------------------------------
        ! Colloid parameters :
        !
        ! num_colloid      : number of colloidal particle.
        ! num_part_colloid : number of colloid boundary 
        !                    particles.       
        !----------------------------------------------------
        
        INTEGER                         :: num_colloid
        TYPE(Colloid), POINTER          :: colloids
        INTEGER                         :: coll_rho_type
        INTEGER                         :: num_part_colloid
        
        !----------------------------------------------------
        ! Cell list, indices / counters :
        !
        ! num_sub          : number of sumdomains on 
        !                    current process
        ! num_cell_dim_sub : number of cells in each direction 
        !                    of each subdomain.
        ! cell list        : cell list
        ! inp              : relative position to center cell 
        !                    of inteacting cell
        ! jnp              : relative position to center cell 
        !                    of inteacted cell
        ! nnp              : number of interaction between
        !                    cell and cell
        ! iinter           : index of interaction between
        !                    cell and cell
        !----------------------------------------------------
        
        INTEGER                                 :: num_sub
        INTEGER, DIMENSION(:,:), POINTER        :: num_cell_dim_sub
        TYPE(ppm_type_ptr_to_clist), DIMENSION(:), POINTER :: cell_list
        INTEGER, DIMENSION(:,:), POINTER        :: sub_bcdef
        INTEGER, DIMENSION(:,:), POINTER        :: inp
        INTEGER, DIMENSION(:,:), POINTER        :: jnp
        INTEGER                                 :: nnp
        INTEGER                                 :: iinter
        
        !----------------------------------------------------
        ! Number of cells in 1 and 2 dimension.
        !----------------------------------------------------
        
        INTEGER                         :: n1
        INTEGER                         :: n2
        
        !----------------------------------------------------
      	! Indices about subdomains and cells
        !
        ! idom  : index of subdomains
        ! i*    : cell index in first dimension
        ! j*    : cell index in secondd dimension
        ! k*    : cell index in third dimension
        ! ccell : center cell
        ! icell : interacting cell
        ! jcell : interacted cell
        !----------------------------------------------------
        
        INTEGER                         :: idom
        INTEGER                         :: icstart, icend, i
        INTEGER                         :: jcstart, jcend, j
        INTEGER                         :: kcstart, kcend, k
        INTEGER				:: ccell,icell,jcell
        
        !----------------------------------------------------
        ! Indices about particles in cells
        !
        ! i*   : indices  particles in icell.
        ! j*   : indices  particles in jcell.
        !----------------------------------------------------
        
        INTEGER				:: istart,iend, ipart
        INTEGER                         :: jstart,jend, jpart
        
        !----------------------------------------------------
        ! Indices about particles in data structure,
        ! i.e. in the sense of array, e.g. this%x(ip)
        !----------------------------------------------------
        
        INTEGER				:: ip,jp
        
        !----------------------------------------------------
      	! dij : Inter-particle distance.
        !----------------------------------------------------
        
        REAL(MK), DIMENSION(3)          :: rij
        REAL(MK)                        :: dij
        
        !----------------------------------------------------
        ! kernel parameters :
        !----------------------------------------------------
        
        REAL(MK)                        :: w 
        
        !----------------------------------------------------
        ! Density contribution on particle i or j :
        !----------------------------------------------------
        
        REAL(MK)                        :: rhoi
        REAL(MK)                        :: rhoj
        
        !----------------------------------------------------
      	! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        NULLIFY(bcdef)
        NULLIFY(tboundary)
        wall_rho_type = 0
        
        NULLIFY(colloids)
        coll_rho_type = 0
        
        NULLIFY(num_cell_dim_sub)
        NULLIFY(cell_list)
        NULLIFY(sub_bcdef)
        NULLIFY(inp)
        NULLIFY(jnp)        
        
        !----------------------------------------------------
        ! Parameters of control components : 
        !
        ! rhs_density_type : density formulation;
        ! symmetry         : indication if we use symmetry
        !                    inter-process communication.
        !----------------------------------------------------
        
        rhs_density_type = &
             control_get_rhs_density_type(this%ctrl,stat_info_sub)
        symmetry         = &
             control_get_symmetry(this%ctrl,stat_info_sub)
        
        !----------------------------------------------------
        ! Physics parameters : 
        !
        ! Get ncut off and initial density
        ! from a object of physics class.
        !
        !----------------------------------------------------
        
        num_dim      = this%num_dim
        cut_off      = &
             physics_get_cut_off(this%phys,stat_info_sub)
        cut_off2     = cut_off * cut_off
        init_density = &
             physics_get_rho(this%phys,stat_info_sub)
        
        !----------------------------------------------------
        ! Boundary parameters :
        !----------------------------------------------------
        
        CALL physics_get_bcdef(this%phys,bcdef,stat_info_sub)
        CALL physics_get_boundary(this%phys,tboundary,stat_info_sub)
        
        num_wall_solid = &
             boundary_get_num_wall_solid(tboundary,stat_info_sub)
        IF ( num_wall_solid > 0 ) THEN
           wall_rho_type = &
                boundary_get_wall_rho_type(tboundary,stat_info_sub)
        END IF

        num_part_wall_solid =  &
             boundary_get_num_part_wall_solid(tboundary,stat_info_sub)
        
        !----------------------------------------------------
        ! Colloid parameters :
        !----------------------------------------------------
        
        num_colloid     = &
             physics_get_num_colloid(this%phys,stat_info_sub)
        
        IF ( num_colloid > 0 ) THEN
           
           CALL physics_get_colloid(this%phys,colloids,stat_info_sub)
           
           coll_rho_type = &
                colloid_get_rho_type(colloids, stat_info_sub)
           num_part_colloid = &
                colloid_get_num_numerical_part_tot(colloids, &
                stat_info_sub)
           
        END IF
        
        !----------------------------------------------------
        ! Allocate memory for density calculation :
        !
        ! for symmetry     : allocate for all
        !                    including ghost particles.
        ! for non-symmetry : allocate for only
        !                    real particles
        !----------------------------------------------------
        
        IF(ASSOCIATED(this%rho)) THEN
           DEALLOCATE(this%rho)
        END IF
        
        IF( symmetry ) THEN
           
           ALLOCATE(this%rho(this%num_part_all), &
                STAT=stat_info_sub) 
           this%rho(1:this%num_part_all) = 0.0_MK
           
        ELSE
           
           ALLOCATE(this%rho(this%num_part_real), &
                STAT=stat_info_sub)
           this%rho(1:this%num_part_real) = 0.0_MK
           
        END IF
        
        !----------------------------------------------------
        ! Reset
        !
        ! this%num_part_wall_solid_real and
        ! this%num_part_colloid to zero.
        !
        ! Allocate memory for saving the indices and species
        ! ID of solid wall real boundary particles on current
        ! process.
        !
        ! num_part_wall_solid : the maximum solid wall boundary
        ! particles each process can have.
        !
        ! num_part_colloid    : the maximum colloid boundary 
        ! particles each process can have.
        !----------------------------------------------------
        
        IF ( num_wall_solid > 0 ) THEN
           
           this%num_part_wall_solid_real  = 0
           
           IF ( ASSOCIATED(this%part_wall_solid_real_list)) THEN
              DEALLOCATE(this%part_wall_solid_real_list)
           END IF
           
           ALLOCATE(this%part_wall_solid_real_list(2,num_part_wall_solid))
           
        END IF
        
        IF ( num_colloid > 0 ) THEN
           
           this%num_part_colloid = 0
           
           IF ( ASSOCIATED(this%part_colloid_list)) THEN
              DEALLOCATE(this%part_colloid_list)
           END IF
           
           ALLOCATE(this%part_colloid_list(2,num_part_colloid))
           
        END IF
        
        
        !----------------------------------------------------
        ! Get value of kernel at zero distance.
        !----------------------------------------------------
        
        CALL kernel_kernel(this%kern,0.0_MK,w,stat_info_sub)
        
        !----------------------------------------------------
        ! Initialize density contribution of fluid particle
        ! from itself to itself and set constant density to
        ! non-fluid particle, i.e.,
        ! colloid boundary particles and
        ! solid wall boundary particles created by MCF.
        !
        ! Note that symmetry boundary particles are fluid 
        ! particles in essence, therefore it doesn't need
        ! to assign them constant density, since their values
        ! are from other processes as ghosts.
        !
        ! Note that initial contribution of density only
        ! on real particles, since the ghost ones are 
        ! done on other processes at this point already.
        !
        ! Meanwhile :
        ! Record  wall_solid, colloid boundary real particles.
        !
        ! 1 : mass density.
        ! 2 : number density.
        !----------------------------------------------------
        
        SELECT CASE( rhs_density_type )
           
        CASE (1) ! mass density
           
           DO j = 1, this%num_part_real
              
              IF ( this%id(this%sid_idx,j) == 0 ) THEN
                 
                 !-------------------------------------------
                 ! Fluid particle.
                 !-------------------------------------------
                 
                 this%rho(j) = w * this%m(j)
                 
              ELSE IF ( this%id(this%sid_idx,j) < 0 ) THEN
                 
                 !-------------------------------------------
                 ! Real particles' speceis ID is less than
                 ! zero, it must be wall solid boundary
                 ! particle.
                 ! Give constant density or as fluid particle
                 ! to solid wall boundary particle, depending
                 ! on wall density type.
                 !
                 ! Record the wall boundary particle.
                 !-------------------------------------------
           
                 SELECT CASE ( wall_rho_type )

                 CASE ( mcf_wall_rho_type_constant ) 
                    
                    this%rho(j) = init_density
                 
                 CASE ( mcf_wall_rho_type_dynamic )
                    
                    this%rho(j) = w * this%m(j)
                    
                 END SELECT
                 
                 this%num_part_wall_solid_real = &
                      this%num_part_wall_solid_real + 1
                 
                 this%part_wall_solid_real_list(1, &
                      this%num_part_wall_solid_real) = j
                 this%part_wall_solid_real_list(2, &
                      this%num_part_wall_solid_real) = &
                      this%id(this%sid_idx,j)
                 
              ELSE IF ( this%id(this%sid_idx,j) > 0 ) THEN
                 
                 !-------------------------------------------
                 ! Give constant density or as fliud particle
                 ! to colloid boundary particle, depending on
                 ! colloid density type.
                 !
                 ! Record the colloid boundary particle.
                 !-------------------------------------------
            
                 SELECT CASE ( coll_rho_type )
                    
                 CASE ( mcf_colloid_rho_type_constant )
                    
                    this%rho(j) = init_density
                    
                 CASE ( mcf_colloid_rho_type_dynamic ) 
                    
                    this%rho(j) = w * this%m(j)
                    
                 END SELECT
                 
                 this%num_part_colloid = &
                      this%num_part_colloid + 1
                 
                 this%part_colloid_list(1, &
                      this%num_part_colloid) = j
                 this%part_colloid_list(2, &
                      this%num_part_colloid) = &
                      this%id(this%sid_idx,j)
                 
              END IF ! id(sid_idx,j)
              
           END DO ! j
           
           
        CASE (2) ! number density
           
           DO j = 1, this%num_part_real
              
              IF ( this%id(this%sid_idx,j) == 0 ) THEN
              
                 !-------------------------------------------
                 ! Fluid particle.
                 !-------------------------------------------
                 
                 this%rho(j) = w
                 
              ELSE IF ( this%id(this%sid_idx,j) < 0 ) THEN
                 
                 !-------------------------------------------
                 ! Real particles' speceis ID is less than
                 ! zero, it must be wall solid boundary
                 ! particle.
                 !-------------------------------------------
                 
                 SELECT CASE ( wall_rho_type )
                    
                 CASE ( mcf_wall_rho_type_constant )
                 
                    IF ( this%m(j) < mcf_machine_zero )  THEN
                       PRINT *, "particles_compute_density : ", &
                            j, "s mass is zero !"
                       stat_info = -1
                       GOTO 9999
                    END IF
                    
                    this%rho(j) = init_density / this%m(j)
                 
                 CASE ( mcf_wall_rho_type_dynamic )
                    
                    this%rho(j) = w
                    
                 END SELECT
                 
                 this%num_part_wall_solid_real = &
                      this%num_part_wall_solid_real + 1
                 
                 this%part_wall_solid_real_list(1, &
                      this%num_part_wall_solid_real) = j
                 this%part_wall_solid_real_list(2, &
                      this%num_part_wall_solid_real) = &
                      this%id(this%sid_idx,j)
                 
              ELSE IF ( this%id(this%sid_idx,j) > 0 ) THEN
                 
                 SELECT CASE ( coll_rho_type )
                    
                 CASE ( mcf_colloid_rho_type_constant )
                    
                    IF ( this%m(j) < mcf_machine_zero )  THEN
                       PRINT *, "particles_compute_density : ", &
                            j, "s mass is zero !"
                       stat_info = -1
                       GOTO 9999
                    END IF
                    
                    this%rho(j) = init_density / this%m(j)
                    
                 CASE ( mcf_colloid_rho_type_dynamic )
                    
                    this%rho(j) = w
                    
                 END SELECT
                 
                 this%num_part_colloid = &
                      this%num_part_colloid + 1
                 
                 this%part_colloid_list(1, &
                      this%num_part_colloid) = j
                 this%part_colloid_list(2, &
                      this%num_part_colloid) = &
                      this%id(this%sid_idx,j)
                 
              END IF ! id(sid_idx,j)
              
           END DO ! j
           
        END SELECT ! rhs_density_type
        
        
        !----------------------------------------------------
        ! Get the cell list etc. 
        ! from a object of technique class.
        !----------------------------------------------------
        
        CALL technique_get_cell_list(this%tech,num_sub, &
             num_cell_dim_sub, cell_list, &
             inp,jnp,nnp,sub_bcdef,stat_info_sub)
        
        
        !----------------------------------------------------
        ! Loop over all the sub-domains on this process.
        ! Currently, each process has one subdomain.
        !----------------------------------------------------
        
        DO idom = 1,num_sub
           
           !-------------------------------------------------
           ! For symmetry, cell indices starts from 0
           ! otherwise from 1. 
           ! Both symmetry and non-symmetry cell indices
           ! end at num_cell_dim_sub(*,idom) - 2
           !-------------------------------------------------
           
           IF ( symmetry ) THEN
              
              icstart = 0
              jcstart = 0
              kcstart = 0
              
           ELSE
              
              icstart = 1
              jcstart = 1
              kcstart = 1 
              
           END IF
           
           icend = num_cell_dim_sub(1,idom)-2
           jcend = num_cell_dim_sub(2,idom)-2           
           
           !-------------------------------------------------
           ! Number of cells in first dimension.
           !-------------------------------------------------
           
           n1 = num_cell_dim_sub(1,idom)
           
           !-------------------------------------------------
           ! For 2D, k loops has only one iteration which
           !         exactly means 2D; and n2 has 0 cells.
           ! For 3D, k loops similar way as i,j in 1st and
           !         2nd dimension.
           !         n2 has Number of cells in first and
           !         second dimension, i.e., each plane.
           !-------------------------------------------------
           
           IF ( num_dim == 2 ) THEN
              
              kcend = kcstart
              n2    = 0
              
           ELSE              
              
              kcend = num_cell_dim_sub(3,idom)-2
              n2    = num_cell_dim_sub(1,idom) * &
                   num_cell_dim_sub(2,idom)
              
           END IF
           
           !-------------------------------------------------
           ! Loop over all real cells.
           !-------------------------------------------------
           
           DO k = kcstart, kcend
              
              DO j = jcstart, jcend
                 
                 DO i = icstart, icend
                    
                    !----------------------------------------
                    ! Get index of the center cell.
                    !----------------------------------------
                    
                    ccell = i + 1 + n1 * j + n2 * k
                    
                    !----------------------------------------
                    ! Loop all interactions between cells
                    !----------------------------------------
                    
                    DO iinter = 1, nnp
                       
                       !-------------------------------------
                       ! Get interacting cells, i.e.,
                       ! icell and jcell indices.
                       !-------------------------------------
                       
                       icell = ccell+inp(1,iinter)+ &
                            n1 * inp(2,iinter) + &
                            n2 * inp(3,iinter)
                       
                       jcell = ccell+jnp(1,iinter)+ &
                            n1 * jnp(2,iinter) + &
                            n2 * jnp(3,iinter)
                       
                       !-------------------------------------
                       ! Get pointers of first and last
                       ! particles in icell and jcell.
                       !-------------------------------------
                       
                       istart = cell_list(idom)%lhbx(icell)
                       iend   = cell_list(idom)%lhbx(icell+1)-1
                       
                       jstart = cell_list(idom)%lhbx(jcell)
                       jend   = cell_list(idom)%lhbx(jcell+1)-1
                       
                       
                       !-------------------------------------
                       ! Loop over particles in icell.
                       !-------------------------------------
                       
                       DO ipart = istart, iend
                          
                          !----------------------------------
                          ! Get index of particle in data
                          ! array, i.e., this%x(:,*). 
                          !----------------------------------
                          
                          ip = cell_list(idom)%lpdx(ipart)
                          
                          
                          !----------------------------------
                          ! If ip wall boundary particle
                          ! with constant density,
                          ! for non-symmetric inteaction,
                          ! cycle the loop.
                          !----------------------------------
                          
                          IF( this%id(this%sid_idx,ip) < 0 .AND. &
                               wall_rho_type==&
                               mcf_wall_rho_type_constant .AND. &
                               ( .NOT. symmetry) ) THEN
                             
                             CYCLE
                             
                          END IF
                          
                          !----------------------------------
                          ! If ip colloidal boundary particle
                          ! with constant density,
                          ! for non-symmetric inteaction,
                          ! cycle the loop.
                          !----------------------------------
                          
                          IF( this%id(this%sid_idx,ip) > 0 .AND. &
                               coll_rho_type==&
                               mcf_colloid_rho_type_constant .AND. &
                               ( .NOT. symmetry) ) THEN
                             
                             CYCLE
                             
                          END IF
                          
                          !----------------------------------
                          ! For symmetry case:
                          ! if icell and jcell are the
                          ! same cell, we make sure pair
                          ! particles interaction happen only
                          ! once, i.e., second particle
                          ! comes after first particle.
                          !----------------------------------
                          
                          IF ( jcell == icell .AND. &
                               symmetry ) THEN
                             
                             jstart = ipart + 1 
                             
                          END IF
                          
                          !----------------------------------
                          ! Loop over particles in jcell.
                          !----------------------------------
                          
                          Do jpart = jstart, jend
                             
                             !-------------------------------
                             ! Exclude two same particle.
                             !-------------------------------
                             
                             IF( jcell == icell .AND. &
                                  jpart == ipart ) THEN
                                
                                CYCLE
                                
                             END IF
                             
                             jp = cell_list(idom)%lpdx(jpart) 
                             
                             !-------------------------------
                             ! boundary particles has constant
                             ! rho, no contribution to each
                             ! other.
                             !-------------------------------
                             
                             IF ( this%id(this%sid_idx,ip) /= 0 .AND. &
                                  this%id(this%sid_idx,jp) /= 0 .AND. &
                                  wall_rho_type==mcf_wall_rho_type_constant .AND. &
                                  coll_rho_type==mcf_colloid_rho_type_constant ) THEN
                                CYCLE
                             END IF
                             
                             !-------------------------------
                             ! Dist. of particles ip and jp.
                             !-------------------------------
                             
                             rij(1:num_dim) = &
                                  this%x(1:num_dim,ip) - &
                                  this%x(1:num_dim,jp)
                             
                             dij = &
                                  DOT_PRODUCT(rij(1:num_dim), rij(1:num_dim))
                             
                             !-------------------------------
                             ! Skip 2 particles beyond cut off.
                             !-------------------------------
                             
                             IF ( dij >= cut_off2 ) THEN
                                CYCLE
                             ELSE
                                dij = SQRT(dij)
                             END IF
                             
                             !-------------------------------
                             ! Get kernel value and
                             ! calculate density contribution.
                             !-------------------------------
                             
                             CALL kernel_kernel(this%kern,dij,&
                                  w,stat_info_sub)
                             CALL rhs_density_ff(this%rhs, w,&
                                  this%m(ip),this%m(jp), &
                                  rhoi,rhoj,stat_info_sub)
                             
                             !-------------------------------
                             ! Density calculated for
                             ! fluid particles, or
                             ! boundary particle without
                             ! constant density type.
                             !-------------------------------
                             
                             IF ( this%id(this%sid_idx,ip) == 0 .OR. &
                                  ( this%id(this%sid_idx,ip) < 0 .AND. &
                                  wall_rho_type==mcf_wall_rho_type_dynamic ) .OR. &
                                  ( this%id(this%sid_idx,ip) > 0 .AND. &
                                  coll_rho_type==mcf_colloid_rho_type_dynamic ) ) THEN
                                
#if 0
                                IF ( this%id(this%sid_idx,ip) /=0 .AND. &
                                     this%id(this%sid_idx,jp) == &
                                     this%id(this%sid_idx,ip) ) THEN
                                   CYCLE
                                END IF
#endif
                                
                                this%rho(ip)= this%rho(ip) + rhoi
                                
                             END IF
                             
                             !-------------------------------
                             ! For symmetry, jp also gets
                             ! calculated(if it is fluid).
                             !-------------------------------
                             
                             IF ( symmetry .AND. &
                                  ( this%id(this%sid_idx,jp) == 0 .OR. &
                                  ( this%id(this%sid_idx,jp) < 0 .AND. &
                                  wall_rho_type ==mcf_wall_rho_type_dynamic ) .OR. &
                                  ( this%id(this%sid_idx,jp) > 0 .AND. &
                                  coll_rho_type==mcf_colloid_rho_type_dynamic ) ) ) THEN
                             
#if 0   
                                IF ( this%id(this%sid_idx,jp) /=0 .AND. &
                                     this%id(this%sid_idx,ip) == &
                                     this%id(this%sid_idx,jp) ) THEN
                                   
                                   CYCLE

                                END IF
#endif                           
                                this%rho(jp)= this%rho(jp) + rhoj
                                
                             END IF
                             
                          END DO ! jpart
                          
                       END DO ! ipart
                       
                    END DO ! iinter : 1, nnp
                    
                 END DO ! i : icstart, icend
                 
              END DO ! j: jcstart, jcend
              
           END DO ! k:  kcstart, kcend
           
        END DO ! idom : 1,num_sub       
        
        
      	!----------------------------------------------------
      	! Return.
      	!----------------------------------------------------
        
        
9999	CONTINUE
        
        IF(ASSOCIATED(bcdef)) THEN
           DEALLOCATE(bcdef)
        END IF
        
        IF(ASSOCIATED(num_cell_dim_sub)) THEN
           DEALLOCATE(num_cell_dim_sub)
        END IF
        
        IF(ASSOCIATED(cell_list)) THEN
           DEALLOCATE(cell_list)
        END IF
        
        IF(ASSOCIATED(sub_bcdef)) THEN
           DEALLOCATE(sub_bcdef)
        END IF
        
        IF(ASSOCIATED(inp)) THEN
           DEALLOCATE(inp)
        END IF
        
        IF(ASSOCIATED(jnp)) THEN
           DEALLOCATE(jnp)
        END IF
        
        RETURN
        
      END SUBROUTINE particles_compute_density

