      SUBROUTINE particles_normalize_density(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particls_normalize_density
        !----------------------------------------------------
        !
        ! Purpose     : Normalize 2D/3D density.
        !	 	      	 
        !                
        ! Reference   :
        !
        ! Remark      : When there are colloidal or solid 
        !               wall boundary particles in the system,
        !               they are created as thin layer at
        !               the boundaries. Therefore,
        !               either density are kept constant
        !               or summations over neighbours are
        !               taken as for fluid particles.
        !               However, for 2nd case, the
        !               neighbours are not suffcient, since
        !               the colloids are hollow inside or
        !               the walls are at the boundary of
        !               domain, renormalization of this
        !               summantion must be taken.
        !
        !               In a word, target of normalization
        !               for density is non-fluid particles,
        !               i.e., colloidal or solid wall
        !               boundary particles.
        !               
        !
        ! Revisions   : V0.1 06.05.2010, original version.
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
        INTEGER                         :: num_wall_solid
        
        !----------------------------------------------------
        ! Colloid parameters :
        !
        ! num_colloid      : number of colloidal particle.
        ! num_part_colloid : number of colloid boundary 
        !                    particles.       
        !----------------------------------------------------
        
        INTEGER                         :: num_colloid
        
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
        ! d   : index for loop over different dimension.
        !----------------------------------------------------
        
        REAL(MK)                        :: dij
        INTEGER                         :: d
        
        !----------------------------------------------------
        ! kernel parameters :
        !----------------------------------------------------
        
        REAL(MK)                        :: w 
        
        !----------------------------------------------------
      	! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        NULLIFY(bcdef)
        NULLIFY(tboundary)
        
        
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
        ! Get cut off and initial density
        ! from a object of physics class.
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
        
        CALL physics_get_boundary(this%phys,tboundary,stat_info_sub)
        
        num_wall_solid      = &
             boundary_get_num_wall_solid(tboundary,stat_info_sub)
        
        !----------------------------------------------------
        ! Colloid parameters :
        !----------------------------------------------------
        
        num_colloid     = &
             physics_get_num_colloid(this%phys,stat_info_sub)
        
        
        !----------------------------------------------------
        ! Get value of kernel at zero distance.
        !----------------------------------------------------
        
        CALL kernel_kernel(this%kern,0.0_MK,w,stat_info_sub)
        
  
        !----------------------------------------------------
        ! Allocate memory for density normalization :
        !
        ! for symmetry     : allocate for all
        !                    including ghost particles.
        ! for non-symmetry : allocate for only
        !                    real particles
        !
        ! Initial value set to one, since for fluid particle
        ! it is 1, i.e., no need for normalization.
        !----------------------------------------------------
        
        IF(ASSOCIATED(this%rho_norm)) THEN
           DEALLOCATE(this%rho_norm)
        END IF
        
        IF( symmetry ) THEN
           
           ALLOCATE(this%rho_norm(this%num_part_all), &
                STAT=stat_info_sub) 
           this%rho_norm(1:this%num_part_all) = 1.0_MK
           
        ELSE
           
           ALLOCATE(this%rho_norm(this%num_part_real), &
                STAT=stat_info_sub)
           this%rho_norm(1:this%num_part_real) = 1.0_MK
           
        END IF

        !----------------------------------------------------
        ! Initialize  normalization contribution of 
        ! non-fluid particle from itself to itself.
        !
        ! Note that symmetry boundary particles are fluid 
        ! particles in essence, therefore it doesn't need
        ! be to assigned anything, since their values
        ! are from other processes as ghosts.
        !
        ! Note that initial contribution of density only
        ! on real particles, since the ghost ones are 
        ! done on other processes at this point already.
        !
        !
        ! 1 : mass density.
        ! 2 : number density.
        !----------------------------------------------------
        
        
        SELECT CASE( rhs_density_type )
           
        CASE (1) ! mass density formulation
           
           !-------------------------------------------------
           ! For non-fluid particle : 
           ! colloid or solid wall boundary particle.
           !-------------------------------------------------
           
           IF ( num_wall_solid > 0 ) THEN
              
              DO i = 1, this%num_part_wall_solid_real
                 
                 ip = this%part_wall_solid_real_list(1,i)
                 
                 !this%rho_norm(ip) = this%m(ip) * w / this%rho(ip)
                 this%rho_norm(ip) = 0.0_MK
                 
              END DO
              
           END IF
           
           IF ( num_colloid > 0 ) THEN
              
              DO i = 1, this%num_part_colloid
                 
                 ip = this%part_colloid_list(1,i)
                 
                 !this%rho_norm(ip) = this%m(ip) * w / this%rho(ip)
                 this%rho_norm(ip) = 0.0_MK

              END DO
              
           END IF
           
        CASE (2) ! number density formulation
           
           !-------------------------------------------------
           ! For non-fluid particle : 
           ! colloid or solid wall boundary particle.
           !-------------------------------------------------
           
           IF ( num_wall_solid > 0 ) THEN
              
              DO i = 1, this%num_part_wall_solid_real
                 
                 ip = this%part_wall_solid_real_list(1,i)
                 
                 !this%rho_norm(ip) = w / this%rho(ip)
                 this%rho_norm(ip) = 0.0_MK

                 
              END DO
              
           END IF
           
           IF ( num_colloid > 0 ) THEN
              
              DO i = 1, this%num_part_colloid
                 
                 ip = this%part_colloid_list(1,i)
                 
                 !this%rho_norm(ip) = w / this%rho(ip)
                 this%rho_norm(ip) = 0.0_MK

                 
              END DO
              
           END IF
           
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
           
           n1    = num_cell_dim_sub(1,idom)
           
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
                          ! If ip is fluid particle,
                          ! for non-symmetric inteaction,
                          ! cycle the loop, since fluid
                          ! particle ip has either no need
                          ! to normalize or desn't contribute
                          ! to other non-fluid particle's 
                          ! normalization.
                          !-----------------------------------
                          
                          IF( this%id(this%sid_idx,ip) == 0 .AND. &
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
                             ! Dist. of particles ip and jp.
                             !-------------------------------
                             
                             dij = 0.0_MK
                             
                             DO  d = 1, num_dim
                                
                                dij = dij + ( this%x(d,ip) - &
                                     this%x(d,jp) ) ** 2
                                
                             END DO
                             
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
                             
                             !-------------------------------
                             ! Density only calculated for
                             ! non-fluid particles.
                             !-------------------------------
                             
                             IF ( this%id(this%sid_idx,ip) /= 0 ) THEN
                                
                                       
                                IF ( this%id(this%sid_idx,jp) /= &
                                     this%id(this%sid_idx,ip) ) THEN
                                   
                                   SELECT CASE (rhs_density_type)
                                      
                                   CASE (1)
                                      
                                      this%rho_norm(ip) = this%rho_norm(ip) + &
                                           this%m(jp) * w / this%rho(jp)
                                      
                                   CASE (2)
                                      
                                      this%rho_norm(ip)= this%rho_norm(ip) + &
                                           w / this%rho(jp)
                                      
                                   END SELECT ! rhs_density_type

                                END IF
                                
                             END IF
                             
                             !-------------------------------
                             ! For symmetry, jp also gets
                             ! calculated(if it is fluid).
                             !-------------------------------
                             
                             IF ( symmetry .AND. &
                                  this%id(this%sid_idx,jp) /= 0 ) THEN
                              
                                IF ( this%id(this%sid_idx,ip) /= &
                                     this%id(this%sid_idx,jp) ) THEN
                                   
                                   SELECT CASE (rhs_density_type)
                                      
                                   CASE (1)
                                      
                                      this%rho_norm(jp) = this%rho_norm(jp) + &
                                           this%m(ip) * w / this%rho(ip)
                                      
                                   CASE (2)
                                      
                                      this%rho_norm(jp)= this%rho_norm(jp) + &
                                           w / this%rho(ip)
                                      
                                   END SELECT ! rhs_density_type
                                
                                END IF
                                
                             END IF ! symmetry and sid /=0
                             
                          END DO ! jpart
                          
                       END DO ! ipart
                       
                    END DO ! iinter : 1, nnp
                    
                 END DO ! i : icstart, icend
                 
              END DO ! j: jcstart, jcend
              
           END DO ! k:  kcstart, kcend
           
        END DO ! idom : 1,num_sub       
        

        !----------------------------------------------------
        ! Normalize density for wall boundary particles
      	!----------------------------------------------------

        IF ( num_wall_solid > 0 ) THEN
           
           DO i = 1, this%num_part_wall_solid_real
              
              ip = this%part_wall_solid_real_list(1,i)
              
              IF ( this%rho_norm(ip) > 0.0_MK ) THEN

                 this%rho(ip) = this%rho(ip) / this%rho_norm(ip)
                 
              END IF
              
           END DO
           
        END IF
        
        !----------------------------------------------------
        ! Normalize density for colloidal boundary particles
      	!----------------------------------------------------

        IF ( num_colloid > 0 ) THEN
           
           DO i = 1, this%num_part_colloid
              
              ip = this%part_colloid_list(1,i)
              
              IF ( this%rho_norm(ip) > 0.0_MK ) THEN

                 this%rho(ip) = this%rho(ip) / this%rho_norm(ip)
              
              END IF
              
           END DO
           
        END IF
     
        
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
        
      END SUBROUTINE particles_normalize_density

