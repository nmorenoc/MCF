      SUBROUTINE particles_compute_interaction(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_compute_interaction
        !----------------------------------------------------
        !
        ! Purpose     : Computing interaction for 2D or 3D, 
        !               using symmetric or non-symmetric 
        !               inter-process communication.
        !
        !               Interaction includes :
        !               force;
        !               velocity gradient tensor in case of
        !               non-Newtonian Oldroyd-B model.
        !                  
        !
        ! Routines    : pp_interaction.inc
        !
        !
        ! References  : Sbalzarini et al.
        !               2006 Journal of Computational Physics.
        !
        ! Remarks     :  
        !
        ! Revisions   : V0.7 22.04.2010, changes in the 
        !               density of the wall particles in 
        !               the case of non newtonian
        !               fluids. For every pairwise f-w, the 
        !               density of the wall particle is 
        !               the same than the density of the
        !               fluid particle. I think this should
        !               be implemented also for the particles
        !               of the colloids, and also for the 
        !               newtonian case. (Adolfo Vázquez-Quesada)
        !
        !               V0.6 03.12 2009, move the allocation
        !               and recording boundary particles
        !               to particles_set_boundary_IDt() and
        !               particles_compute_density().
        !
        !               V0.5 09.10 2009, merge 2D,3D, 
        !               symmetry, non-symmetry together. 
        !               
        !               V0.4 29.09 2009, include wall boundary
        !               particle interaction.
        !                 
        !               V0.3 28.07 2009, merge vgt calculation
        !               with force in this loop.
        !
        !               V0.2 28.09 2009, merge force and vgt
        !               calculation together. 
        !               In general, all particle interactions
        !               should be done here.
        !
        !               V0.1 16.03 2009, original
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
        ! symmetry      : indicate if we use symmetric
        !                 inter-process communication or not.
        ! Newtonian     : if fluid is Newtonian.
        ! Browonian     : if fluid is Brownian.        
        ! p_energy      : if potential energy is needed.  
        !----------------------------------------------------
        
        INTEGER                         :: rhs_density_type
        LOGICAL                         :: symmetry
        LOGICAL                         :: Newtonian
        LOGICAL                         :: Brownian
        LOGICAL                         :: p_energy
        
        !----------------------------------------------------
        ! Physics parameters :
        !
        ! num_species    : number of species.
        ! num_dim        : number of dimension.
        ! cut_off        : compact support domain.
        ! cut_off2       : cut_off * cut_off.
        !
        ! bcdef          : boundary condition definition.
        ! boundary       : boundary object pointer.
        ! num_wall_sym   : number of wall boundaries,
        !                  created by PPM using symmetry.
        ! num_wall_solid : number of wall boundaries,
        !                  created by MCF using solid.
        !
        ! num_wall       : number of wall boundaries,
        !                  in general, either symmetry
        !                  or solid.
        ! num_shear      : number of shear boundaries.
        !
        ! num_colloid    : number of colloidal particle.
        ! colloids       : Colloids object pointer.
        !----------------------------------------------------
        
        INTEGER                         :: num_species
        INTEGER                         :: num_dim
        REAL(MK)                        :: cut_off
        REAL(MK)                        :: cut_off2
        REAL(MK), DIMENSION(:),POINTER  :: dx
        REAL(MK)                        :: eta
        REAL(MK), DIMENSION(3)          :: uij
        REAL(MK)                        :: h

        !----------------------------------------------------
        ! Boundary parameters :
        !
        ! num_wall_solid : number of solid wall boundaries.
        !----------------------------------------------------
       
        INTEGER, DIMENSION(:), POINTER  :: bcdef
        TYPE(Boundary), POINTER         :: tboundary
        INTEGER                         :: wall_noslip
        INTEGER                         :: num_sym
        INTEGER                         :: num_wall_sym
        INTEGER                         :: num_wall_solid
        INTEGER                         :: num_le
        INTEGER                         :: num_shear
        INTEGER                         :: num_colloid
        TYPE(Colloid), POINTER          :: colloids
        INTEGER                         :: coll_noslip
        
        !----------------------------------------------------
        ! Counters / Indices
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
      	! Inter-particle distance.
        !----------------------------------------------------
        
        REAL(MK), DIMENSION(3)          :: rij
        REAL(MK)                        :: dij
        
        !----------------------------------------------------
        ! kernel parameters.
        !----------------------------------------------------
        
        REAL(MK)                        :: w 
        REAL(MK)                        :: gradW
        
        
        !----------------------------------------------------
        ! For Non-Newtonian viscoelastic
        ! Oldroyd-B model
        !
        ! a      : index
        ! b      : index
        !
        !----------------------------------------------------
        
        INTEGER                         :: a, b
        
        !----------------------------------------------------
        ! fip, fjp : pair-wise force
        !----------------------------------------------------
        
        REAL(MK), DIMENSION(3)          :: fip
        REAL(MK), DIMENSION(3)          :: fjp
#ifdef __PARTICLES_FORCE_SEPARATE
        REAL(MK), DIMENSION(3)          :: fpip
        REAL(MK), DIMENSION(3)          :: fpjp
        REAL(MK), DIMENSION(3)          :: fvip
        REAL(MK), DIMENSION(3)          :: fvjp
        REAL(MK), DIMENSION(3)          :: frip
        REAL(MK), DIMENSION(3)          :: frjp        
#endif
        
        !----------------------------------------------------
        ! First try for colloid-colloid interaction.
        !
        ! v_ip, v_jp : temparory velocity of boundary 
        !              particles, used of no-slip interpolation.
        ! rho_ip,rho_jp: desnity of boundary 
        !              particles, used for two boundary 
        !              particles interaction.      
        ! p_ip, p_jp : pressure of boundary 
        !              particles, used for two boundary 
        !              particles interaction.      
        !----------------------------------------------------
        
        REAL(MK), DIMENSION(3)          :: v_ip
        REAL(MK), DIMENSION(3)          :: v_jp
        REAL(MK)                        :: rho_ip
        REAL(MK)                        :: rho_jp        
!        REAL(MK)                        :: p_ip
!        REAL(MK)                        :: p_jp

  	!----------------------------------------------------
      	! Initialization of variables.
      	!----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0    
        
        NULLIFY(bcdef)
        NULLIFY(tboundary)
        wall_noslip = 1
        
        NULLIFY(colloids)
        coll_noslip = 1
        
        NULLIFY(num_cell_dim_sub)
        NULLIFY(cell_list)
        NULLIFY(sub_bcdef)        
        NULLIFY(inp)
        NULLIFY(jnp)
        
        !----------------------------------------------------
        ! Control parameters :
        !
        ! Get control variables.
        !----------------------------------------------------
        
        rhs_density_type = &
             control_get_rhs_density_type(this%ctrl,stat_info_sub)
        symmetry  = &
             control_get_symmetry(this%ctrl,stat_info_sub)
        Newtonian = &
             control_get_Newtonian(this%ctrl,stat_info_sub)
        Brownian = &
             control_get_Brownian(this%ctrl,stat_info_sub)      
        p_energy  = &
             control_get_p_energy(this%ctrl,stat_info_sub)
        
        !----------------------------------------------------
        ! Physics parameters :
        !
        ! from a object of Physics class.
        !
        !----------------------------------------------------
        
        num_species = &
             physics_get_num_species(this%phys,stat_info_sub)
        num_dim     = &
             physics_get_num_dim(this%phys,stat_info_sub)
        cut_off     = &
             physics_get_cut_off(this%phys,stat_info_sub)
        cut_off2 = cut_off * cut_off
        NULLIFY(dx)
        CALL physics_get_dx(this%phys,dx,stat_info_sub)
        eta         = &
             physics_get_eta(this%phys,stat_info_sub)
        !----------------------------------------------------
        ! Boundary parameters :
        !
        ! Get boundary conditions.
        !----------------------------------------------------
        
        CALL physics_get_bcdef(this%phys,bcdef,stat_info_sub)
        CALL physics_get_boundary(this%phys,&
             tboundary,stat_info_sub)
        
        wall_noslip    = &
             boundary_get_wall_noslip_type(tboundary,stat_info_sub)
        num_sym        = &
             boundary_get_num_sym(tboundary,stat_info_sub)
        num_wall_sym   = &
             boundary_get_num_wall_sym(tboundary,stat_info_sub)
        num_wall_solid = &
             boundary_get_num_wall_solid(tboundary,stat_info_sub)
        num_le         = &
             boundary_get_num_le(tboundary,stat_info_sub)
        num_shear      = &
             boundary_get_num_shear(tboundary,stat_info_sub)
        
        !----------------------------------------------------
        ! Colloid parameters :
        !----------------------------------------------------
        
        num_colloid = &
             physics_get_num_colloid(this%phys,stat_info_sub)
        
        IF ( num_colloid > 0 ) THEN
           
           CALL physics_get_colloid(this%phys,&
                colloids,stat_info_sub)
           coll_noslip = &
                colloid_get_noslip_type(colloids,stat_info_sub)
           
        END IF
        
        !----------------------------------------------------
        ! Allocate memory for force,
        !
        ! for symmetry     : allocate num_part_all
        ! for non-symmetry : allocate num_part_real
        !
        ! For symmetry pair interaction/inter-process
        ! communication, we need to allocate for 
        ! ghosts also.
        !
        ! Note that :
        ! For Wall using PPM symmetry boundaries, wall using
        ! solid boundary particles or Lees-Edwards boundaries, 
        ! we need the forces on ghosts particles,
        ! in order calculate ensemble later.
        !----------------------------------------------------
        
        IF (ASSOCIATED(this%f)) THEN
           DEALLOCATE(this%f,STAT=stat_info_sub)
        END IF
     
        IF(  symmetry .OR. &
             num_wall_sym > 0 .OR. num_le > 0 ) THEN
           
           ALLOCATE(this%f(num_dim,this%num_part_all), &
                STAT=stat_info_sub)
           
        ELSE
           
           ALLOCATE(this%f(num_dim,this%num_part_real), &
                STAT=stat_info_sub)
           
        END IF
        
        this%f(:,:) = 0.0_MK

#ifdef __PARTICLES_FORCE_SEPARATE
        IF (ASSOCIATED(this%fp)) THEN
           DEALLOCATE(this%fp,STAT=stat_info_sub)
        END IF
        IF (ASSOCIATED(this%fv)) THEN
           DEALLOCATE(this%fv,STAT=stat_info_sub)
        END IF
        IF (ASsOCIATED(this%fr)) THEN
           DEALLOCATE(this%fr,STAT=stat_info_sub)
        END IF
        
        IF(  symmetry .OR. &
             num_wall_sym > 0 .OR. num_le > 0 ) THEN
           
           ALLOCATE(this%fp(num_dim,this%num_part_all), &
                STAT=stat_info_sub)
           ALLOCATE(this%fv(num_dim,this%num_part_all), &
                STAT=stat_info_sub)
           ALLOCATE(this%fr(num_dim,this%num_part_all), &
                STAT=stat_info_sub)
           
        ELSE
           
           ALLOCATE(this%fp(num_dim,this%num_part_real), &
                STAT=stat_info_sub)
           ALLOCATE(this%fv(num_dim,this%num_part_real), &
                STAT=stat_info_sub)
           ALLOCATE(this%fr(num_dim,this%num_part_real), &
                STAT=stat_info_sub)
           
        END IF
        
        this%fp(:,:) = 0.0_MK
        this%fv(:,:) = 0.0_MK
        this%fr(:,:) = 0.0_MK
#endif 
        
        !----------------------------------------------------
        ! Allocate memory for veocity gradient 
        ! tensor in case of non-Newtonian Oldroyd-B
        ! viscoelastic model being used.
        ! 
        ! for symmetry     : allocate num_part_all
        ! for non-symmetry : allocate num_part_real
        !----------------------------------------------------
        
        IF ( .NOT. Newtonian ) THEN
           
           IF(ASSOCIATED(this%vgt)) THEN
              DEALLOCATE(this%vgt)
           END IF
           
           IF(  symmetry ) THEN
              
              ALLOCATE(this%vgt(num_dim**2,this%num_part_all),&
                   STAT=stat_info_sub)
              
           ELSE
              
              ALLOCATE(this%vgt(num_dim**2,this%num_part_real),&
                   STAT=stat_info_sub)
              
           END IF
           
           this%vgt(:,:) = 0.0_MK
           
        END IF
        
        !----------------------------------------------------
        ! Allocate memory for acceleration of
        ! potential energy in case need potential
        ! energy.
        ! 
        ! for symmetry     : allocate num_part_all
        ! for non-symmetry : allocate num_part_real
        !----------------------------------------------------
    
        IF( p_energy ) THEN
        
           IF (ASSOCIATED(this%au)) THEN
              DEALLOCATE(this%au,STAT=stat_info_sub)
           END IF
           
           IF(  symmetry ) THEN
              
              ALLOCATE(this%au(this%num_part_all), &
                   STAT=stat_info_sub)
              
           ELSE
              
              ALLOCATE(this%au(this%num_part_real), &
                   STAT=stat_info_sub)
              
           END IF
           
           this%au(:) = 0.0_MK
           
        END IF
        
        
        !----------------------------------------------------
        ! Get the cell list from a object of technique class.
        !----------------------------------------------------
        
        CALL technique_get_cell_list(this%tech,&
             num_sub,num_cell_dim_sub,&
             cell_list,inp,jnp,nnp,sub_bcdef,stat_info_sub)
        

        !----------------------------------------------------
        ! Loop over all the sub-domains on this process.
        !----------------------------------------------------
        
        DO idom = 1, num_sub
           
           !-------------------------------------------------
           ! For symmetry, cell indices starts from 0
           ! otherwise from 1 for real particles.
           ! Both symmetry and non-symmetry cell indices
           ! end at num_cell_dim_sub(*,idom) - 2
           !
           ! However, for symmetry inter-process communication,
           ! if the subdomain is at lower boundary and
           ! we have symmetry, wall using symmetry,
           ! or Lees-Edwars boundary, cell indices starting
           ! from 0 can be a problem.
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
           
           n1  = num_cell_dim_sub(1,idom)
           
           !-------------------------------------------------
           ! For 2D, k loops has only one iteration which
           !         means 2D; and n2 has 0 cells.
           ! For 3D, k loops similar way as i,j in 1st and
           !         2nd dimension.
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
           ! Loop over all real cells
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
                          ! array, i.e., this%x. 
                          !----------------------------------
                          
                          ip = cell_list(idom)%lpdx(ipart)
                          
                          !----------------------------------
                          ! For symmetry case:
                          ! if icell and jcell are the
                          ! same cell, we make sure pair
                          ! particles interaction happen only
                          ! once. 
                          !----------------------------------
                          
                          IF ( jcell == icell .AND. &
                               symmetry ) THEN
                             
                             jstart = ipart + 1 
                             
                          END IF
                          
                          !----------------------------------
                          ! Loop over particles in jcell.
                          !----------------------------------
                          
                          DO jpart = jstart, jend
                             
                             !-------------------------------
                             ! Exclude 2 same particles.
                             !-------------------------------
                             
                             IF( jcell == icell .AND. &
                                  jpart == ipart ) THEN
                                
                                CYCLE
                                
                             END IF
                             
                             jp = cell_list(idom)%lpdx(jpart) 
                             
                             !-------------------------------
                             ! Since two particles index
                             ! are known,
                             ! their interactions are 
                             ! calculated in routine/file
                             ! pp_interactoin.inc.
                             !
                             ! Initialize force.
                             !-------------------------------
                             
                             fip (1:num_dim) = 0.0_MK
                             fjp (1:num_dim) = 0.0_MK
#ifdef __PARTICELS_FORCE_SEPARATE
                             fpip (1:num_dim) = 0.0_MK
                             fpjp (1:num_dim) = 0.0_MK
                             fvip (1:num_dim) = 0.0_MK
                             fvjp (1:num_dim) = 0.0_MK
                             frip (1:num_dim) = 0.0_MK
                             frjp (1:num_dim) = 0.0_MK
                             
#endif
                             
#include "pp_interaction.inc"
                             
                             !----------------------------------
                             ! Add up force of jp acting on ip.
                             !----------------------------------
                             
                             this%f(1:num_dim,ip) = &
                                  this%f(1:num_dim,ip) + fip(1:num_dim)
                             
                             !----------------------------------
                             ! 1)In symmetry inter-communiction,
                             ! add force on jp, no matter the 
                             ! particle is fluid, symmetry, 
                             ! wall using symmetry, solid wall,
                             ! Lees-Edwards or
                             ! colloid boundary particle.
                             !
                             ! 2)In non-symmetry inter-communication,
                             ! add force on jp when it is
                             ! wall_sym.
                             !
                             ! Add fj to particle jp.
                             !----------------------------------
                             
                             IF ( symmetry .OR. &
                                  ( this%id(this%sid_idx,jp) < 0  .AND. &
                                  num_wall_sym > 0 ) ) THEN
                                
                                this%f(1:num_dim,jp) = &
                                     this%f(1:num_dim,jp) + fjp(1:num_dim)
                                
                             END IF
                             
#ifdef __PARTICLES_FORCE_SEPARATE
                             this%fp(1:num_dim,ip) = &
                                  this%fp(1:num_dim,ip) + fpip(1:num_dim)
                             this%fv(1:num_dim,ip) = &
                                  this%fv(1:num_dim,ip) + fvip(1:num_dim)
                             this%fr(1:num_dim,ip) = &
                                  this%fr(1:num_dim,ip) + frip(1:num_dim)
                            
                             IF ( symmetry .OR. &
                                  ( this%id(this%sid_idx,jp) < 0  .AND. &
                                  num_wall_sym > 0 ) ) THEN
                                
                                this%fp(1:num_dim,jp) = &
                                     this%fp(1:num_dim,jp) + fpjp(1:num_dim)
                                this%fv(1:num_dim,jp) = &
                                     this%fv(1:num_dim,jp) + fvjp(1:num_dim)
                                this%fr(1:num_dim,jp) = &
                                     this%fr(1:num_dim,jp) + frjp(1:num_dim)
                                
                             END IF
#endif
                             
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
        
        IF( ASSOCIATED(bcdef)) THEN
           DEALLOCATE(bcdef)
        END If
        
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
        
        IF(ASSOCIATED(dx)) THEN
           DEALLOCATE(dx)
        END IF

        RETURN
        
      END SUBROUTINE particles_compute_interaction
      
      
      
