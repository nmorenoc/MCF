      SUBROUTINE particles_reset_boundary_ghost_interaction(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_reset_boundary_interaction
        !----------------------------------------------------
        !
        ! Purpose     : Reset symmetry, wall using symmetry
        !               force to zero.
        !               Reset velocity gradient tensor
        !               if it is Oldrody-B model.
        
        ! Refernece   :
        !
        ! Remark      : Interaction means force, velocity
        !               gradient tensor, etc, i.e.,
        !               general pair-wise interactions.
        !
        !               A For periodic boundary particles,
        !               we do nothing.
        !
        !               B If we are using symmetry 
        !               inter-process communication :
        !
        !               Symmetry or Wall using symmetry
        !               boundary, which consist of only
        !               ghost particles.
        !               
        !               For symmetry boundary particles, 
        !               wall using symmetry/mirror boundary
        !               particles, they are all created by
        !               PPM. They are outside of physical 
        !               domain and all of them are treated
        !               completely as ghost particles, we 
        !               don't need to integrate their 
        !               position at all. Moreover, in 
        !               symmetry inter-process communication
        !               case, the forces on them don't 
        !               contribute to their hosts at all.
        !               Therefore, we reset their force to
        !               zero here, in order to have no effect
        !               on their hosts.
        !
        !               C If there is Wall boundary using solid 
        !               particles, no matter if it is
        !               symmetry inter-process communication:
        !
        !               Since solid wall are modeled using 
        !               solid particles inside the total 
        !               physical domain, i.e., there is no
        !               relative motion among solid wall 
        !               boundary particles if they are
        !               inside the same wall and PPM treat
        !               the solid wall particles as inside
        !               physical domain.
        !
        !               Moreover, the solid wall boundary 
        !               particles are created by client
        !               MCF, thus we have to integrate 
        !               their(real particles) positions
        !               every step.
        !
        !               However, during the interaction(force)
        !               calculation, the each solid wall 
        !               boundary particle  was  given 
        !               an individual force,
        !               therefore we have to reset this
        !               force to zero, in order to make the
        !               solid wall boundary particle have 
        !               a prescribed velocity as rigid body,
        !               i.e., no individual acceleartion.
        !
        !               D Lees-Edward boundary.
        !
        !                We don't do anything here, in
        !                symmetry inter-process communication
        !                case, forces on these boundary 
        !                particles have contributions to
        !                their host particles.
        !
        ! Revision    : V0.2 04.12 2009, considering all 
        !               available boundary conditions.
        !
        !               V0.1 22.10 2009, original version.
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
        !
        ! this       : an object of Particles Class.
        ! stat_info  : return flag of status.
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)  :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        
	!----------------------------------------------------
        ! Local variables
	!----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        LOGICAL                         :: symmetry
        LOGICAL                         :: Newtonian
        INTEGER                         :: num_dim
        INTEGER                         :: num_dim2
        TYPE(Boundary), POINTER         :: tboundary
        INTEGER                         :: num_sym
        INTEGER                         :: num_wall_sym
        INTEGER                         :: num_wall_solid
        INTEGER                         :: j
        INTEGER                         :: ip
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        NULLIFY(tboundary)

        !----------------------------------------------------
        ! Control parameters :
        !----------------------------------------------------
        
        symmetry  = &
             control_get_symmetry(this%ctrl,stat_info_sub)
        Newtonian = &
             control_get_Newtonian(this%ctrl,stat_info_sub)
        
        !----------------------------------------------------
        ! Physics parameters.
        !----------------------------------------------------
     
        num_dim  = this%num_dim
        num_dim2 = num_dim * num_dim
        
        CALL physics_get_boundary(this%phys,tboundary,&
             stat_info_sub)
        
        num_sym      = &
             boundary_get_num_sym(tboundary,stat_info_sub)      
        num_wall_sym = &
             boundary_get_num_wall_sym(tboundary,stat_info_sub)
        num_wall_solid = &
             boundary_get_num_wall_solid(tboundary,stat_info_sub)
        
        
        IF ( symmetry ) THEN
           
           !-------------------------------------------------
           ! Reset ghost particles which constitute symmetry
           ! boundary with zero force.
           !-------------------------------------------------
           
           IF ( num_sym > 0 ) THEN
              
              DO j =1, this%num_part_sym
                 
                 ip  = this%part_sym_list(1,j)
                 this%f(1:num_dim,ip) = 0.0_MK       
                 
              END DO
              
           END IF ! num_sym
           
           !-------------------------------------------------
           ! Reset ghost particles which constitute wall 
           ! using symmetry/mirror boundary with zero force.
           !-------------------------------------------------
           
           IF ( num_wall_sym > 0 ) THEN
              
              DO j =1, this%num_part_wall_sym
                 
                 ip  = this%part_wall_sym_list(1,j)
                 this%f(1:num_dim,ip) = 0.0_MK       
                 
              END DO
              
           END IF ! num_wall_sym
        
           !-------------------------------------------------
           ! Reset ghost particles which constitute wall 
           ! using solid boundary with zero force.
           !-------------------------------------------------
           
           IF ( num_wall_solid > 0 ) THEN
              
              DO j =1, this%num_part_wall_solid_ghost
                 
                 ip  = this%part_wall_solid_ghost_list(1,j)
                 this%f(1:num_dim,ip) = 0.0_MK       
                 
              END DO
              
           END IF ! num_wall_solid
           
        END IF ! symmetry
        
        
        !----------------------------------------------------
        ! In non-Newtonian fluid, we have to reset velocity
        ! gradient tensor on boundary particles.
        !----------------------------------------------------

        IF ( .NOT. Newtonian ) THEN
           
            
           IF ( symmetry ) THEN
              
              IF ( num_sym > 0 ) THEN
                 
                 DO j =1, this%num_part_sym
                    
                    ip  = this%part_sym_list(1,j)
                    
                    this%vgt(1:num_dim2,ip) = 0.0_MK       
                    
                 END DO
                 
              END IF ! num_sym

              IF ( num_wall_sym > 0 ) THEN
                 
                 DO j = 1, this%num_part_wall_sym
                    
                    ip  = this%part_wall_sym_list(1,j)
                    
                    this%vgt(1:num_dim2,ip) = 0.0_MK
                    
                 END DO

              END IF ! num_wall_sym
              
           END IF ! symmetry
           
        END IF ! non-Newtonian
        
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE particles_reset_boundary_ghost_interaction
      
