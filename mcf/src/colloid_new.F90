      SUBROUTINE colloid_init_default(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_init_default
        !----------------------------------------------------
        !
        ! Purpose     : Default construtor of colloid Class.
        !
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revisions   : V0.1 15.07.2009, original version.
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
        ! Arguments
        !----------------------------------------------------

        TYPE(Colloid),INTENT(OUT)       :: this
        INTEGER,INTENT(OUT)             :: stat_info
        
        INTEGER                         :: stat_info_sub
        INTEGER                         :: dim, num
        
        !----------------------------------------------------
        ! For default colloid, there is only one,
        ! i.e., 2D flow around cylinder problem.
        ! Morris et al. J. Comput. Phys. 1997. 
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        NULLIFY(this%tech)
        
        dim = 2
        num = 1

        this%num_dim     = dim
        this%num_colloid = num
        
        this%adapt_t_coef = 1.0_MK
        this%rho         = 1.e3_MK
        this%rho_type    = 0
        this%translate   = .FALSE.
        this%rotate      = .FALSE.
        this%place       = 1
        this%noslip_type = 2
        this%body_force_type = 0
        NULLIFY(this%body_force)
        ALLOCATE(this%body_force(dim))
        this%body_force(:) = 0.0_MK
        
        this%cc_lub_type    = 0
        this%cc_repul_type  = 0
        this%cc_lub_cut_off = 0.0_MK
        this%cc_lub_cut_on  = 0.0_MK
        this%cc_repul_cut_off = 0.0_MK
        this%cc_repul_cut_on  = 0.0_MK
        this%cc_repul_F0      = 0.0_MK
        
        this%cw_lub_type    = 0
        this%cw_repul_type  = 0
        this%cw_lub_cut_off = 0.0_MK
        this%cw_lub_cut_on  = 0.0_MK
        this%cw_repul_cut_off = 0.0_MK
        this%cw_repul_cut_on  = 0.0_MK        
        this%cw_repul_F0      = 0.0_MK
        
        this%h                = 0.0_MK
        this%dt_f             = -1.0_MK
        
        NULLIFY(this%shape)
        ALLOCATE(this%shape(num))
        this%shape(:) = 1
        
        NULLIFY(this%radius)
        ALLOCATE(this%radius(dim,num))
        this%radius(1,:) = 0.02_MK
        
        NULLIFY(this%freq)
        ALLOCATE(this%freq(num))
        this%freq(1) = 0
        
        NULLIFY(this%m)
        ALLOCATE(this%m(num))
        this%m(:) = 0.0_MK

        NULLIFY(this%mmi)
        ALLOCATE(this%mmi(3,num))
        this%mmi(:,:) = 0.0_MK
    
        NULLIFY(this%x)
        ALLOCATE(this%x(dim,num))
        this%x(:,1) = 0.05_MK
        
        NULLIFY(this%v)
        ALLOCATE(this%v(dim,num))
        this%v(:,:) = 0.0_MK
        
#if __DRAG_PART
        NULLIFY(this%drag_lub)
        ALLOCATE(this%drag_lub(dim,num))
        this%drag_lub(:,:) = 0.0_MK
        
        NULLIFY(this%drag_repul)
        ALLOCATE(this%drag_repul(dim,num))
        this%drag_repul(:,:) = 0.0_MK     
#endif
    
        NULLIFY(this%drag)
        ALLOCATE(this%drag(dim,num))
        this%drag(:,:) = 0.0_MK
        
        NULLIFY(this%rot_vector)
        ALLOCATE(this%rot_vector(4,num))
        this%rot_vector(:,:) = 0.0_MK
    
        NULLIFY(this%acc_vector)
        ALLOCATE(this%acc_vector(4,num))
        this%acc_vector(:,:) = 0.0_MK
   
        NULLIFY(this%rot_matrix)
        ALLOCATE(this%rot_matrix(3,3,num))
        this%rot_matrix(:,:,:) = 0.0_MK
        this%rot_matrix(1,1,:) = 1.0_MK
        this%rot_matrix(2,2,:) = 1.0_MK
        this%rot_matrix(3,3,:) = 1.0_MK
        
        NULLIFY(this%acc_matrix)
        ALLOCATE(this%acc_matrix(3,3,num))
        this%acc_matrix(:,:,:) = this%rot_matrix(:,:,:)
        
        NULLIFY(this%theta)
        ALLOCATE(this%theta(3,num))
        this%theta(:,:) = 0.0_MK
    
        NULLIFY(this%omega)
        ALLOCATE(this%omega(3,num))
        this%omega(:,:) = 0.0_MK
        
        NULLIFY(this%torque)
        ALLOCATE(this%torque(3,num))
        this%torque(:,:) = 0.0_MK
        
        NULLIFY(this%num_physical_part)
        ALLOCATE(this%num_physical_part(num))
        this%num_physical_part(:)  = 0
        
        NULLIFY(this%num_numerical_part)
        ALLOCATE(this%num_numerical_part(num))
        this%num_numerical_part(:)  = 0
         
        NULLIFY(this%f)
        ALLOCATE(this%f(dim,num))
        this%f(:,:) = 0.0_MK
        
        this%fa_min = 0.0_MK
        this%fa_max = 0.0_MK

        NULLIFY(this%alpha)
        ALLOCATE(this%alpha(3,num))
        this%alpha(:,:) = 0.0_MK
        
        NULLIFY(this%k_energy)
        ALLOCATE(this%k_energy(num))
        this%k_energy(:) = 0.0_MK

        NULLIFY(this%mom)
        ALLOCATE(this%mom(dim,num))
        this%mom(:,:) = 0.0_MK
        
        NULLIFY(this%mom_tot)
        ALLOCATE(this%mom_tot(dim))
        this%mom_tot(:) = 0.0_MK
        
        this%num_physical_part_tot = 0
        this%num_numerical_part_tot = 0
        
        NULLIFY(this%min_phys)
        ALLOCATE(this%min_phys(dim))
        this%min_phys(:) = 0.0_MK
        
        NULLIFY(this%max_phys)
        ALLOCATE(this%max_phys(dim))
        this%max_phys(:) = 0.0_MK
        
        NULLIFY(this%min_phys_t)
        ALLOCATE(this%min_phys_t(dim))
        this%min_phys_t(:) = 0.0_MK
        
        NULLIFY(this%max_phys_t)
        ALLOCATE(this%max_phys_t(dim))
        this%max_phys_t(:) = 0.0_MK
        
        NULLIFY(this%bcdef)
        ALLOCATE(this%bcdef(2*dim))
        this%bcdef(:)  = ppm_param_bcdef_periodic
        NULLIFY(this%boundary)

        this%cut_off = 0.0_MK
        this%dout    = 0.0_MK
        this%din     = 0.0_MK
        this%eta     = 0.0_MK

        !----------------------------------------------------
        ! Images of center: position and velocity.
        !----------------------------------------------------
        
        this%num_image = 0
        NULLIFY(this%x_image)
        NULLIFY(this%v_image)
        
        CALL tool_new(this%tool,stat_info_sub)
        
        RETURN          
        
      END SUBROUTINE colloid_init_default
      
      
      SUBROUTINE colloid_init(this,d_num_dim,d_num_colloid,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_init
        !----------------------------------------------------
        !
        ! Purpose     : Construtor of colloid Class.
        !
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revisions   : V0.1 15.07.2009, original version.
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
        
        TYPE(Colloid),INTENT(OUT)       :: this
        INTEGER, INTENT(IN)             :: d_num_dim
        INTEGER, INTENT(IN)             :: d_num_colloid
        INTEGER,INTENT(OUT)             :: stat_info
        
        INTEGER                         :: stat_info_sub

        !----------------------------------------------------
        ! For non default colloid(s).
        !----------------------------------------------------

        stat_info     = 0
        stat_info_sub = 0

        NULLIFY(this%tech)

        this%num_dim     = d_num_dim
        this%num_colloid = d_num_colloid
        
        this%adapt_t_coef= 1.0_MK
        this%rho         = 1.e3_MK
        this%rho_type    = 0
        this%translate   = .FALSE.
        this%rotate      = .FALSE.
        this%place       = 1
        this%noslip_type = 1
        this%body_force_type = 0
        NULLIFY(this%body_force)
        ALLOCATE(this%body_force(d_num_dim))
        this%body_force(:) = 0.0_MK
      
        this%cc_lub_type = 0
        this%cc_lub_cut_off = 0.0_MK
        this%cc_lub_cut_on  = 0.0_MK

        this%cc_repul_type = 0
        this%cc_repul_cut_off = 0.0_MK
        this%cc_repul_cut_on  = 0.0_MK
        this%cc_repul_F0      = 0.0_MK
        
        this%cw_lub_type = 0
        this%cw_lub_cut_off = 0.0_MK
        this%cw_lub_cut_on  = 0.0_MK
        
        this%cw_repul_type = 0
        this%cw_repul_cut_off = 0.0_MK
        this%cw_repul_cut_on  = 0.0_MK
        this%cw_repul_F0      = 0.0_MK
        
        this%h                = 0.0_MK
        this%dt_f             = -1.0_MK
        
        !----------------------------------------------------
        ! Default shape is sphereical.
        !----------------------------------------------------

        NULLIFY(this%shape)
        ALLOCATE(this%shape(d_num_colloid))
        this%shape(:) = 1
        
        NULLIFY(this%radius)
        ALLOCATE(this%radius(d_num_dim,d_num_colloid))
        this%radius(:,:) = 0.0_MK
        
        NULLIFY(this%freq)
        ALLOCATE(this%freq(d_num_colloid))
        this%freq(:) = 0
        
        NULLIFY(this%m)
        ALLOCATE(this%m(d_num_colloid))
        this%m(:) = 0.0_MK

        NULLIFY(this%mmi)
        ALLOCATE(this%mmi(3,d_num_colloid))
        this%mmi(:,:) = 0.0_MK
        
        NULLIFY(this%x)
        ALLOCATE(this%x(d_num_dim,d_num_colloid))
        this%x(:,:) = 0.0_MK
        
        NULLIFY(this%v)
        ALLOCATE(this%v(d_num_dim,d_num_colloid))
        this%v(:,:) = 0.0_MK
        
#if __DRAG_PART
        NULLIFY(this%drag_lub)
        ALLOCATE(this%drag_lub(d_num_dim,d_num_colloid))
        this%drag_lub(:,:) = 0.0_MK
        
        NULLIFY(this%drag_repul)
        ALLOCATE(this%drag_repul(d_num_dim,d_num_colloid))
        this%drag_repul(:,:) = 0.0_MK
#endif

        NULLIFY(this%drag)
        ALLOCATE(this%drag(d_num_dim,d_num_colloid))
        this%drag(:,:) = 0.0_MK
        
        NULLIFY(this%rot_vector)
        ALLOCATE(this%rot_vector(4,d_num_colloid))
        this%rot_vector(:,:) = 0.0_MK
    
        NULLIFY(this%acc_vector)
        ALLOCATE(this%acc_vector(4,d_num_colloid))
        this%acc_vector(:,:) = 0.0_MK
    
        NULLIFY(this%rot_matrix)
        ALLOCATE(this%rot_matrix(3,3,d_num_colloid))
        this%rot_matrix(:,:,:) = 0.0_MK
        
        NULLIFY(this%acc_matrix)
        ALLOCATE(this%acc_matrix(3,3,d_num_colloid))
        this%acc_matrix(:,:,:) = 0.0_MK
        
        NULLIFY(this%theta)
        ALLOCATE(this%theta(3,d_num_colloid))
        this%theta(:,:) = 0.0_MK
    
        NULLIFY(this%omega)
        ALLOCATE(this%omega(3,d_num_colloid))
        this%omega(:,:) = 0.0_MK
        
        NULLIFY(this%torque)
        ALLOCATE(this%torque(3,d_num_colloid))
        this%torque(:,:) = 0.0_MK
        
        NULLIFY(this%num_physical_part)
        ALLOCATE(this%num_physical_part(d_num_colloid))
        this%num_physical_part(:)  = 0
        
        NULLIFY(this%num_numerical_part)
        ALLOCATE(this%num_numerical_part(d_num_colloid))
        this%num_numerical_part(:)  = 0
         

        !----------------------------------------------------
        ! Derived quantities.                               
        !----------------------------------------------------
        
        NULLIFY(this%f)
        ALLOCATE(this%f(d_num_dim,d_num_colloid))
        this%f(:,:) = 0.0_MK
        
        this%fa_min = 0.0_MK
        this%fa_max = 0.0_MK
        
        NULLIFY(this%alpha)
        ALLOCATE(this%alpha(3,d_num_colloid))
        this%alpha(:,:) = 0.0_MK
        
        NULLIFY(this%k_energy)
        ALLOCATE(this%k_energy(d_num_colloid))
        this%k_energy(:) = 0.0_MK

        NULLIFY(this%mom)
        ALLOCATE(this%mom(d_num_dim,d_num_colloid))
        this%mom(:,:) = 0.0_MK
        
        NULLIFY(this%mom_tot)
        ALLOCATE(this%mom_tot(d_num_dim))
        this%mom_tot(:) = 0.0_MK
        
        this%num_physical_part_tot = 0
        this%num_numerical_part_tot = 0
        
        !----------------------------------------------------
        ! Physics parameters, boundaries.
        !----------------------------------------------------
        
        NULLIFY(this%min_phys)
        ALLOCATE(this%min_phys(d_num_dim))
        this%min_phys(:) = 0.0_MK
        
        NULLIFY(this%max_phys)
        ALLOCATE(this%max_phys(d_num_dim))
        this%max_phys(:) = 0.0_MK
        
        NULLIFY(this%min_phys_t)
        ALLOCATE(this%min_phys_t(d_num_dim))
        this%min_phys_t(:) = 0.0_MK
        
        NULLIFY(this%max_phys_t)
        ALLOCATE(this%max_phys_t(d_num_dim))
        this%max_phys_t(:) = 0.0_MK
        
        NULLIFY(this%bcdef)
        ALLOCATE(this%bcdef(2*d_num_dim))
        this%bcdef(:)  = ppm_param_bcdef_periodic
        NULLIFY(this%boundary)

        this%cut_off = 0.0_MK
        this%dout    = 0.0_MK
        this%din     = 0.0_MK
        this%eta     = 0.0_MK

        !----------------------------------------------------
        ! Images of center: position and velocity.
        !----------------------------------------------------
        
        this%num_image = 0
        NULLIFY(this%x_image)
        NULLIFY(this%v_image)

        CALL tool_new(this%tool,stat_info_sub)
                
        RETURN          
        
      END SUBROUTINE colloid_init
      
      
      SUBROUTINE colloid_display_parameters(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_display_parameters
        !----------------------------------------------------
        !
        ! Purpose     : To display colloid paramters.
        !
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revisions   : V0.1 15.07.2009, original version.
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
        
        TYPE(Colloid),INTENT(IN)        :: this
        INTEGER,INTENT(OUT)             :: stat_info
        
        INTEGER                         :: dim
        INTEGER                         :: i
        
        stat_info = 0
        dim = this%num_dim
        
        PRINT *, '---***************Start***************---'
        PRINT *, '     Colloid parameters'
        PRINT *, '---***********************************---'
        
        PRINT *, "num_colloid        : ", this%num_colloid
        PRINT *, "adapt_t_coef       : ", this%adapt_t_coef
        PRINT *, "rho                : ", this%rho
        PRINT *, "rho type           : ", this%rho_type
        PRINT *, "translate          : ", this%translate
        PRINT *, "rotate             : ", this%rotate
        PRINT *, "particle placement : ", this%place
        PRINT *, "no slip            : ", this%noslip_type
        PRINT *, "body force type    : ", this%body_force_type
        PRINT *, "body force         : ", this%body_force(1:dim)
        
        PRINT *, "cc_lub_type        : ", this%cc_lub_type
        PRINT *, "cc_lub_cut_off     : ", this%cc_lub_cut_off
        PRINT *, "cc_lub_cut_on      : ", this%cc_lub_cut_on

        PRINT *, "cc_repul_type      : ", this%cc_repul_type
        PRINT *, "cc_repul_cut_off   : ", this%cc_repul_cut_off
        PRINT *, "cc_repul_cut_on    : ", this%cc_repul_cut_on
        PRINT *, "cc_repul_F0        : ", this%cc_repul_F0
        
        PRINT *, "cw_lub_type        : ", this%cw_lub_type
        PRINT *, "cw_lub_cut_off     : ", this%cw_lub_cut_off
        PRINT *, "cw_lub_cut_on      : ", this%cw_lub_cut_on

        PRINT *, "cw_repul_type      : ", this%cw_repul_type
        PRINT *, "cw_repul_cut_off   : ", this%cw_repul_cut_off
        PRINT *, "cw_repul_cut_on    : ", this%cw_repul_cut_on
        PRINT *, "cw_repul_F0        : ", this%cw_repul_F0
        
        PRINT *, "num_image          : ", this%num_image
        
        
        DO i = 1, this%num_colloid
           
           PRINT *, '   ----------------------------------'
           PRINT *, "colloid index      : ", i
           PRINT *, "shape              : ", this%shape(i)
           PRINT *, "radius             : ", this%radius(1:dim,i)
           PRINT *, "freq               : ", this%freq(i)
           PRINT *, "m                  : ", this%m(i)
           PRINT *, "mmi                : ", this%mmi(1:3,i)
           PRINT *, "x                  : ", this%x(1:dim,i)
           PRINT *, "v                  : ", this%v(1:dim,i)
           PRINT *, "rotation vector    : ", this%acc_vector(1:4,i)
           !PRINT *, "theta              : ", this%theta(1:3,i)
           PRINT *, "omega              : ", this%omega(1:3,i)
           PRINT *, "num_physical_part  : ",&
                this%num_physical_part(i)
           PRINT *, "num_numerical_part : ", &
                this%num_numerical_part(i)
           PRINT *, '   ----------------------------------'
           
        END DO
        
        PRINT *, '---****************End****************---'
        

9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE colloid_display_parameters
      
      
