      SUBROUTINE particles_integrate_ct(this,&
           num,dt,lambda,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_integrate_ct
        !----------------------------------------------------
        !
        ! Purpose     : Integrate the conformation tensor of
        !               particles with required accuracy.
        !
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revision    :  V0.1  29.07.2009, original version.
        !
        !----------------------------------------------------
        ! Author    : Xin Bian
        ! Contact   : xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        !----------------------------------------------------
        ! Arguments
        !
        ! this        : an object of Particles Class.
        ! num         : number of particles updated,
        !               i.e. first num particles in this%x 
        !               are operated.
        ! dt          : time step.
        ! lambda      : coefficient required.
        ! stat_info   : return flag of status.
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)          :: this
        INTEGER, INTENT(IN)                     :: num
        REAL(MK), INTENT(IN)                    :: dt
        REAL(MK), INTENT(IN)                    :: lambda
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------------------------------------
        ! Local variables
	!----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: dim
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        dim = physics_get_num_dim(this%phys,stat_info_sub)
        
        dim = dim**2
        
        !----------------------------------------------------
        ! Update conformation tensor.
        !----------------------------------------------------
        
        this%ct(1:dim,1:num) = &
             this%ct(1:dim,1:num) + &
             this%act(1:dim,1:num) * dt * lambda
        
        
9999    CONTINUE      
        
        RETURN
        
      END SUBROUTINE particles_integrate_ct
      
      
