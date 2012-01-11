      SUBROUTINE  particles_init_partial_exter(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : particles_init_partial_exter
        !----------------------------------------------------
        !
        ! Purpose      :  Create the rest quntities besides
        !                 the ones read from external files.
        !                 Currently no need.e
        !                  
        ! Routines     :
        !
        ! Remarks      :  Create the rest quntities,
        !                 besides position, velocity,
        !                 rho,mass, IDs, for fluid particles 
        !                 in sub-domains of each process.
        !
        !                 In case of Non-Newtonian fluid,
        !                 conformation tensor C or eigenvalues
        !                 and eigenvectors have been
        !                 read already.
        !
        ! References   :
        !
        ! Revisions    :  V0.1 30.07 2009, original version.
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
        ! Arguments
        !----------------------------------------------------
        
        TYPE(Particles), INTENT(INOUT)          :: this
        INTEGER,INTENT(OUT)	                :: stat_info
    	
	!----------------------------------------------------
    	! Local variables
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        
        
        !---------------------------------------------------
    	! Initialization of variables.
    	!---------------------------------------------------
        
        stat_info     = this%num_part_real
        stat_info     = 0
        stat_info_sub = 0
        
        
9999	CONTINUE

        
	RETURN
 
    END SUBROUTINE particles_init_partial_exter
      
      
     
