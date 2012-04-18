      SUBROUTINE colloid_adjust_parameters(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_adjust_parameters
        !----------------------------------------------------
        !
        ! Purpose     : Adjust colloid parameters resonably.
        !
        ! Reference   :
        !
        ! Remark      :
        !
        ! Revisions   :  V0.1 21.05.2010, original version.
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
        
        TYPE(colloid), INTENT(INOUT)    :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info = 0
        
        !----------------------------------------------------
        ! If colloid doesn't translate or rotate,
        ! we set its velocity to zero.
        !----------------------------------------------------
        
        IF ( .NOT. this%translate ) THEN
              
           this%v(:,:,:) = 0.0_MK
           
        END IF
           
        IF ( .NOT. this%rotate ) THEN
              
           this%omega(:,:,:) = 0.0_MK
           
        END IF
        
        !----------------------------------------------------
        ! Set h to threshold gap, i.e., cc_lub_cut_on
        !----------------------------------------------------
        
        !this%h = this%cc_lub_cut_on
        
        !----------------------------------------------------
        ! Set h to threshold gap, i.e., cc_lub_cut_off
        !----------------------------------------------------
        
        IF ( this%cc_lub_type > 0 .AND. &
             this%cc_lub_cut_off > 0.0_MK ) THEN
           
           this%h = this%cc_lub_cut_off / 3.0_MK
           
        END IF
        
        IF ( this%cc_repul_type > 0 .AND. &
             this%cc_repul_cut_off < this%cc_lub_cut_off) THEN
           
           this%h = this%cc_repul_cut_off / 3.0_MK
           
        END IF
        
        RETURN
        
      END SUBROUTINE  colloid_adjust_parameters
      
