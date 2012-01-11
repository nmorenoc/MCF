      SUBROUTINE colloid_apply_body_force(this,stat_info)
        !----------------------------------------------------
        ! Subroutine  : colloid_apply_body_force
        !----------------------------------------------------
        ! Purpose     : Apply body force on each colloid.
        !
        ! Routines    :
        !
        ! References  :
        !
        ! Remarks     : body force by default here means
        !               body force per unit mass, i.e.,
        !               acceleration.
        !
        ! Revisions   :  V0.1 28.07.2009, original version
        !
        !----------------------------------------------------
        ! Author      :  Xin Bian 
        ! Contact     :  xin.bian@aer.mw.tum.de
        !
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
       
        
        !----------------------------------------------------
        ! Arguments
        !----------------------------------------------------
        
        TYPE(Colloid), INTENT(OUT)              :: this
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !---------------------------------------------------
        ! Local variables.
        !----------------------------------------------------
        
        INTEGER                                 :: dim,num,i
        
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info = 0
        
        dim = this%num_dim
        num = this%num_colloid
        
        !----------------------------------------------------
        ! Add up body force.
        !----------------------------------------------------
        
        SELECT CASE (this%body_force_type)

        CASE (1)
           
           DO i = 1, dim
              
              this%drag(i,1:num) = this%drag(i,1:num) +&
                   this%body_force(i) * this%m(1:num)
              
           END DO
           
        END SELECT ! body_force_type
        
        
9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE colloid_apply_body_force
      
      
