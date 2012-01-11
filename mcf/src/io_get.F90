!--------------------------------------------------
! Subroutine  :  io_get_*
!--------------------------------------------------
!
! Purpose     : Get routines of Class IO.
!
! Reference   :
!
! Remark      :
!
! Revisions   : V0.1 10.03.2010, original version.
!
!--------------------------------------------------
! Author      : Xin Bian
! Contact     : xin.bian@aer.mw.tum.de
!
! Dr. Marco Ellero's Emmy Noether Group,
! Prof. Dr. N. Adams' Chair of Aerodynamics,
! Faculty of Mechanical Engineering,
! Technische Universitaet Muenchen, Germany.
!-------------------------------------------------  

      INTEGER FUNCTION io_get_output_particles_relax_freq_step(this,&
           stat_info)
        
        !----------------------------------------------------
        ! Arguments.
        !----------------------------------------------------
        
        TYPE(IO), INTENT(IN)            :: this             
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0

        io_get_output_particles_relax_freq_step = &
             this%output_particles_relax_freq_step
        
        RETURN
        
      END FUNCTION io_get_output_particles_relax_freq_step
      
      
      INTEGER FUNCTION io_get_statistic_relax_freq_step(this,&
           stat_info)
        
        !----------------------------------------------------
        ! Arguments.
        !----------------------------------------------------
        
        TYPE(IO), INTENT(IN)            :: this             
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0

        io_get_statistic_relax_freq_step = &
             this%statistic_relax_freq_step
        
        RETURN
        
      END FUNCTION io_get_statistic_relax_freq_step
