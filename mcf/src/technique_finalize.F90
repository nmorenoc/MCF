      SUBROUTINE technique_finalize(this,error_info,stat_info)
        !--------------------------------------------------------------
        ! Subroutine  :  technique_finalize
        !--------------------------------------------------------------
        !
        ! Purpose     : 
        !
        ! Input       :
        !
        ! Output      : stat_info   (I) return status. 0 on success.
        !
        ! Routines    : ppm_finalize
        !
        ! Remarks     :
        !
        ! References  :
        !
        ! Revisions   : V0.1 01.02.2009
        ! 
        !----------------------------------------------------
        ! Author       : Xin Bian
        ! Contact      : xin.bian@aer.mw.tum.de
        ! Dr. Marco Ellero's Emmy Noether Group,
        ! Prof. Dr. N. Adams' Chair of Aerodynamics,
        ! Faculty of Mechanical Engineering,
        ! Technische Universitaet Muenchen, Germany.
        !----------------------------------------------------
        
        TYPE(Technique), INTENT(INOUT)  :: this
        INTEGER, INTENT(IN)             :: error_info
        INTEGER, INTENT(OUT)            :: stat_info
        
        INTEGER                         :: stat_info_sub
        
        
        stat_info     = 0
        stat_info_sub = 0
        
        IF (ASSOCIATED(this%min_phys_t)) THEN
           DEALLOCATE(this%min_phys_t)
        END IF

        IF (ASSOCIATED(this%max_phys_t)) THEN
           DEALLOCATE(this%max_phys_t)
        END IF
        
        IF (ASSOCIATED(this%bcdef)) THEN
           DEALLOCATE(this%bcdef)
        END IF
        
        IF (ASSOCIATED(this%min_sub)) THEN
           DEALLOCATE(this%min_sub)
        END IF
       
        IF (ASSOCIATED(this%max_sub)) THEN
           DEALLOCATE(this%max_sub)
        END IF

        IF (ASSOCIATED(this%sub_cost)) THEN
           DEALLOCATE(this%sub_cost)
        END IF
        
        IF (ASSOCIATED(this%sub2proc)) THEN
           DEALLOCATE(this%sub2proc)
        END IF

        IF (ASSOCIATED(this%proc2sub)) THEN
           DEALLOCATE(this%proc2sub)
        END IF
        
        IF (ASSOCIATED(this%sub_list)) THEN
           DEALLOCATE(this%sub_list)
        END IF
       
        IF (ASSOCIATED(this%sub_bcdef)) THEN
           DEALLOCATE(this%sub_bcdef)
        END IF
        
        IF (ASSOCIATED(this%num_cell_dim_sub)) THEN
           DEALLOCATE(this%num_cell_dim_sub)
        END IF

        IF (ASSOCIATED(this%cell_list)) THEN
           DEALLOCATE(this%cell_list)
        END IF
        
        IF (ASSOCIATED(this%inp)) THEN
           DEALLOCATE(this%inp)
        END IF
       
        IF (ASSOCIATED(this%jnp)) THEN
           DEALLOCATE(this%jnp)
        END IF
        
        IF (ASSOCIATED(this%num_cell_dim_sub_c)) THEN
           DEALLOCATE(this%num_cell_dim_sub_c)
        END IF

        IF (ASSOCIATED(this%cell_list_c)) THEN
           DEALLOCATE(this%cell_list_c)
        END IF
        
        IF (ASSOCIATED(this%inp_c)) THEN
           DEALLOCATE(this%inp_c)
        END IF
       
        IF (ASSOCIATED(this%jnp_c)) THEN
           DEALLOCATE(this%jnp_c)
        END IF

        PRINT *, "technique_finalize : ", "Finished!"

        CALL technique_finalize_parallelization(this,error_info,stat_info_sub)

        RETURN
        
      END SUBROUTINE technique_finalize
      
      
      SUBROUTINE technique_finalize_parallelization(this,d_error_info,stat_info)
        !--------------------------------------------------------------
        ! Subroutine  :  technique_finalize_parallelization
        !--------------------------------------------------------------
        !
        ! Purpose     : 
        !
        ! Input       :
        !
        ! Output      : stat_info   (I) return status. 0 on success.
        !
        ! Routines    : ppm_finalize
        !
        ! Remarks     :
        !
        ! References  :
        !
        ! Revisions   : V0.1 01.02.2009
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
        
        TYPE(Technique),INTENT(IN)      :: this
        INTEGER, INTENT(IN)             :: d_error_info
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Local variables 
        !----------------------------------------------------
        
        INTEGER                   :: stat_info_sub
        INTEGER                   :: error_info
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        error_info    = d_error_info
        stat_info     = 0
        stat_info_sub = 0
        
	
        !----------------------------------------------------
        ! Finalize PPM
        !----------------------------------------------------
        
        CALL technique_finalize_ppm(this,stat_info_sub)
        
        IF ( error_info == 0 ) THEN
           error_info = stat_info_sub
        END IF
        
        !--------------------------------
        ! Finalize MPI
        !--------------------------------
        
        IF ( error_info == 0 ) THEN
           
           CALL MPI_Finalize(stat_info_sub)
           
           IF (stat_info_sub /= 0) THEN
              PRINT *, "technique_finalize_parallelization : ", &
                   "FAILED TO FINALIZE MPI!"
              stat_info = -1
              GOTO 9999
           END IF
           
        ELSE
           
           !-----------------------------
           ! Abort MPI
           !-----------------------------
           
           CALL MPI_Abort(this%comm,error_info,stat_info_sub)
           
           IF (stat_info_sub /= 0) THEN
              PRINT * , "technique_finalize_parallelization : ", &
                   "FAILED TO ABORT MPI!"
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
    
        !----------------------------------------------------
        ! Return 
        !----------------------------------------------------
        
9999    CONTINUE
        
        
        PRINT *, "technique_finalize_parallelization : ", "Finished!"

        RETURN
        
      END SUBROUTINE technique_finalize_parallelization



      SUBROUTINE technique_finalize_ppm(this,stat_info)
        
        !----------------------------------------------------
        !  Subroutine   :    technique_finalize_pppm
        !----------------------------------------------------
        !
        !  Purpose      : 
        !
        !  Input        :     
        !
        !  Input/output : 
        !
        !  Output       : stat_info (I) return status. 0 on success.
        !
        !  Routines     : 
        !
        !  Remarks      :
        !
        !  References   :
        !
        !  Revisions    : V0.1 01.02.2009
        ! 
        !----------------------------------------------------
        !  Author       : Xin Bian
        !  Contact      : xin.bian@aer.mw.tum.de
        !----------------------------------------------------
        
        
        !----------------------------------------------------
        !  Modules 
        !----------------------------------------------------
        USE ppm_module_user_util        
        
        !----------------------------------------------------
        !  Includes 
        !----------------------------------------------------
                
        !----------------------------------------------------
        !  Arguments     
        !----------------------------------------------------
        
        TYPE(Technique),INTENT(IN)      :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        !  Local variables 
        !----------------------------------------------------
        
        INTEGER                   :: stat_info_sub
        
        !----------------------------------------------------
        !  Initialization of variables
        !----------------------------------------------------
        
        stat_info = this%num_dim
        stat_info = 0
        stat_info_sub = 0
        
        !----------------------------------------------------
        !  Finalize PPM
        !-----------------------------------------------------
        
        CALL ppm_finalize(stat_info_sub)
        
        IF( stat_info_sub /= 0 ) THEN
           
           PRINT *, "technique_finalize_ppm : ", " error !"
           stat_info=stat_info_sub
        END IF               
        
        
        RETURN
        
      END SUBROUTINE technique_finalize_ppm
      
      
