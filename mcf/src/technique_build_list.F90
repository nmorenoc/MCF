      SUBROUTINE technique_build_list(this,x,num_part,symmetry,stat_info)
        !----------------------------------------------------
        !  Program      :  technique_build_list
        !----------------------------------------------------
        !
        !  Purpose      : Building cell list using particles'
        !                 positions.
        !
        !  Reference    :
        !
        !  Remark       :
        !
        !  Revisions    : V0.2 08.07.2009, 
        !                 check again the work flow is correct and
        !                 supply with more comments for code.
        !                 Change cell_size(1:num_dim) = 
        !                 this%ghost_size - EPSILON(this%ghost_size)
        !                 to cell_size(1:num_dim) = this%ghost_size
        !
        !                 V0.1 03.03.2009, original version.
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
        ! Module.
        !----------------------------------------------------
        
        USE ppm_module_neighlist
        
        !----------------------------------------------------
        ! Arguments
        !
        ! this       : an object of Marching Class.
        ! stat_info  : return flag of status.
        !----------------------------------------------------
        
        TYPE(Technique), INTENT(INOUT)          :: this
        REAL(MK),DIMENSION(:,:),POINTER         :: x
        INTEGER, INTENT(IN)                     :: num_part
        LOGICAL, INTENT(IN)                     :: symmetry
        INTEGER, INTENT(out)                    :: stat_info
        
        
        !----------------------------------------------------
        ! Local variables
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: num_dim
        REAL(MK), DIMENSION(:),POINTER          :: cell_size
        
        
        !--------------------------------
        ! Initialization of variables.
        !--------------------------------

        stat_info = 0
        stat_info_sub = 0        
        
        NULLIFY(cell_size)
        num_dim = this%num_dim
        ALLOCATE(cell_size(num_dim))
        
        !------------------------------------------
        ! Build neighbor lists.
        !------------------------------------------
        
        IF (this%ghost_size <= 0.0_MK) THEN
           PRINT *, "technique_build_list : ", &
                "ghost size should be positive !"
           stat_info = -1
           GOTO 9999         
        END IF
        
        IF (this%neighbor_list /= 2) THEN
           PRINT *, "technique_build_list : ",&
                "neighbor list besides 2 is not available currently !"
           stat_info = -1
           GOTO 9999           
        END IF

        !------------------------------------------
        ! Boxes size is cutoff in all directions.
        !------------------------------------------
        
        cell_size(1:num_dim) = this%ghost_size
        
        !------------------------------------------
        ! Generate cell lists.
        !------------------------------------------
        
        CALL ppm_neighlist_clist(x,num_part,&
             cell_size,symmetry,this%cell_list,&
             this%num_cell_dim_sub,stat_info_sub)        
        
        !------------------------------------------
        ! Generate cell neighbor lists.
        !------------------------------------------
        
        CALL ppm_neighlist_MkNeighIdx(symmetry,&
             this%inp,this%jnp,this%nnp,stat_info_sub)      
        
        IF (stat_info_sub /= 0) THEN
           PRINT *,'technique_build_list : ',&
                'Failed to build neighborlists!'
           stat_info = -1
           GOTO 9999           
        END IF
        
        
9999    CONTINUE
        
        IF(ASSOCIATED(cell_size)) THEN
           DEALLOCATE(cell_size)
        END IF
        
        RETURN
        
      END SUBROUTINE technique_build_list

      
      SUBROUTINE technique_build_list_c(this,x,num_part,symmetry,ghost_size,stat_info)
        !----------------------------------------------------
        !  Program      :  technique_build_list_c
        !----------------------------------------------------
        !
        !  Purpose      : Building cell list using colloids'
        !                 positions.
        !
        !  Reference    :
        !
        !  Remark       :
        !
        !  Revisions    : V0.1 15.11.2010, original version.
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
        ! Module.
        !----------------------------------------------------
        
        USE ppm_module_neighlist
        
        !----------------------------------------------------
        ! Arguments
        !
        ! this       : an object of Marching Class.
        ! stat_info  : return flag of status.
        !----------------------------------------------------
        
        TYPE(Technique), INTENT(INOUT)          :: this
        REAL(MK),DIMENSION(:,:),POINTER         :: x
        INTEGER, INTENT(IN)                     :: num_part
        LOGICAL, INTENT(IN)                     :: symmetry
        REAL(MK),INTENT(IN)                     :: ghost_size
        INTEGER, INTENT(out)                    :: stat_info
        
        
        !----------------------------------------------------
        ! Local variables
        !----------------------------------------------------
        
        INTEGER                                 :: stat_info_sub
        INTEGER                                 :: num_dim
        REAL(MK), DIMENSION(:),POINTER          :: cell_size
        
        
        !--------------------------------
        ! Initialization of variables.
        !--------------------------------

        stat_info = 0
        stat_info_sub = 0        
        
        NULLIFY(cell_size)
        num_dim = this%num_dim
        ALLOCATE(cell_size(num_dim))
        
        !------------------------------------------
        ! Build neighbor lists.
        !------------------------------------------
        
        IF ( ghost_size <= 0.0_MK ) THEN
           PRINT *, "technique_build_list_c : ", &
                "ghost size should be positive !"
           stat_info = -1
           GOTO 9999         
        END IF
        
        IF (this%neighbor_list_c /= 2) THEN
           PRINT *, "technique_build_list_c : ",&
                "neighbor list besides 2 is not available currently !"
           stat_info = -1
           GOTO 9999           
        END IF
        
        !------------------------------------------
        ! Boxes size is cutoff in all directions.
        !------------------------------------------
        
        cell_size(1:num_dim) = ghost_size
        
        !------------------------------------------
        ! Generate cell lists.
        !------------------------------------------
        
        CALL ppm_neighlist_clist(x,num_part,&
             cell_size,symmetry,this%cell_list_c,&
             this%num_cell_dim_sub_c,stat_info_sub)        
        
        !------------------------------------------
        ! Generate cell neighbor lists.
        !------------------------------------------
        
        CALL ppm_neighlist_MkNeighIdx(symmetry,&
             this%inp_c,this%jnp_c,this%nnp_c,stat_info_sub)      
        
        IF (stat_info_sub /= 0) THEN
           PRINT *,'technique_build_list_c : ',&
                'Failed to build neighborlists!'
           stat_info = -1
           GOTO 9999           
        END IF
        
        
9999    CONTINUE
        
        IF(ASSOCIATED(cell_size)) THEN
           DEALLOCATE(cell_size)
        END IF
        
        RETURN
        
      END SUBROUTINE technique_build_list_c
      
