      INTEGER FUNCTION technique_get_comm(this,stat_info)

        TYPE(Technique),INTENT(IN)              :: this
        INTEGER,INTENT(OUT)                     :: stat_info
        
        stat_info = 0
        
        technique_get_comm = this%comm
        
        RETURN
        
      END FUNCTION technique_get_comm


      INTEGER FUNCTION technique_get_rank(this,stat_info)

        TYPE(Technique),INTENT(IN)              :: this
        INTEGER,INTENT(OUT)                     :: stat_info
        
        
        stat_info = 0
        
        technique_get_rank = this%rank
        
        RETURN
        
      END FUNCTION technique_get_rank
      
      
      INTEGER FUNCTION technique_get_MPI_PREC(this,stat_info)

        TYPE(Technique),INTENT(IN)              :: this
        INTEGER,INTENT(OUT)                     :: stat_info
        
        
        stat_info = 0
        
        technique_get_MPI_PREC = this%MPI_PREC
        
        RETURN
        
      END FUNCTION technique_get_MPI_PREC

      
      REAL(MK) FUNCTION technique_get_ghost_size(this,stat_info)
    
        TYPE(Technique),INTENT(IN)              :: this
        INTEGER,INTENT(OUT)                     :: stat_info
        
        
        stat_info = 0
        
        technique_get_ghost_size = this%ghost_size
    
        RETURN
    
      END FUNCTION technique_get_ghost_size
      
      
      INTEGER FUNCTION technique_get_topo_id(this,stat_info)

        TYPE(Technique),INTENT(IN)              :: this
        INTEGER,INTENT(OUT)                     :: stat_info
        
        
        stat_info = 0
        
        technique_get_topo_id = this%topo_id
        
        RETURN
        
      END FUNCTION technique_get_topo_id
      

      SUBROUTINE technique_get_min_sub(this,d_min,stat_info)
        TYPE(Technique),INTENT(IN)      :: this
        REAL(MK),DIMENSION(:), POINTER  :: d_min
        INTEGER,INTENT(OUT)             :: stat_info
        
        INTEGER                         :: rank, dim
        
        stat_info = 0
        
        IF (ASSOCIATED(d_min)) THEN
           DEALLOCATE(d_min)
        END IF
        
        rank = this%rank
        dim  = this%num_dim
        
        ALLOCATE(d_min(dim))
        d_min(1:dim) = this%min_sub(1:dim,this%proc2sub(rank))
        
        RETURN
      END SUBROUTINE technique_get_min_sub
      

      SUBROUTINE technique_get_max_sub(this,d_max,stat_info)
        TYPE(Technique),INTENT(IN)      :: this
        REAL(MK),DIMENSION(:), POINTER  :: d_max
        INTEGER,INTENT(OUT)             :: stat_info
        
        INTEGER                         :: rank,dim
        
        stat_info = 0
        
        IF (ASSOCIATED(d_max)) THEN
           DEALLOCATE(d_max)
        END IF
        
        rank = this%rank
        dim  = this%num_dim
        
        ALLOCATE(d_max(dim))
        d_max(1:dim) = this%max_sub(1:dim,this%proc2sub(rank))
        
        RETURN
      END SUBROUTINE technique_get_max_sub
      
      
      SUBROUTINE technique_get_sub_bcdef(this,d_bcdef,stat_info)
        TYPE(Technique),INTENT(IN)              :: this
        INTEGER,DIMENSION(:,:), POINTER         :: d_bcdef
        INTEGER,INTENT(OUT)                     :: stat_info
        
        INTEGER                                 :: dim
        
        stat_info = 0
        
        IF (ASSOCIATED(d_bcdef)) THEN
           DEALLOCATE(d_bcdef)
        END IF

        dim = this%num_dim
        
        ALLOCATE(d_bcdef(2*dim,this%num_sub))
        d_bcdef(:,:) = this%sub_bcdef(:,:)
        
        RETURN
      END SUBROUTINE technique_get_sub_bcdef
      
      
      SUBROUTINE technique_get_cell_list(this,&
           d_num_sub,d_num_cell_dim_sub,&
           d_cell_list,d_inp,d_jnp,d_nnp,d_sub_bcdef,stat_info)
        
        TYPE(Technique),INTENT(IN)              :: this
        INTEGER, INTENT(OUT)                    :: d_num_sub
        INTEGER,DIMENSION(:,:),POINTER          :: d_num_cell_dim_sub
        TYPE(ppm_type_ptr_to_clist),DIMENSION(:),POINTER      :: d_cell_list
        INTEGER,DIMENSION(:,:),POINTER          :: d_inp
        INTEGER,DIMENSION(:,:),POINTER          :: d_jnp
        INTEGER, INTENT(OUT)                    :: d_nnp
        INTEGER,DIMENSION(:,:),POINTER          :: d_sub_bcdef
        
        INTEGER,INTENT(OUT)                     :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_num_cell_dim_sub)) THEN
           DEALLOCATE(d_num_cell_dim_sub)
        END IF

        IF(ASSOCIATED(d_sub_bcdef)) THEN
           DEALLOCATE(d_sub_bcdef)
        END IF
        
        IF(ASSOCIATED(d_cell_list)) THEN
           DEALLOCATE(d_cell_list)
        END IF
        
        IF(ASSOCIATED(d_inp)) THEN
           DEALLOCATE(d_inp)
        END IF
        
        IF(ASSOCIATED(d_jnp)) THEN
           DEALLOCATE(d_jnp)
        END IF
        
        ALLOCATE(d_num_cell_dim_sub(SIZE(this%num_cell_dim_sub,1),&
             SIZE(this%num_cell_dim_sub,2)))
        ALLOCATE(d_cell_list(SIZE(this%cell_list,1)))
        ALLOCATE(d_sub_bcdef(SIZE(this%sub_bcdef,1),&
             SIZE(this%sub_bcdef,2)))
        ALLOCATE(d_inp(SIZE(this%inp,1),SIZE(this%inp,2)))
        ALLOCATE(d_jnp(SIZE(this%jnp,1),SIZE(this%jnp,2)))
        
        d_num_sub               = this%num_sub
        d_num_cell_dim_sub(:,:) = this%num_cell_dim_sub(:,:)
        d_cell_list(:)          = this%cell_list(:)
        d_inp(:,:)              = this%inp(:,:)
        d_jnp(:,:)              = this%jnp(:,:)
        d_nnp                   = this%nnp
        
        RETURN
        
      END SUBROUTINE technique_get_cell_list

      
      INTEGER FUNCTION technique_get_nnp(this,stat_info)
        TYPE(Technique),INTENT(IN)              :: this
        INTEGER,INTENT(OUT)                     :: stat_info

        stat_info = 0
        
        technique_get_nnp = this%nnp

        RETURN
      END FUNCTION  Technique_get_nnp
        
        
      SUBROUTINE technique_get_cell_list_c(this,&
           d_num_sub,d_num_cell_dim_sub,&
           d_cell_list,d_inp,d_jnp,d_nnp,d_sub_bcdef,stat_info)
        
        TYPE(Technique),INTENT(IN)              :: this
        INTEGER, INTENT(OUT)                    :: d_num_sub
        INTEGER,DIMENSION(:,:),POINTER          :: d_num_cell_dim_sub
        TYPE(ppm_type_ptr_to_clist),DIMENSION(:),POINTER      :: d_cell_list
        INTEGER,DIMENSION(:,:),POINTER          :: d_inp
        INTEGER,DIMENSION(:,:),POINTER          :: d_jnp
        INTEGER, INTENT(OUT)                    :: d_nnp
        INTEGER,DIMENSION(:,:),POINTER          :: d_sub_bcdef
        
        INTEGER,INTENT(OUT)                     :: stat_info
        
        stat_info = 0
        
        IF(ASSOCIATED(d_num_cell_dim_sub)) THEN
           DEALLOCATE(d_num_cell_dim_sub)
        END IF

        IF(ASSOCIATED(d_sub_bcdef)) THEN
           DEALLOCATE(d_sub_bcdef)
        END IF
        
        IF(ASSOCIATED(d_cell_list)) THEN
           DEALLOCATE(d_cell_list)
        END IF
        
        IF(ASSOCIATED(d_inp)) THEN
           DEALLOCATE(d_inp)
        END IF
        
        IF(ASSOCIATED(d_jnp)) THEN
           DEALLOCATE(d_jnp)
        END IF
        
        ALLOCATE(d_num_cell_dim_sub(SIZE(this%num_cell_dim_sub_c,1),&
             SIZE(this%num_cell_dim_sub_c,2)))
        ALLOCATE(d_cell_list(SIZE(this%cell_list_c,1)))
        ALLOCATE(d_sub_bcdef(SIZE(this%sub_bcdef,1),&
             SIZE(this%sub_bcdef,2)))
        ALLOCATE(d_inp(SIZE(this%inp_c,1),SIZE(this%inp_c,2)))
        ALLOCATE(d_jnp(SIZE(this%jnp_c,1),SIZE(this%jnp_c,2)))
        
        d_num_sub               = this%num_sub
        d_num_cell_dim_sub(:,:) = this%num_cell_dim_sub_c(:,:)
        d_cell_list(:)          = this%cell_list_c(:)
        d_inp(:,:)              = this%inp_c(:,:)
        d_jnp(:,:)              = this%jnp_c(:,:)
        d_nnp                   = this%nnp_c
        
        RETURN
        
      END SUBROUTINE technique_get_cell_list_c

      
      INTEGER FUNCTION technique_get_nnp_c(this,stat_info)
        TYPE(Technique),INTENT(IN)              :: this
        INTEGER,INTENT(OUT)                     :: stat_info

        stat_info = 0
        
        technique_get_nnp_c = this%nnp_c

        RETURN
      END FUNCTION  Technique_get_nnp_c
      
