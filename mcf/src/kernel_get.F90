      INTEGER FUNCTION kernel_get_num_dim(this,stat_info)

        TYPE(Kernel), INTENT(IN)        :: this
        INTEGER, INTENT(OUT)            :: stat_info

        stat_info = 0

        kernel_get_num_dim = this%num_dim

        RETURN

      END FUNCTION kernel_get_num_dim


     INTEGER FUNCTION kernel_get_kernel_type(this,stat_info)

        TYPE(Kernel), INTENT(IN)        :: this
        INTEGER, INTENT(OUT)            :: stat_info

        stat_info = 0

        kernel_get_kernel_type = this%kernel_type

        RETURN

      END FUNCTION kernel_get_kernel_type


      REAL(MK) FUNCTION kernel_get_h(this,stat_info)

        TYPE(Kernel), INTENT(IN)        :: this
        INTEGER, INTENT(OUT)            :: stat_info

        stat_info = 0

        kernel_get_h = this%h

        RETURN

      END FUNCTION kernel_get_h


      REAL(MK) FUNCTION kernel_get_cut_off(this,stat_info)

        TYPE(Kernel), INTENT(IN)        :: this
        INTEGER, INTENT(OUT)            :: stat_info

        stat_info = 0

        kernel_get_cut_off = this%cut_off

        RETURN

      END FUNCTION kernel_get_cut_off
