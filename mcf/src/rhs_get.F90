      INTEGER FUNCTION rhs_get_rhs_density_type(this,stat_info)

        TYPE(Rhs), INTENT(INOUT)        :: this
        INTEGER, INTENT(OUT)            :: stat_info

        stat_info = 0
        
        rhs_get_rhs_density_type = this%rhs_density_type
        
        RETURN
        
      END FUNCTION  rhs_get_rhs_density_type
      
      
      INTEGER FUNCTION rhs_get_rhs_force_type(this,stat_info)

        TYPE(Rhs), INTENT(INOUT)        :: this
        INTEGER, INTENT(OUT)            :: stat_info

        stat_info = 0
        
        rhs_get_rhs_force_type = this%rhs_force_type
        
        RETURN
        
      END FUNCTION  rhs_get_rhs_force_type
      

      REAL(MK) FUNCTION rhs_get_dt(this,stat_info)
        
        TYPE(Rhs), INTENT(INOUT)        :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        rhs_get_dt = this%dt
        
        RETURN
        
      END FUNCTION  rhs_get_dt
      
      
      REAL(MK) FUNCTION rhs_get_kt(this,stat_info)
        
        TYPE(Rhs), INTENT(INOUT)        :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        
        rhs_get_kt = this%kt
        
        RETURN
        
      END FUNCTION  rhs_get_kt
      
      
      SUBROUTINE rhs_get_random(this,d_random,stat_info)

        TYPE(Rhs), INTENT(INOUT)        :: this
        TYPE(Random), POINTER           :: d_random
        INTEGER, INTENT(OUT)            :: stat_info

        stat_info = 0
        
        NULLIFY(d_random)
        d_random => this%random
        
        RETURN
        
      END SUBROUTINE  rhs_get_random
