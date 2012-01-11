        SUBROUTINE debug_set_flag(this,d_flag,stat_info)
          
          TYPE(Debug), INTENT(INOUT)       :: this
          INTEGER, INTENT(IN)           :: d_flag
          INTEGER, INTENT(OUT)          :: stat_info


          stat_info = 0
          this%flag = d_flag
          

          RETURN           

        END SUBROUTINE debug_set_flag
