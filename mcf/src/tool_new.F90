      SUBROUTINE tool_init(this,stat_info)

        TYPE(Tool), INTENT(OUT)         :: this
        INTEGER, INTENT(OUT)            :: stat_info
        
        stat_info = 0
        this%flag = 0

        RETURN

      ENd SUBROUTINE tool_init
