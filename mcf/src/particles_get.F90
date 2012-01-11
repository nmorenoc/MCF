!------------------------------------------------------------
! All the public "get" subroutines of Class Particles,
! which return member variables of Particles.
!
! Reference   :
!
! Remark      :
!
! Revisions   : V0.1 18.08.2010
!               Added the subroutine particles_get_eval.
!               (Adolfo)
!
!               V0.1 03.03.2009, original version.
!
!------------------------------------------------------------
! Author       : Xin Bian
! Contact      : xin.bian@aer.mw.tum.de
!
! Dr. Marco Ellero's Emmy Noether Group,
! Prof. Dr. N. Adams' Chair of Aerodynamics,
! Faculty of Mechanical Engineering,
! Technische Universitaet Muenchen, Germany.
!------------------------------------------------------------

      SUBROUTINE particles_get_ctrl(this,d_ctrl,stat_info)

        TYPE(Particles), INTENT(IN)             :: this
        TYPE(Control), POINTER                  :: d_ctrl
        INTEGER, INTENT(OUT)                    :: stat_info

        stat_info = 0
        NULLIFY(d_ctrl)
        d_ctrl => this%ctrl

        RETURN
        
      END SUBROUTINE particles_get_ctrl
      
      
      SUBROUTINE particles_get_phys(this,d_phys,stat_info)
        
        TYPE(Particles), INTENT(IN)             :: this
        TYPE(Physics), POINTER                  :: d_phys
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = 0
        NULLIFY(d_phys)
        d_phys => this%phys

        RETURN
        
      END SUBROUTINE particles_get_phys
      
      
      SUBROUTINE particles_get_rhs(this,d_rhs,stat_info)
        
        TYPE(Particles), INTENT(IN)             :: this
        TYPE(Rhs), POINTER                      :: d_rhs
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = 0
        NULLIFY(d_rhs)
        d_rhs => this%rhs

        RETURN
        
      END SUBROUTINE particles_get_rhs
      

      SUBROUTINE particles_get_stateEquation(this,d_sE,stat_info)
        
        TYPE(Particles), INTENT(IN)             :: this
        TYPE(StateEquation), POINTER            :: d_sE
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = 0
        NULLIFY(d_sE)
        d_sE => this%stateEquation
        
        RETURN
        
      END SUBROUTINE particles_get_stateEquation

      
      SUBROUTINE particles_get_kernel(this,d_kern,stat_info)
        
        TYPE(Particles), INTENT(IN)             :: this
        TYPE(Kernel), POINTER                  :: d_kern
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = 0
        NULLIFY(d_kern)
        d_kern => this%kern
        
        RETURN
        
      END SUBROUTINE particles_get_kernel
      

      SUBROUTINE particles_get_tech(this,d_tech,stat_info)
         
        TYPE(Particles), INTENT(IN)             :: this
        TYPE(Technique), POINTER                :: d_tech
        INTEGER, INTENT(OUT)                    :: stat_info
         
         
        stat_info = 0
        NULLIFY(d_tech)        
        d_tech => this%tech
         
        RETURN
         
      END SUBROUTINE particles_get_tech
    
      
      INTEGER FUNCTION particles_get_num_dim(this,stat_info)
        
        TYPE(Particles), INTENT(IN)             :: this
        INTEGER, INTENT(OUT)                    :: stat_info
        
        stat_info = 0
        
        particles_get_num_dim = this%num_dim

        RETURN
        
      END FUNCTION particles_get_num_dim
      
      
      SUBROUTINE particles_get_x(this,x,num,stat_info)

        TYPE(Particles),INTENT(IN)              :: this
        REAL(MK),DIMENSION(:,:),POINTER         :: x
        INTEGER, INTENT(IN)                     :: num
        INTEGER,INTENT(OUT)                     :: stat_info
        
        !----------------------
        ! Local variables
        !----------------------
        
        INTEGER                         :: dim
        
        stat_info = 0
        
        IF ( num > this%num_part_all ) THEN
           PRINT *, "particles_get_x : ", &
                "Required particles is more than existed !"
           stat_info  = -1
           GOTO 9999
        END IF
        
        IF(ASSOCIATED(x)) THEN 
           DEALLOCATE(x)
        END IF
        
        IF ( num > 0 ) THEN
           
           dim = SIZE(this%x,1)
           
           ALLOCATE(x(dim,num))
           
           x(1:dim,1:num) = this%x(1:dim,1:num)
           
        END IF
        
9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE particles_get_x
      
      
      SUBROUTINE particles_get_v(this,v,num,stat_info)
          
        TYPE(Particles),INTENT(IN)              :: this
        REAL(MK),DIMENSION(:,:),POINTER         :: v
        INTEGER, INTENT(IN)                     :: num
        INTEGER, INTENT(OUT)                    :: stat_info
        
        !----------------------
        ! Local variables
        !----------------------
        
        INTEGER                         :: dim

        stat_info = 0
        
        IF ( num > this%num_part_all ) THEN
           PRINT *, "particles_get_v : ", &
                "Required particles is more than existed !"
           stat_info  = -1
           GOTO 9999
        END IF
        
        IF(ASSOCIATED(v)) THEN 
           DEALLOCATE(v)
        END IF
        
        IF ( num > 0 ) THEN
           
           dim = SIZE(this%x,1)
           
           ALLOCATE(v(dim,num))
           
           v(1:dim,1:num) = this%v(1:dim,1:num)
           
        END IF
        
9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE particles_get_v
      
      
      SUBROUTINE particles_get_rho(this,rho,num,stat_info)

        TYPE(Particles),INTENT(IN)              :: this
        REAL(MK),DIMENSION(:),POINTER           :: rho
        INTEGER, INTENT(IN)                     :: num
        INTEGER,INTENT(OUT)                     :: stat_info
        
        stat_info = 0
        
        IF( num > this%num_part_all) THEN
           PRINT *, "particles_get_rho :  ",&
                "Required particles is more than existed!"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF(ASSOCIATED(rho)) THEN 
           DEALLOCATE(rho)
        END IF

        IF ( num > 0 ) THEN
           
           ALLOCATE(rho(num))
           
           rho(1:num) = this%rho(1:num)
           
        END IF
        
9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE particles_get_rho

      
      SUBROUTINE particles_get_rho_norm(this,rho_norm,num,stat_info)

        TYPE(Particles),INTENT(IN)              :: this
        REAL(MK),DIMENSION(:),POINTER           :: rho_norm
        INTEGER, INTENT(IN)                     :: num
        INTEGER,INTENT(OUT)                     :: stat_info
        
        stat_info = 0
        
        IF( num > this%num_part_all) THEN
           PRINT *, "particles_get_rho_norm :  ",&
                "Required particles is more than existed!"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF(ASSOCIATED(rho_norm)) THEN 
           DEALLOCATE(rho_norm)
        END IF

        IF ( num > 0 ) THEN
           
           ALLOCATE(rho_norm(num))
           
           rho_norm(1:num) = this%rho_norm(1:num)
           
        END IF
        
9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE particles_get_rho_norm

      
      REAL(MK) FUNCTION particles_get_rho_min(this,stat_info)
        
        TYPE(Particles),INTENT(IN)              :: this
        INTEGER,INTENT(OUT)                     :: stat_info
        
        stat_info = 0
        
        particles_get_rho_min = this%rho_min
        
        RETURN          
        
      END FUNCTION particles_get_rho_min
      

      REAL(MK) FUNCTION particles_get_rho_max(this,stat_info)
        
        TYPE(Particles),INTENT(IN)              :: this
        INTEGER,INTENT(OUT)                     :: stat_info
        
        stat_info = 0
        
        particles_get_rho_max = this%rho_max
        
        RETURN          
        
      END FUNCTION particles_get_rho_max
      

      SUBROUTINE particles_get_m(this,m,num,stat_info)

        TYPE(Particles),INTENT(IN)              :: this
        REAL(MK),DIMENSION(:),POINTER           :: m
        INTEGER, INTENT(IN)                     :: num
        INTEGER,INTENT(OUT)                     :: stat_info
        
        
        stat_info = 0
        
        IF( num > this%num_part_all) THEN
           PRINT *, "particles_get_m :  ",&
                "Required particles is more than existed !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF(ASSOCIATED(m)) THEN 
           DEALLOCATE(m)
        END IF
        
        IF ( num > 0 ) THEN
           
           ALLOCATE( m(num))
           
           m(1:num) = this%m(1:num)
           
        END IF
           
9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE particles_get_m
      
      
      SUBROUTINE particles_get_p(this,p,num,stat_info)

        TYPE(Particles),INTENT(IN)              :: this
        REAL(MK),DIMENSION(:),POINTER           :: p
        INTEGER, INTENT(IN)                     :: num
        INTEGER,INTENT(OUT)                     :: stat_info
        
        stat_info = 0
        
        IF( num > this%num_part_all) THEN
           PRINT *, "particles_get_p :  ",&
                "Required particles is more than existed !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF(ASSOCIATED(p)) THEN 
           DEALLOCATE(p)
        END IF
        
        IF ( num > 0 ) THEn
           
           ALLOCATE( p(num) )
           
           p(1:num) = this%p(1:num)
           
        END IF
           
9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE particles_get_p
      
      
      SUBROUTINE particles_get_id(this,id,num,stat_info)
        
        TYPE(Particles),INTENT(IN)              :: this
        INTEGER,DIMENSION(:,:),POINTER          :: id
        INTEGER, INTENT(IN)                     :: num
        INTEGER,INTENT(OUT)                     :: stat_info
        
        !-----------------------
        !Local variables
        !------------------------
        
        INTEGER                         :: dim
        
        stat_info = 0
        
        IF( num > this%num_part_all) THEN
           PRINT *, "particles_get_id :  ",&
                "Required particles is more than existed !"
           stat_info = - 1
           GOTO 9999
        END IF
        
        IF(ASSOCIATED(id)) THEN 
           DEALLOCATE(id)
        END IF
        
        IF ( num > 0 ) THEN
           
           dim =this%num_id
           
           ALLOCATE(id(dim,num))
           
           id(1:dim,1:num) = this%id(1:dim,1:num)
           
        END IF
          
9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE particles_get_id


      SUBROUTINE particles_get_pid(this,pid,num,stat_info)

        TYPE(Particles),INTENT(IN)              :: this
        INTEGER,DIMENSION(:),POINTER            :: pid
        INTEGER, INTENT(IN)                     :: num
        INTEGER,INTENT(OUT)                     :: stat_info
        
        stat_info = 0
        
        IF( num > this%num_part_all) THEN
           PRINT *, "particles_get_pid :  ",&
                "Required particles is more than existed !"
           stat_info = -1 
           GOTO 9999
        END IF
        
        IF(ASSOCIATED(pid)) THEN 
           DEALLOCATE(pid)
        END IF
        
        ALLOCATE(pid(num))
        
        IF ( num > 0 ) THEN
           
           pid(1:num) = this%id(this%pid_idx,1:num)
           
        END IF
        
9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE particles_get_pid
      
       
      SUBROUTINE particles_get_sid(this,sid,num,stat_info)

        TYPE(Particles),INTENT(IN)              :: this
        INTEGER,DIMENSION(:),POINTER            :: sid
        INTEGER, INTENT(IN)                     :: num
        INTEGER,INTENT(OUT)                     :: stat_info
        
        stat_info = 0
        
        IF( num > this%num_part_all) THEN
           PRINT *, "particles_get_sid :  ",&
                "Required particles is more than existed !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF(ASSOCIATED(sid)) THEN 
           DEALLOCATE(sid)
        END IF
        
        IF ( num > 0 ) THEN
           
           ALLOCATE(sid(num))
           
           sid(1:num) = this%id(this%sid_idx,1:num)
           
        END IF

9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE particles_get_sid
      
      
      INTEGER FUNCTION particles_get_num_id(this,stat_info)
        
        TYPE(Particles),INTENT(IN)              :: this
        INTEGER,INTENT(OUT)                     :: stat_info

        stat_info = 0
        particles_get_num_id = this%num_id

        RETURN

      END FUNCTION particles_get_num_id
      
      
      SUBROUTINE particles_get_f(this,f,num,stat_info)

        TYPE(Particles),INTENT(IN)              :: this
        REAL(MK),DIMENSION(:,:),POINTER         :: f
        INTEGER, INTENT(IN)                     :: num
        INTEGER,INTENT(OUT)                     :: stat_info
        
        !-----------------------
        !Local variables
        !------------------------
        
        INTEGER                         :: dim

        stat_info = 0
        
        IF( num > this%num_part_all) THEN
           PRINT *, "particles_get_f :  ",&
                "Required particles is more than existed !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF(ASSOCIATED(f)) THEN 
           DEALLOCATE(f)
        END IF
        
        IF( num > 0 ) THEN
           
           dim = SIZE(this%f,1)
           
           ALLOCATE(f(dim,num))
           
           f(1:dim,1:num) = this%f(1:dim,1:num)
           
        END IF
        
9999    CONTINUE
        
        RETURN
        
      END SUBROUTINE particles_get_f
      
      
      REAL(MK) FUNCTION particles_get_fa_min(this,stat_info)
        
        TYPE(Particles),INTENT(IN)              :: this
        INTEGER,INTENT(OUT)                     :: stat_info
        
        stat_info = 0
        
        particles_get_fa_min = this%fa_min
        
        RETURN          
        
      END FUNCTION particles_get_fa_min
      

      REAL(MK) FUNCTION particles_get_fa_max(this,stat_info)
        
        TYPE(Particles),INTENT(IN)              :: this
        INTEGER,INTENT(OUT)                     :: stat_info
        
        stat_info = 0
        
        particles_get_fa_max = this%fa_max
        
        RETURN          
        
      END FUNCTION particles_get_fa_max
      
      
      REAL(MK) FUNCTION particles_get_dt_f(this,stat_info)
        
        TYPE(Particles),INTENT(IN)              :: this
        INTEGER,INTENT(OUT)                     :: stat_info
        
        stat_info = 0
        
        particles_get_dt_f = this%dt_f
        
        RETURN          
        
      END FUNCTION particles_get_dt_f
      

      SUBROUTINE particles_get_vgt(this,vgt,num,stat_info)
        
        TYPE(Particles),INTENT(IN)              :: this
        REAL(MK),DIMENSION(:,:),POINTER         :: vgt
        INTEGER, INTENT(IN)                     :: num
        INTEGER,INTENT(OUT)                     :: stat_info
        
        !-----------------------
        ! Local variables
        !------------------------
        
        INTEGER                       :: dim
        
        stat_info = 0
        
        IF( num > this%num_part_all) THEN
           PRINT *, "particles_get_vgt :  ",&
                "Required particles is more than existed !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF(ASSOCIATED(vgt)) THEN 
           DEALLOCATE(vgt)
        END IF
        
        IF ( num > 0 ) THEN
           
           dim = SIZE(this%vgt,1)
           
           ALLOCATE( vgt(dim,num) )
           
           vgt(1:dim,1:num) = this%vgt(1:dim,1:num)
           
        END IF
        
9999    CONTINUE
   
        RETURN          
        
      END SUBROUTINE particles_get_vgt
      

      SUBROUTINE particles_get_evgt(this,evgt,num,stat_info)
        
        TYPE(Particles),INTENT(IN)              :: this
        REAL(MK),DIMENSION(:,:),POINTER         :: evgt
        INTEGER, INTENT(IN)                     :: num
        INTEGER,INTENT(OUT)                     :: stat_info
        
        !-----------------------
        ! Local variables
        !------------------------
        
        INTEGER                       :: dim
        
        stat_info = 0
        
        IF( num > this%num_part_all) THEN
           PRINT *, "particles_get_evgt :  ",&
                "Required particles is more than existed !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF(ASSOCIATED(evgt)) THEN 
           DEALLOCATE(evgt)
        END IF
        
        IF ( num > 0 ) THEN
           
           dim = SIZE(this%evgt,1)
           
           ALLOCATE( evgt(dim,num) )
           
           evgt(1:dim,1:num) = this%evgt(1:dim,1:num)
           
        END IF
        
9999    CONTINUE
   
        RETURN          
        
      END SUBROUTINE particles_get_evgt

      
      SUBROUTINE particles_get_evec(this,evec,num,stat_info)

        TYPE(Particles),INTENT(IN)              :: this
        REAL(MK),DIMENSION(:,:),POINTER         :: evec
        INTEGER, INTENT(IN)                     :: num
        INTEGER,INTENT(OUT)                     :: stat_info
        
        !-----------------------
        ! Local variables
        !------------------------
        
        INTEGER                       :: dim
        
        stat_info = 0
        
        IF( num > this%num_part_all) THEN
           PRINT *, "particles_get_evec :  ",&
                "Required particles is more than existed !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF(ASSOCIATED(evec)) THEN 
           DEALLOCATE(evec)
        END IF
        
        IF ( num > 0 ) THEN
           
           dim = SIZE(this%evec,1)
           
           ALLOCATE(evec(dim,num) )
           
           evec(1:dim,1:num) = this%evec(1:dim,1:num)
           
        END IF
        
9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE particles_get_evec
      
      SUBROUTINE particles_get_eval(this,eval,num,stat_info)

        TYPE(Particles),INTENT(IN)              :: this
        REAL(MK),DIMENSION(:,:),POINTER         :: eval
        INTEGER, INTENT(IN)                     :: num
        INTEGER,INTENT(OUT)                     :: stat_info
        
        !-----------------------
        ! Local variables
        !------------------------
        
        INTEGER                       :: dim
        
        stat_info = 0
        
        IF( num > this%num_part_all) THEN
           PRINT *, "particles_get_eval :  ",&
                "Required particles is more than existed !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF(ASSOCIATED(eval)) THEN 
           DEALLOCATE(eval)
        END IF
        
        IF ( num > 0 ) THEN
           
           dim = SIZE(this%eval,1)
           
           ALLOCATE(eval(dim,num) )
           
           eval(1:dim,1:num) = this%eval(1:dim,1:num)
           
        END IF
        
9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE particles_get_eval


      SUBROUTINE particles_get_aevec(this,aevec,num,stat_info)

        TYPE(Particles),INTENT(IN)              :: this
        REAL(MK),DIMENSION(:,:),POINTER         :: aevec
        INTEGER, INTENT(IN)                     :: num
        INTEGER,INTENT(OUT)                     :: stat_info
        
        !-----------------------
        ! Local variables
        !------------------------
        
        INTEGER                       :: dim
        
        stat_info = 0
        
        IF( num > this%num_part_all) THEN
           PRINT *, "particles_get_aevec :  ",&
                "Required particles is more than existed !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF(ASSOCIATED(aevec)) THEN 
           DEALLOCATE(aevec)
        END IF
        
        IF ( num > 0 ) THEN
           
           dim = SIZE(this%aevec,1)
           
           ALLOCATE(aevec(dim,num) )
           
           aevec(1:dim,1:num) = this%aevec(1:dim,1:num)
           
        END IF
        
9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE particles_get_aevec
      

      SUBROUTINE particles_get_pt(this,pt,num,stat_info)

        TYPE(Particles),INTENT(IN)              :: this
        REAL(MK),DIMENSION(:,:,:),POINTER       :: pt
        INTEGER, INTENT(IN)                     :: num
        INTEGER,INTENT(OUT)                     :: stat_info
        
        !-----------------------
        ! Local variables
        !------------------------
        
        INTEGER                                 :: dim
        
        stat_info = 0
        
        IF( num > this%num_part_all) THEN
           PRINT *, "particles_get_pt :  ",&
                "Required particles is more than existed !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF(ASSOCIATED(pt)) THEN 
           DEALLOCATE(pt)
        END IF
        
        IF ( num > 0 ) THEN
           
           dim = SIZE(this%pt,1)
           
           ALLOCATE(pt(dim,dim,num))
           
           pt(1:dim,1:dim,1:num) = this%pt(1:dim,1:dim,1:num)
           
        END IF
        
9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE particles_get_pt


      SUBROUTINE particles_get_ct(this,ct,num,stat_info)

        TYPE(Particles),INTENT(IN)              :: this
        REAL(MK),DIMENSION(:,:),POINTER         :: ct
        INTEGER, INTENT(IN)                     :: num
        INTEGER,INTENT(OUT)                     :: stat_info
        
        !-----------------------
        ! Local variables
        !------------------------
        
        INTEGER                       :: dim
        
        stat_info = 0
        
        IF( num > this%num_part_all) THEN
           PRINT *, "particles_get_ct :  ",&
                "Required particles is more than existed !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF(ASSOCIATED(ct)) THEN 
           DEALLOCATE(ct)
        END IF
        
        IF ( num > 0 ) THEN
           
           dim = SIZE(this%ct,1)
           
           ALLOCATE(ct(dim,num) )
           
           ct(1:dim,1:num) = this%ct(1:dim,1:num)
           
        END IF
        
9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE particles_get_ct
      
      
      SUBROUTINE particles_get_act(this,act,num,stat_info)
        
        TYPE(Particles),INTENT(IN)              :: this
        REAL(MK),DIMENSION(:,:),POINTER         :: act
        INTEGER, INTENT(IN)                     :: num
        INTEGER,INTENT(OUT)                     :: stat_info
        
        !-----------------------
        ! Local variables
        !------------------------
        
        INTEGER                       :: dim
        
        stat_info = 0
        
        IF( num > this%num_part_all) THEN
           PRINT *, "particles_get_act :  ",&
                "Required particles is more than existed !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF(ASSOCIATED(act)) THEN 
           DEALLOCATE(act)
        END IF
        
        IF ( num > 0 )  THEN

           dim = SIZE(this%act,1)
           
           ALLOCATE( act(dim,num) )
           
           act(1:dim,1:num) = this%act(1:dim,1:num)
           
        END IF

9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE particles_get_act
      
      
      SUBROUTINE particles_get_u(this,u,num,stat_info)
        
        TYPE(Particles),INTENT(IN)              :: this
        REAL(MK),DIMENSION(:),POINTER           :: u
        INTEGER, INTENT(IN)                     :: num
        INTEGER,INTENT(OUT)                     :: stat_info
        
        
        stat_info = 0
        
        IF( num > this%num_part_all) THEN
           PRINT *, "particles_get_u :  ",&
                "Required particles is more than existed !"
           stat_info = -1
           GOTO 9999
        END IF
        
        IF(ASSOCIATED(u)) THEN 
           DEALLOCATE(u)
        END IF
        
        IF ( num > 0 ) THEN
           
           ALLOCATE( u(num) )
           
           u(1:num) = this%u(1:num)
           
        END IF
        
9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE particles_get_u
      
      
      SUBROUTINE particles_get_au(this,au,num,stat_info)
        
        TYPE(Particles),INTENT(IN)              :: this
        REAL(MK),DIMENSION(:),POINTER           :: au
        INTEGER, INTENT(IN)                     :: num
        INTEGER,INTENT(OUT)                     :: stat_info
        
        
        stat_info = 0
        
        IF( num > this%num_part_all) THEN
           PRINT *, "particles_get_au :  ",&
                "Required particles is more than existed !"
        END IF
        
        IF(ASSOCIATED(au)) THEN 
           DEALLOCATE(au)
        END IF
        
        IF ( num > 0 ) THEN
           
           ALLOCATE( au(num) )
           
           au(1:num) = this%au(1:num)
           
        END IF
        
9999    CONTINUE
        
        RETURN          
        
      END SUBROUTINE particles_get_au

      
      INTEGER FUNCTION particles_get_num_part_real(this,stat_info)
        
        TYPE(Particles),INTENT(IN)           :: this
        INTEGER,INTENT(OUT)                  :: stat_info
        
        stat_info = 0
        
        particles_get_num_part_real = &
             this%num_part_real
        
        RETURN          
        
      END FUNCTIOn particles_get_num_part_real
      
      
      INTEGER FUNCTION particles_get_num_part_all(this,stat_info)
        
        TYPE(Particles),INTENT(IN)           :: this
        INTEGER,INTENT(OUT)                  :: stat_info
        
        
        stat_info = 0
        
        particles_get_num_part_all = &
             this%num_part_all
                
        RETURN          

      END FUNCTIOn particles_get_num_part_all


      INTEGER FUNCTION particles_get_num_part_ghost(this,stat_info)
        
        TYPE(Particles),INTENT(IN)           :: this
        INTEGER,INTENT(OUT)                  :: stat_info
        
        
        stat_info = 0
        
        particles_get_num_part_ghost = &
             this%num_part_ghost
        
        RETURN          
        
      END FUNCTION particles_get_num_part_ghost
      
      
      INTEGER FUNCTION particles_get_num_part_colloid(this,stat_info)
        
        TYPE(Particles),INTENT(IN)           :: this
        INTEGER,INTENT(OUT)                  :: stat_info
        
        
        stat_info = 0
        
        particles_get_num_part_colloid = &
             this%num_part_colloid
                
        RETURN          

      END FUNCTIOn particles_get_num_part_colloid
      
      
      SUBROUTINE particles_get_part_colloid_list(this,part_list,stat_info)
        
        TYPE(Particles),INTENT(IN)      :: this
        INTEGER, DIMENSION(:,:),POINTER :: part_list
        INTEGER,INTENT(OUT)             :: stat_info
        
        !----------------------
        ! Local variables.
        !----------------------

        INTEGER                         :: dim1, dim2
        
        stat_info = 0
        
        IF(ASSOCIATED(part_list)) THEN
           DEALLOCATE(part_list)
        END IF
        
        !------------------------------------------
        ! Copy particle list of colloid boundary 
        ! particles only if there are some.
        !------------------------------------------
        
        IF (this%num_part_colloid > 0 .AND. &
             ASSOCIATED(this%part_colloid_list)) THEN
           
           dim1 = SIZE(this%part_colloid_list,1)
           dim2 = SIZE(this%part_colloid_list,2)
           
           ALLOCATE(part_list(dim1,dim2))
           part_list(1:dim1,1:dim2) = &
                this%part_colloid_list(1:dim1,1:dim2)
           
        END IF
        
        RETURN
        
      END SUBROUTINE  particles_get_part_colloid_list
      
      
