
        !-----------------------------------------------------
        ! $Log: ppm_ode_rhsfunc_macro.h,v $
        ! Revision 1.1  2004/07/26 07:52:28  michaebe
        ! Macro to make ppm_ode_step shorter and easier to read.
        !
        !-----------------------------------------------------
        IF(PRESENT(ipackdata)) THEN
           IF(PRESENT(lpackdata)) THEN
              IF(PRESENT(rpackdata)) THEN
                 throwaway = rhsfunc(xp,up,dup,lda,Npart,&
                      & ipack=ipackdata,lpack=lpackdata,rpack=rpackdata,&
                      & info=info)
              ELSE
                 throwaway = rhsfunc(xp,up,dup,lda,Npart,&
                      & ipack=ipackdata,lpack=lpackdata,info=info)
              END IF
           ELSE
              IF(PRESENT(rpackdata)) THEN
                 throwaway = rhsfunc(xp,up,dup,lda,Npart,&
                      & ipack=ipackdata,rpack=rpackdata,info=info)
              ELSE
                 throwaway = rhsfunc(xp,up,dup,lda,Npart,&
                      & ipack=ipackdata,info=info)
              END IF
           END IF
        ELSE
           IF(PRESENT(lpackdata)) THEN
              IF(PRESENT(rpackdata)) THEN
                 throwaway = rhsfunc(xp,up,dup,lda,Npart,&
                      & lpack=lpackdata,rpack=rpackdata,info=info)
              ELSE
                 throwaway = rhsfunc(xp,up,dup,lda,Npart,&
                      & lpack=lpackdata,info=info)
              END IF
           ELSE
              IF(PRESENT(rpackdata)) THEN
                 throwaway = rhsfunc(xp,up,dup,lda,Npart,&
                      & rpack=rpackdata,info=info)
              ELSE
                 throwaway = rhsfunc(xp,up,dup,lda,Npart,&
                      & info=info)
              END IF
           END IF
        END IF
        
