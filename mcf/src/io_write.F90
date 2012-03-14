      SUBROUTINE io_write(this,rank,step,time,&
           parts,num_part,statis,stat_info)
        !----------------------------------------------------
        !  Subroutine   :  io_write
        !----------------------------------------------------
        !
        !  Purpose      :  Write information into output files.
        !
        !  Reference    :
        !
        !  Remark       :
        !
        !  Revisions    : 
        !                 V0.1 13.10.2009, original version,
        !                 merge all io_write_* together here.
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
        
        TYPE(IO), INTENT(IN)            :: this
        INTEGER, INTENT(IN)             :: rank
        INTEGER, INTENT(IN)             :: step
        REAL(MK), INTENT(IN)            :: time
        TYPE(Particles),INTENT(IN)      :: parts
        INTEGER, INTENT(IN)             :: num_part
        TYPE(Statistic), INTENT(IN)     :: statis
        INTEGER, INTENT(OUT)            :: stat_info
        
        !----------------------------------------------------
        ! Local variables.
        !----------------------------------------------------
        
        INTEGER                         :: stat_info_sub
        TYPE(Colloid),POINTER           :: colloids
        TYPE(Boundary),POINTER          :: tboundary

        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        NULLIFY(colloids)
        NULLIFY(tboundary)        

        !----------------------------------------------------
        ! Write particles' to file
        ! i.e., position, velocity, mass/number density, 
        ! mass, species ID(if more than one species).
        !----------------------------------------------------

        IF ( this%write_particles ) THEN
           
           CALL io_write_particles(this,&
                rank,step,parts,num_part,stat_info_sub)

           IF (stat_info_sub /= 0) THEN
              PRINT *, 'io_write : ',&
                   'Writing particles failed !'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
                
        !----------------------------------------------------
	! Write conformation tensor to file in case
        ! of Non-Newtonian fluids.
      	!----------------------------------------------------
        
        IF ( this%write_conformation ) THEN
           
           CALL io_write_conformation(this,&
                rank,step,parts,num_part,&
                stat_info_sub)
           
           IF (stat_info_sub /= 0) THEN
              PRINT *, 'io_write : ',&
                   'Writing conformation tensor failed !'
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        !----------------------------------------------------
        ! Done by rank=0 process.
        !----------------------------------------------------
        
        IF ( rank == 0 ) THEN
           
           !-------------------------------------------------
           ! Write colloid file.
           !-------------------------------------------------
           
           IF( this%write_colloid ) THEN
              
              CALL physics_get_colloid(this%phys,&
                   colloids,stat_info_sub)
              
#ifdef __IO_COLLOID_SEPARATE
              
              CALL io_write_colloid_separate(this,rank,&
                   step,time,colloids,stat_info_sub)
              
#else
              
              CALL io_write_colloid(this,rank,&
                   step,colloids,stat_info_sub)
              
              
#endif         
              
              IF(stat_info_sub /=0 ) THEN
                 PRINT *,"io_write : ",&
                      "Step writing colloids failed !"
                 stat_info = -1
                 GOTO 9999
              END IF
                         
           END IF
                      
           !-------------------------------------------------
           ! Write statistics file.
           !-------------------------------------------------
           

           IF ( this%write_statistic ) THEN
              
              CALL io_write_statistic(this,rank,&
                   step,time,statis,stat_info_sub)
              
              IF(stat_info_sub /=0 ) THEN
                 PRINT *,"io_write : ",&
                      "Step writing statistic failed !"
                 stat_info = -1
                 GOTO 9999
              END IF
              
           END IF

                      
           !-------------------------------------------------
           ! Write boundary file.
           !-------------------------------------------------
           
           IF( this%write_boundary) THEN
              
              CALL physics_get_boundary(this%phys,tboundary,stat_info_sub)
              CALL io_write_boundary(this,rank,&
                   step,time,tboundary,stat_info_sub)
              
              IF(stat_info_sub /=0 ) THEN
                 PRINT *,"io_write : ",&
                      "Step writing boundary failed !"
                 stat_info = -1
                 GOTO 9999
              END IF
              
           END IF
           
           !-------------------------------------------------
           ! Write restart physics file.
           !-------------------------------------------------
      
           IF ( this%write_restart_physics) THEN
              
              CALL io_write_restart_physics(this,&
                   rank,step,time,&
                   this%phys,stat_info_sub)
              
              IF(stat_info_sub /=0 ) THEN
                 PRINT *,"io_write : ",&
                      "Step writing restart physics failed !"
                 stat_info = -1
                 GOTO 9999
              END IF
              
           END IF ! wrte_restart, physics
           
        END IF ! rank == 0
        
        !----------------------------------------------------
        ! Write restart particles file.
        !----------------------------------------------------
      
        IF ( this%write_restart_particles) THEN
           
           CALL io_write_restart_particles(this,&
                rank,step,parts,num_part,stat_info_sub)
           
           IF(stat_info_sub /=0 ) THEN
              PRINT *,"io_write : ",&
                   "Step writing restart particles failed !"
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF
        
        !----------------------------------------------------
        ! Write restart conformation tensor file.
        !----------------------------------------------------
        
        IF ( this%write_restart_conformation ) THEN
           
           CALL io_write_restart_conformation(this,&
                rank,step,parts,num_part,stat_info_sub)
           
           IF(stat_info_sub /=0 ) THEN
              PRINT *,"io_write : ",&
                   "Step writing restart conformation failed !"
              stat_info = -1
              GOTO 9999
           END IF
        END IF


9999    CONTINUE
        
        
        RETURN
        
      END SUBROUTINE io_write
      
#include "io_write_particles_relax.F90"      
#include "io_write_particles.F90"
#include "io_write_conformation.F90"
#include "io_write_statistic_relax.F90"
#include "io_write_statistic.F90"
#include "io_write_boundary.F90"
#include "io_write_colloid.F90"
#include "io_write_colloid_separate.F90"
#include "io_write_restart_physics.F90"
#include "io_write_restart_particles_relax.F90"
#include "io_write_restart_particles.F90"
#include "io_write_restart_conformation.F90"

      
