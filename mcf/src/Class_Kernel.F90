      MODULE Class_Kernel
        !----------------------------------------------------
      	!  Class      : Kernel
	!----------------------------------------------------
      	!
      	!  Purpose    :
	!> \brief       Variables and corresponding operations
        !>              for kernel quantities.
      	!>	   	
        !>              The variable memebers are public 
        !
        !  Remarks    : Since there is a 'SAVE' after 'IMPLICIT NONE',
        !               all procedures using this Module sharing the
        !               same copy of this Module.
        !
      	!
      	!  References :
     	!
      	!  Revisions  : 0.1 03.03.2009, original version.
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
   
        USE mcf_header

        IMPLICIT NONE
        SAVE 

        !----------------------------------------------------
        !  num_dim     : number of dimension.
        !  kernel_type : type of kernel.
        !  cut_off     : cut off.
        !  h           : smoothing length.
        !----------------------------------------------------
        
        TYPE Kernel
           
           PRIVATE
           
           INTEGER              :: num_dim
           INTEGER              :: kernel_type
           REAL(MK)             :: h
           REAL(MK)             :: cut_off
           REAL(MK)             :: coef
           REAL(MK)             :: coef_grad
           
        END TYPE Kernel
        
        INTERFACE kernel_new
           MODULE PROCEDURE kernel_init
        END INTERFACE

        INTERFACE kernel_kernel
           MODULE PROCEDURE kernel_kernel_w
           MODULE PROCEDURE kernel_kernel_w_gradw
        END INTERFACE
        
        INTERFACE kernel_kernel_quintic_spline
           MODULE PROCEDURE kernel_kernel_quintic_spline_w
           MODULE PROCEDURE kernel_kernel_quintic_spline_w_gradw
        END INTERFACE
        
        INTERFACE kernel_kernel_Lucy
           MODULE PROCEDURE kernel_kernel_Lucy_w
           MODULE PROCEDURE kernel_kernel_Lucy_w_gradw
        END INTERFACE
        

      CONTAINS

#include "kernel_new.F90"
#include "kernel_finalize.F90"
#include "kernel_kernel.F90"
#include "kernel_get.F90"


      END MODULE Class_Kernel

