#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Subroutine   :                    ppm_fmm_expansion
      !-------------------------------------------------------------------------
      !
      !  Purpose      : Compute the expansions of the leaf boxes of the created
      !                 tree structure.
      !                 Calls ppm_fmm_traverse, which traverses the tree and 
      !                 shifts the expansions up the tree.
      !                 
      !
      !  Input        : xpunord(:,:)      (F) particle positions
      !                 wpunord(:,:)      (F) particle strengths
      !                 Np                (I) number of particles.
      !                 lda               (I) leading dimension of vector case
      !
      !  Input/output :     
      !
      !  Output       :
      !                 info         (I) return status. 0 upon success.
      !
      !  Remarks      :  
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_fmm_expansion.f,v $
      !  Revision 1.20  2006/09/04 18:34:45  pchatela
      !  Fixes and cleanups to make it compile for
      !  release-1-0
      !
      !  Revision 1.19  2006/06/29 10:28:34  pchatela
      !  Added vector strengths support
      !
      !  Revision 1.18  2006/06/20 15:17:05  hiebers
      !  BUGFIX: adjusted indices of ppm_boxid and ppm_subid
      !
      !  Revision 1.17  2006/06/16 07:52:20  hiebers
      !  Added a new list of topo IDs (topoidlist) to prevent overwriting user defined
      !  topologies
      !
      !  Revision 1.16  2006/02/03 09:39:12  ivos
      !  Changed the expansion loop to allow vectorization.
      !
      !  Revision 1.15  2005/09/19 13:03:27  polasekb
      !  code cosmetics
      !
      !  Revision 1.14  2005/09/11 11:43:59  polasekb
      !  moved mapping and second tree call to init
      !
      !  Revision 1.13  2005/08/25 13:52:28  polasekb
      !  mapping of boxpart added
      !
      !  Revision 1.12  2005/08/23 14:22:43  polasekb
      !  corrected to xpunord and wpunord
      !
      !  Revision 1.11  2005/08/23 14:11:38  polasekb
      !  fixed the formula for the Cnm
      !
      !  Revision 1.10  2005/08/11 15:13:15  polasekb
      !  now using maxboxcost from the data file
      !
      !  Revision 1.9  2005/08/08 13:33:41  polasekb
      !  init some more variables
      !
      !  Revision 1.8  2005/08/04 16:01:01  polasekb
      !  moved some data allocation to init
      !
      !  Revision 1.7  2005/07/29 12:35:27  polasekb
      !  changed diagonal to radius
      !
      !  Revision 1.6  2005/07/27 15:01:27  polasekb
      !  changed tree variables
      !
      !  Revision 1.5  2005/07/25 14:39:54  polasekb
      !  adapted to new constants saved in the
      !  ppm_module_data_fmm file
      !  adapted function call to ppm_fmm_traverse
      !
      !  Revision 1.4  2005/07/21 13:20:25  polasekb
      !  nullify lpdx, lhbx pointers
      !
      !  Revision 1.3  2005/06/02 14:18:39  polasekb
      !  removed variable totalmass
      !  changed some comments
      !
      !  Revision 1.2  2005/05/27 12:45:29  polasekb
      !  removed debug output
      !
      !  Revision 1.1  2005/05/27 07:56:40  polasekb
      !  initial implementation
      !
      !  Revision 0  2004/11/16 15:59:14 polasekb
      !  start
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

#if   (__KIND == __SINGLE_PRECISION && __DIM == __SFIELD)
      SUBROUTINE ppm_fmm_expansion_s_sf(xpunord,wpunord,Np,info)
#elif (__KIND == __DOUBLE_PRECISION && __DIM == __SFIELD)
      SUBROUTINE ppm_fmm_expansion_d_sf(xpunord,wpunord,Np,info)
#elif (__KIND == __SINGLE_PRECISION && __DIM == __VFIELD)
      SUBROUTINE ppm_fmm_expansion_s_vf(xpunord,wpunord,lda,Np,info)
#elif (__KIND == __DOUBLE_PRECISION && __DIM == __VFIELD)
      SUBROUTINE ppm_fmm_expansion_d_vf(xpunord,wpunord,lda,Np,info)
#endif


      !-------------------------------------------------------------------------
      !  Modules 
      !-------------------------------------------------------------------------
      USE ppm_module_data
      USE ppm_module_data_fmm
      USE ppm_module_alloc
      USE ppm_module_error
      USE ppm_module_fmm_traverse
      USE ppm_module_substart
      USE ppm_module_substop 
      USE ppm_module_util_cart2sph
      USE ppm_module_write
      IMPLICIT NONE
      
      !-------------------------------------------------------------------------
      !  Includes
      !-------------------------------------------------------------------------
#include "ppm_define.h"

#ifdef HAVE_MPI
      INCLUDE 'mpif.h'
#else
#include "fakempi.h"
#endif

      !-------------------------------------------------------------------------
      !  Precision
      !-------------------------------------------------------------------------      
#if    __KIND == __SINGLE_PRECISION
      INTEGER, PARAMETER :: MK = ppm_kind_single

#else
      INTEGER, PARAMETER :: MK = ppm_kind_double

#endif
      
      !-------------------------------------------------------------------------
      !  Arguments     
      !-------------------------------------------------------------------------
      REAL(MK),DIMENSION(:,:),POINTER       :: xpunord 
      INTEGER                ,INTENT(INOUT) :: Np
      INTEGER                ,INTENT(  OUT) :: info

#if   __DIM == __SFIELD
      REAL(MK),DIMENSION(:)  ,POINTER       :: wpunord

#else
      REAL(MK),DIMENSION(:,:),POINTER       :: wpunord
      INTEGER                               :: lda

#endif

      !-------------------------------------------------------------------------
      !  Local variables 
      !-------------------------------------------------------------------------
      ! auxiliary variables
      LOGICAL ,DIMENSION(3)                 :: fixed
      INTEGER                               :: iopt,i,j,k,m,n,l
      INTEGER                               :: first,last,box,nrpart
      INTEGER ,DIMENSION(2)                 :: ldu2 
      REAL(MK)                              :: sine,cosine,val,prod 
      REAL(MK)                              :: angle,reci 
      REAL(MK)                              :: t0,x0,y0,z0
      REAL(MK)                              :: dx,dy,dz,dist
      REAL(MK),DIMENSION(:,:), POINTER      :: min_box,max_box
      REAL(MK),DIMENSION(:,:), POINTER      :: min_sub,max_sub
      COMPLEX(MK),PARAMETER                 :: CI=(0.0_MK,1.0_MK)
      COMPLEX(MK)                           :: temp
      CHARACTER(LEN=ppm_char)               :: cbuf
      
      ! parallelisation
      INTEGER                               :: isub,topoid
      INTEGER                               :: nsublist,in_topoid
      REAL(MK),DIMENSION(:)  , POINTER      :: boxcost 
      
      ! traversing
      INTEGER                               :: root
      
      ! fmm
      REAL(MK),DIMENSION(:)  , POINTER      :: rho,theta,phi,radius 
      REAL(MK),DIMENSION(:)  , POINTER      :: fac,fracfac
      REAL(MK),DIMENSION(:,:), POINTER      :: Anm,Pnm,sqrtfac,xp 
      REAL(MK),DIMENSION(:,:), POINTER      :: centerofbox
      COMPLEX(MK),DIMENSION(:,:)  ,POINTER  :: Ynm
      
#if   __DIM == __SFIELD
      COMPLEX(MK),DIMENSION(:,:)  ,POINTER  :: Cnm  
      COMPLEX(MK),DIMENSION(:,:,:),POINTER  :: expansion    

#else
      COMPLEX(MK),DIMENSION(:,:,:)  ,POINTER:: Cnm    
      COMPLEX(MK),DIMENSION(:,:,:,:),POINTER:: expansion  
#endif
      

      !-------------------------------------------------------------------------
      !  Initialize 
      !-------------------------------------------------------------------------
      CALL substart('ppm_fmm_expansion',t0,info)
      



      !-------------------------------------------------------------------------
      !  Check arguments
      !-------------------------------------------------------------------------
      IF (ppm_debug.GT.0) THEN  
            IF (Np.LT. 0) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_fmm_expansion',   &
      &               'number of particles must be >= 0 !',__LINE__,info)
               GOTO 9999
            ENDIF
            IF (.NOT. ppm_fmm_initialized) THEN
               info = ppm_error_error
               CALL ppm_error(ppm_err_argument,'ppm_fmm_expansion',   &
      &               'Please call ppm_fmm_init first',__LINE__,info)
               GOTO 9999
            ENDIF

      ENDIF
      
      !-------------------------------------------------------------------------
      ! Checking precision and pointing tree data to correct variables
      !-------------------------------------------------------------------------
#if   __KIND == __SINGLE_PRECISION
      min_box      => min_box_s
      max_box      => max_box_s
      boxcost      => boxcost_s
      centerofbox  => centerofbox_s
      radius       => radius_s
#if   __DIM == __SFIELD
      expansion    => expansion_s_sf
      Cnm          => Cnm_s_sf
#else
      expansion    => expansion_s_vf 
      Cnm          => Cnm_s_vf           
#endif
      Anm          => Anm_s
      sqrtfac      => sqrtfac_s
      fracfac      => fracfac_s
      fac          => fac_s
      Ynm          => Ynm_s
      Pnm          => Pnm_s
      rho          => rho_s
      theta        => theta_s
      phi          => phi_s

#else
      min_box      => min_box_d
      max_box      => max_box_d
      boxcost      => boxcost_d
      centerofbox  => centerofbox_d
      radius       => radius_d
#if   __DIM == __SFIELD
      expansion    => expansion_d_sf
      Cnm          => Cnm_d_sf
#else
      expansion    => expansion_d_vf
      Cnm          => Cnm_d_vf            
#endif
      Anm          => Anm_d
      sqrtfac      => sqrtfac_d
      fracfac      => fracfac_d
      fac          => fac_d
      Ynm          => Ynm_d
      Pnm          => Pnm_d
      rho          => rho_d
      theta        => theta_d
      phi          => phi_d

#endif

      !-------------------------------------------------------------------------
      ! Allocating temporary array for sorted particle positions 
      !     according to tree
      !-------------------------------------------------------------------------
      iopt = ppm_param_alloc_fit
      ldu2(1) = ppm_dim
      ldu2(2) = Np
      CALL ppm_alloc(xp,ldu2,iopt,info)
      IF (info .NE. 0) THEN
         info = ppm_error_fatal
         CALL ppm_error(ppm_err_alloc,'ppm_fmm_expansion', &
      &       'error allocating xp',__LINE__,info)
      GOTO 9999
      ENDIF  
      
      !-------------------------------------------------------------------------
      ! Computing the expansions of the leafs 
      !-------------------------------------------------------------------------
      topoid = topoidlist(nlevel)
      in_topoid = ppm_internal_topoid(topoid)
      
      DO i=1,ppm_nsublist(in_topoid)
         
	 !----------------------------------------------------------------------
         ! initialize arrays
         !----------------------------------------------------------------------
         DO j=0,order
            DO k= -order,order
#if               __DIM == __SFIELD
               Cnm(j,k) = 0.0_MK
#else
               DO l=1,lda
                  Cnm(l,j,k) = 0.0_MK
               ENDDO
#endif
               Ynm(j,k) = 0.0_MK
               Pnm(j,k) = 0.0_MK
            ENDDO
         ENDDO

	 !----------------------------------------------------------------------
         ! Store the subid and boxid and the indices of the first and the last 
         ! particle in the box
         !----------------------------------------------------------------------
         isub = ppm_isublist(i,in_topoid)
         box = ppm_boxid(isub,nlevel)
         first = lhbx(1,box)
         last  = lhbx(2,box)
	 
	 !----------------------------------------------------------------------
         ! Sort xp (particle order from tree)
         !----------------------------------------------------------------------
         IF(ppm_dim.EQ.2)THEN
            DO j=first,last
               xp(1,j) = xpunord(1,lpdx(j))
               xp(2,j) = xpunord(2,lpdx(j))
            ENDDO
         ENDIF
         IF(ppm_dim.EQ.3)THEN
            DO j=first,last
               xp(1,j) = xpunord(1,lpdx(j))
               xp(2,j) = xpunord(2,lpdx(j))
               xp(3,j) = xpunord(3,lpdx(j))
            ENDDO
         ENDIF
         
	 nrpart = last-first + 1

	 !----------------------------------------------------------------------
         ! Compute the spherical coord. of the particles
         !----------------------------------------------------------------------
         IF(nrpart.GT.0)THEN
            CALL ppm_util_cart2sph(xp(1,first:last), &
                & xp(2,first:last),xp(3,first:last),nrpart, &
                & centerofbox(1,box),centerofbox(2,box),centerofbox(3,box), &
                & rho(first:last),theta(first:last),phi(first:last),info)
                IF (info .NE. 0) THEN
                   CALL ppm_error(ppm_err_sub_failed,'ppm_fmm_expansion', &
                   & 'Failed calling util_cart2sph',__LINE__,info)
                ENDIF
	 ENDIF
         
	 !----------------------------------------------------------------------
         ! Loop over the particles in the leaves
         !----------------------------------------------------------------------
	 DO j=first,last 
            
	    cosine = COS(theta(j))
            sine   = SIN(theta(j))
            
	    !-------------------------------------------------------------------
            !  Recurrence for Pnm
            !-------------------------------------------------------------------
	    val  = -sine
            prod = 1.0_MK
            
	    DO m=0,order
               Pnm(m,m) = fracfac(m)*prod
               prod     = prod * val
            ENDDO
            
	    DO m=0,order-1
               Pnm(m+1,m) = cosine*REAL(2*m + 1,MK)*Pnm(m,m)
            ENDDO
            
	    DO n=2,order
               val = cosine*REAL(2*n-1,MK)
               DO m=0,n-2
                  Pnm(n,m)=(val*Pnm(n-1,m) - & 
                  &        DBLE(n+m-1)*Pnm(n-2,m))/REAL(n-m,MK)
               ENDDO
            ENDDO
	    
	    !-------------------------------------------------------------------
            !  Compute Ynm(n,m) and Ynm(n,-m)
            !-------------------------------------------------------------------
	    DO n=0,order
               m = 0
               angle    = REAL(m,MK)*phi(j)
               Ynm(n,m) = sqrtfac(n,m)*Pnm(n,m)* &
               &          CMPLX(COS(angle),SIN(angle))
               
	       DO m=1,n
                  angle     = REAL(m,MK)*phi(j)
                  Ynm(n,m)  = sqrtfac(n,m)* &
                  &           Pnm(n,m)*CMPLX(COS(angle),SIN(angle))
                  Ynm(n,-m) = CONJG(Ynm(n,m))
               ENDDO
            ENDDO
	    
	    !-------------------------------------------------------------------
            !  Computing Cnm(n,m) - the expansion coefficients
            !-------------------------------------------------------------------
	    prod = 1.0_MK
            val  = rho(j)
            
	    DO n=0,order
               DO m=0,n
#if               __DIM == __SFIELD
                  Cnm(n,m) = Cnm(n,m) + (wpunord(lpdx(j))*prod*Ynm(n,m)* &
                  &          Anm(n,m))/((-1)**n)
#else
                  DO l =1,lda
                     Cnm(l,n,m) = Cnm(l,n,m) + (wpunord(l,lpdx(j))*prod*Ynm(n,m)* &
                  &          Anm(n,m))/((-1)**n)
                  ENDDO
#endif
               ENDDO
	       
               prod = prod * val
            ENDDO
         ENDDO !particles in one sub

	 !----------------------------------------------------------------------
         ! Computing Cnm(n,-m)
         !----------------------------------------------------------------------
	 DO n=0,order
            DO m=1,n
#if            __DIM == __SFIELD
               Cnm(n,-m) = CONJG(Cnm(n,m))

#else
               DO l =1,lda
                  Cnm(l,n,-m) = CONJG(Cnm(l,n,m))
               ENDDO
#endif 
            ENDDO
         ENDDO
         
	 DO m=1,order
            temp = CI**(-m)
            DO n=m,order
#if            __DIM == __SFIELD
               Cnm(n,m) = Cnm(n,m)*temp
               Cnm(n,-m)= Cnm(n,-m)*temp
#else
               DO l=1,lda
                  Cnm(l,n,m) = Cnm(l,n,m)*temp
                  Cnm(l,n,-m)= Cnm(l,n,-m)*temp
               ENDDO
#endif
            
	    ENDDO
         ENDDO
	 
         !----------------------------------------------------------------------
         ! Save expansion in fmm_module_data file
         !----------------------------------------------------------------------
#if      __DIM == __SFIELD
         DO n=0,order
            DO m=-order,order
               expansion(box,n,m) = Cnm(n,m)
            ENDDO
         ENDDO
#else
         DO n=0,order
            DO m=-order,order
               DO l=1,lda
                  expansion(l,box,n,m) = Cnm(l,n,m)
               ENDDO
            ENDDO
         ENDDO
#endif
      
      ENDDO !loop over all subs
      
      !-------------------------------------------------------------------------
      !  Nullify data pointers
      !-------------------------------------------------------------------------
      NULLIFY(min_box)
      NULLIFY(max_box)
      NULLIFY(boxcost)
      NULLIFY(centerofbox)
      NULLIFY(radius)
      NULLIFY(expansion)
      NULLIFY(Anm)
      NULLIFY(sqrtfac)
      NULLIFY(fracfac)
      NULLIFY(Ynm)
      NULLIFY(Pnm)
      NULLIFY(Cnm)
      
      !-------------------------------------------------------------------------
      ! Traverse the tree and shift the expansions upwards
      !-------------------------------------------------------------------------
      IF (parent(1) .EQ. ppm_param_undefined) THEN
         root = 1
      ELSE
         DO i=2,nbox
            IF (parent(i) .EQ. ppm_param_undefined) THEN
               root = i
               EXIT
            ENDIF
         ENDDO
      ENDIF

#if  __DIM == __SFIELD
      CALL ppm_fmm_traverse(root,t0,info)
#else
      CALL ppm_fmm_traverse(root,lda,t0,info)
#endif

      IF (info.NE.0) THEN
         CALL ppm_error(ppm_err_sub_failed,'ppm_fmm_expansion', &
              &         'traversing tree failed',__LINE__,info)
      ENDIF

      IF (ppm_debug .GT. 0) THEN
        CALL ppm_write(ppm_rank,'ppm_fmm_expansion','traversed tree',info)
      ENDIF     
      
      !-------------------------------------------------------------------------
      ! deallocate local data
      !-------------------------------------------------------------------------
      CALL ppm_alloc(xp,ldu2,ppm_param_dealloc,info)
      IF (info .NE. 0) THEN
          WRITE(cbuf,'(A,I3,A)') 'for ',info,'error while dealloc'
          info = ppm_error_error
          CALL ppm_error(ppm_err_dealloc,'ppm_fmm_expansion',cbuf,__LINE__,&
     &                                                                 info)
          GOTO 9999
      ENDIF

      
      !-------------------------------------------------------------------------
      !  Return
      !-------------------------------------------------------------------------
9999  CONTINUE
      
      CALL substop('ppm_fmm_expansion',t0,info)
      
      RETURN

#if   (__KIND == __SINGLE_PRECISION && __DIM == __SFIELD)
      END SUBROUTINE ppm_fmm_expansion_s_sf
#elif (__KIND == __DOUBLE_PRECISION && __DIM == __SFIELD)
      END SUBROUTINE ppm_fmm_expansion_d_sf
#elif (__KIND == __SINGLE_PRECISION && __DIM == __VFIELD)
      END SUBROUTINE ppm_fmm_expansion_s_vf
#elif (__KIND == __DOUBLE_PRECISION && __DIM == __VFIELD)
      END SUBROUTINE ppm_fmm_expansion_d_vf
#endif
