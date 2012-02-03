#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Module       :                ppm_module_io_open
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module contains all data structures and
      !                 definitions that are PRIVATE to the IO routines.
      !                 It also included those routines and provides
      !                 INTERFACEs.
      !                
      !  Remarks      : 
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_io_open.f,v $
      !  Revision 1.1  2004/07/26 07:29:39  ivos
      !  First commit after spitting the old modules into single-interface
      !  units.
      !
      !-------------------------------------------------------------------------
      !  Perallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschengraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------

      MODULE ppm_module_io_open

         !----------------------------------------------------------------------
         !  Define interface to ppm_io_open
         !----------------------------------------------------------------------
         INTERFACE ppm_io_open
             MODULE PROCEDURE ppm_io_open
         END INTERFACE

         !----------------------------------------------------------------------
         !  Include the source
         !----------------------------------------------------------------------
         CONTAINS

#include "ppm_io_open.F90"

      END MODULE ppm_module_io_open
