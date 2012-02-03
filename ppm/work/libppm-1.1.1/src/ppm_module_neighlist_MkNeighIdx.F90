#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

      !-------------------------------------------------------------------------
      !  Module       :           ppm_module_neighlist_MkNeighIdx
      !-------------------------------------------------------------------------
      !
      !  Purpose      : This module includes the source code for neighbor
      !                 search routines (cell lists, Verlet lists).
      !
      !  Remarks      :
      !
      !  References   :
      !
      !  Revisions    :
      !-------------------------------------------------------------------------
      !  $Log: ppm_module_neighlist_MkNeighIdx.f,v $
      !  Revision 1.1  2004/07/26 07:30:02  ivos
      !  First commit after spitting the old modules into single-interface
      !  units.
      !
      !-------------------------------------------------------------------------
      !  Parallel Particle Mesh Library (PPM)
      !  Institute of Computational Science
      !  ETH Zentrum, Hirschebngraben 84
      !  CH-8092 Zurich, Switzerland
      !-------------------------------------------------------------------------
     
      MODULE ppm_module_neighlist_MkNeighIdx

         !----------------------------------------------------------------------
         !  Define interface to ppm_neighlist_MkNeighIdx
         !----------------------------------------------------------------------
         INTERFACE ppm_neighlist_MkNeighIdx
            MODULE PROCEDURE ppm_neighlist_MkNeighIdx
         END INTERFACE

         !----------------------------------------------------------------------
         !  include the source 
         !----------------------------------------------------------------------
         CONTAINS

#include "ppm_neighlist_MkNeighIdx.F90"

      END MODULE ppm_module_neighlist_MkNeighIdx
