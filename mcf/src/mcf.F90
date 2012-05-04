!------------------------------------------------------------
! Program : Multiscale Complex Fluids simulation (MCF)
!------------------------------------------------------------
!
!	 Purpose :
!
!>	\brief Main program of 2D/3D MCF(Multiscale Complex 
!>       Fluids simulation) using PPM (Parallel Particle Mesh
!>       library).
!>
!>	Object Oriented Programming implementation for
!>      2/3D simulation of micro-<->meso-<->macroscopic 
!>	flow phenomena,	using particle method,
!>      i.e., SPH (Smoothed Particle Hydrodynamics ) and
!>            SDPD (Smoothed Dissipative Particle Dynamics).
!>	< according to Monaghan 1992;
!>          Morris et al., J. Comput. Phys. 1997;
!>          Espanol and Revenga, Phys. Rev. E 2003;
!>          Hu and Adams, J. Comput. Phys. 2006.>
!>
!>	The PPM (Parallel Particle Mesh library) has been used
!>	for the domain decomposition and inter-communication
!>      between processes for particles.
!>	< Sbalzarini et al. J. Comput. Phys. 2006 >
!
!  Routines     :
!
!  References   :
!
!  Remarks      : The code has been designed and written 
!                 in the way it appears. In any case,
!                 if you want to modify the structure or
!                 try to improve efficiency, please
!                 stop for a moment and think if you will
!                 really achieve your goal without harming
!                 the other properties of the code.
!                 It is always very helpful and meaningful
!                 to talk with the author before you make
!                 big changes, since your new ideas may 
!                 miss the big picture of the code. And
!                 suppose you have better ideas, the author
!                 would appreaciate a lot to know from you.
!                 
!
!
!> \author	  M.Sc. Bian, Xin	 
!  contact	  xin.bian@aer.mw.tum.de
!> \version	  V0.4
!> \date	  24.08.2010
!
!  Revision     : V0.4 24.08 2010, check the workflow.
!
!                 V0.3 04.12 2009, check the workflow.
!           
!                 V0.2 12.08 2009, check the workflow.
!
!                 V0.1 03.03.2009 original version.
!
!------------------------------------------------------------
! Institute	:
! Dr. Marco Ellero's Emmy Noether Group,
! at Prof. Dr. Adams' Chair of Aerodynamics, 
! Faculty of Mechanical Engineering,
! Technische Universitaet Muenchen, Germany.
!
! Copyright must be authorized by Xin Bian and Marco Ellero.
! Abuse or illegal redistribution without authorization will
! be very disappointing to the author:<( and cause legal
! matters.
!------------------------------------------------------------


!-------------------------------------------------------------
! Remark on Fortran:
! The MCF, although uses Fortran 90's grammer,
! tries to use object oriented way of programming(OOP),
! at least mimics OOP.
! For example, the file Class_Control.F90 defines the module
! Class_Control, which itself includes a derived type
! called Control. 
! This derived type Control is the name of one class.
! All the subroutines of this Class Control, start with prefix
! "control_", to be more distinguished from other subroutines.
! In main program mcf, a object of Class Control is defined with
! name mcf_ctrl. 
! Once the memeber functions/subroutiens of object
! mcf_ctrl need to be called, e.g., control_new(), the constructor,
! mcf_ctrl must be the first parameter, in order to indicate
! control_new is a member routine of object mcf_ctrl, i.e.,
! CALL control_new(mcf_ctrl, ...).
! For the rest of Classes, same rules apply.
! For OOP using Fotran 90, please refer to 
! V.K. Decyk, C.D. Norton and B.K.Szymanski,
! Introduction to Object-Orirented Concepts using Fortran 90.
!-------------------------------------------------------------


      PROGRAM mcf
        
        !----------------------------------------------------
        ! Class/Moduels used :
        !----------------------------------------------------
        
        USE mcf_header
        USE Class_Debug        
        USE CLass_Control
        USE Class_Physics
        USE Class_Technique
        USE Class_IO
        USE Class_Random
        USE Class_Rhs
        USE Class_StateEquation
        USE Class_Kernel
        USE Class_Particles
        USE Class_Marching
        
        
        IMPLICIT NONE
        
        !----------------------------------------------------
        ! Status flag of simulation.        
        !----------------------------------------------------
        
        INTEGER                         :: stat_info
        INTEGER                         :: stat_info_sub
        
        !----------------------------------------------------
        ! Create objects of different classes.
        !----------------------------------------------------
        
        TYPE(Control)                   :: mcf_ctrl
        TYPE(Physics)                   :: mcf_phys
        TYPE(Technique)                 :: mcf_tech
        TYPE(IO)                        :: mcf_io
        TYPE(Random)                    :: mcf_random
        TYPE(Rhs)                       :: mcf_rhs
        TYPE(StateEquation)             :: mcf_stateEquation
        TYPE(Kernel)                    :: mcf_kern
        TYPE(Particles)                 :: mcf_particles
        TYPE(Marching)                  :: mcf_marching
        
        !----------------------------------------------------
        ! Command line arguments :
        !
        ! n_arg      : number of command line argumetns.
        ! ctrl_file  : name of control file
        ! buf_string : buffer of characters.
        ! igroup     : index of group of processes
        ! ngroup     : total number of groups of processes.
        !----------------------------------------------------
        
        INTEGER                         :: n_arg
        CHARACTER(LEN=MAX_CHAR)         :: ctrl_file
        CHARACTER(LEN=MAX_CHAR)         :: buf_string
        INTEGER                         :: igroup
        INTEGER                         :: ngroup
        
        !----------------------------------------------------
        ! Here define the control variables, 
        ! which are read from control config file,
        ! e.g., ctrl.mcf and
        ! determine components of the program.
        !----------------------------------------------------
 
        INTEGER                         :: debug_flag
        LOGICAL                         :: relax_run
        LOGICAL                         :: read_external
        INTEGER                         :: kernel_type
        LOGICAL                         :: symmetry        
        INTEGER                         :: rhs_density_type
        INTEGER                         :: stateEquation_type
        LOGICAL                         :: Newtonian
        LOGICAL                         :: Brownian
        INTEGER                         :: random_seed
        INTEGER                         :: rhs_force_type
        LOGICAL                         :: p_energy
        INTEGER                         :: integrate_type
        INTEGER                         :: write_restart
        
        !----------------------------------------------------
        ! flags to check if control, physical and io
        ! parameters are given reasonably.
        !----------------------------------------------------
        
        LOGICAL                         :: check_ctrl
        LOGICAL                         :: check_physics
        LOGICAL                         :: check_io
        
        !----------------------------------------------------
        ! The rank in parallel programming context;
        ! Rank is 0 in serial version.
        !----------------------------------------------------

        INTEGER                         :: rank
        
        !----------------------------------------------------
        ! Physics parameters:
        !
        ! num_dim     : number of dimension
        ! min_phys_t  : minimum boundary of total physical domain.
        ! max_phys_t  : maximum boundary of total physical domain.
        ! cut_off     : cut off for the kernel compact support.
        ! h           : smoothing length.
        ! c           : sound speed.
        ! rho         : initial density.
        ! rho_ref     : reference density.
        ! gamma       : power of the state equation.
        ! bcdef       : boundary conditions.
        ! num_colloid : number of colloids, if there is.
        !----------------------------------------------------
        
        INTEGER                               :: num_dim
        REAL(MK), DIMENSION(:), POINTER       :: min_phys_t
        REAL(MK), DIMENSION(:), POINTER       :: max_phys_t
        REAL(MK)                              :: cut_off
        REAL(MK)                              :: h
        REAL(MK)                              :: c        
        REAL(MK)                              :: rho,rho_ref,gamma
        INTEGER, DIMENSION(:), POINTER        :: bcdef
        INTEGER                               :: num_colloid
        INTEGER                               :: comm
        REAL(MK), DIMENSION(2,2)              :: A
        REAL(MK), DIMENSION(2,1)              :: B
        INTEGER, DIMENSION(2)                 :: IPIV
        
#if __DEBUG
        !----------------------------------------------------
        ! Starting time of the simulation,
        ! for debug purpose, record 
        ! starting time of mcf.
        !----------------------------------------------------
        REAL(MK)    :: time_start
#endif  
#if 0
        A(1,1)=2
        A(1,2)=3
        A(2,1)=4
        A(2,2)=9
        B(1,1)=6
        B(2,1)=15
        PRINT *, "A, B:", A(:,:),B(:,:)
        CALL DGESV(2,1,A,2,IPIV(:), B,2,stat_info_sub)
        PRINT *, "X: ", B(:,:)
        PRINT *, "IPIV: ", IPIV(:)
        STOP
#endif
        !----------------------------------------------------
        ! By default, we use all the processes given as 
        ! in the same group.
        !
        ! Note that it would be helpful to divide the
        ! processes given into different groups, since
        ! some supercomputer favors "big jobs", i.e., 
        ! easier to grant high number of CPUs. 
        ! After division, several "small jobs" can share
        ! the total number of processors granted and run 
        ! simutaneously. 
        !----------------------------------------------------
        
        igroup        = 1
        ngroup        = 1
        
        !----------------------------------------------------
        ! Initialization of variables.
        !----------------------------------------------------
        
        stat_info     = 0
        stat_info_sub = 0
        
        check_ctrl    = .FALSE.
        check_physics = .FALSE.
        check_io      = .FALSE.
        
        NULLIFY(min_phys_t)
        NULLIFY(max_phys_t)
        NULLIFY(bcdef)
        
        !----------------------------------------------------
        ! Reading in arguments from command-line, 
        ! if there is any.
        !
        ! If there is no argument,
        ! default control file's name is used.
        !
        ! If there is one argument,
        ! it must be control file's name;
        !
        ! If there are 3 arguments,
        ! first is control file's name,
        ! second is igroup,
        ! third is ngroup.
        !----------------------------------------------------

        n_arg = IARGC()
        
        IF ( n_arg > 0 ) THEN
           
           CALL GETARG(1,ctrl_file)
           
           IF ( n_arg == 3 ) THEN
              CALL GETARG(2,buf_string)
              READ(buf_string,*) igroup
              CALL GETARG(3,buf_string)
              READ(buf_string,*) ngroup
              
              IF ( igroup > ngroup ) THEN
                 PRINT *, "mcf : ", &
                      "igroup must not be bigger than nroup !"
                 stat_info =  -1
                 GOTO 9999
              END IF
              
              IF ( igroup < 1 ) THEN
                 PRINT *, "mcf : ", &
                      "igroup must be positive !"
                 stat_info =  -1
                 GOTO 9999
              END IF
              
           END IF
           
        ELSE
           
           ctrl_file = "ctrl.mcf"
           
        END IF
        
        !----------------------------------------------------
        ! Initialize default values for three basic objects
        ! by calling default constructors.
        ! Control, Physics, Technique objects.
        ! (Starts MPI, if in parallel)
        !----------------------------------------------------
        
        CALL control_new(mcf_ctrl,stat_info_sub)
        CALL physics_new(mcf_phys,mcf_ctrl,stat_info_sub)
        CALL technique_new(mcf_tech,stat_info_sub)
        CALL technique_init_parallelization(mcf_tech,&
             igroup,ngroup,stat_info_sub)

        !----------------------------------------------------
        ! Get the communicator and 
        ! the rank of current process.
        !----------------------------------------------------
        
        comm = technique_get_comm(mcf_tech,stat_info_sub)
        rank = technique_get_rank(mcf_tech,stat_info_sub)
        
        !----------------------------------------------------
        ! Default run : Initialize IO object
        !
        ! NOTE : only one of default io_new() and 
        ! non-default io_new() should be used, 
        ! not at same time!
        !----------------------------------------------------
        
        !CALL io_new(mcf_io,stat_info_sub)
        
        !----------------------------------------------------
        ! Non-default run: 
        !
        ! all parameters need to be read from 
        ! configuration files. 
        ! (Control file e.g., ctrl.mcf(default), 
        ! which contains the names of 
        ! physics and io configuration files, 
        ! by default they are named as 
        ! physics_config.mcf and io_config.mcf.)
        !
        ! Initialize IO object, 
        ! which will read parameters into control and 
        ! physics objects, according to input config files.
        !
        ! NOTE : only one of default io_new() and 
        ! non-default io_new() should be used, 
        ! not same time
        !----------------------------------------------------
        
        CALL io_new(mcf_io,mcf_ctrl,ctrl_file,&
             mcf_phys,stat_info_sub)
        
        IF(stat_info_sub /= 0 ) THEN
           PRINT *, "mcf : io_new() has problem "           
           stat_info = -1
           GOTO 9999
        END IF

        !----------------------------------------------------
        ! Check if control parameters given are resonable.
        !----------------------------------------------------
        
        check_ctrl = &
             control_check_parameters(mcf_ctrl,stat_info_sub)
        
        IF( stat_info_sub /= 0 ) THEN
           PRINT *, "mcf : ", &
                "Checking control parameters has problem ! "
           stat_info = -1
           GOTO 9999
        END IF
        
        IF( .NOT. check_ctrl ) THEN
           PRINT *, "mcf : ", &
                "control parameters are not reasonable ! "
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Check if physics parameters given are resonable.
        !----------------------------------------------------
        
        check_physics= &
             physics_check_parameters(mcf_phys,stat_info_sub) 
        
        IF( stat_info_sub /= 0 ) THEN
           PRINT *, "mcf : ", &
                "Checking physics parameters has problem ! "
           stat_info = -1
           GOTO 9999
        END IF
        
        IF( .NOT. check_physics ) THEN
           PRINT *, "mcf : ", &
                "Physics parameters are not resonable !"
           stat_info = -1
           GOTO 9999
        END IF

        !----------------------------------------------------
        ! Check if io parameters given are reasonable.
        !----------------------------------------------------
        
        check_io = &
             io_check_parameters(mcf_io,stat_info_sub)
        
        IF(stat_info_sub /= 0 ) THEN
           PRINT *, "mcf : ", &
                "Check IO parameters has problem !" 
           stat_info = -1
           GOTO 9999
        END IF
        
        IF( .NOT. check_io ) THEN
           PRINT *, "mcf : ", &
                "IO parameters are not reasonable ! " 
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Get all the control variables which are read 
        ! from ctrl.mcf file.
        ! They decide the components of the simulation.
        !----------------------------------------------------
        
        debug_flag         = &
             control_get_debug_flag(mcf_ctrl,stat_info_sub)
        relax_run          = &
             control_get_relax_run(mcf_ctrl,stat_info_sub)
        read_external      = &
             control_get_read_external(mcf_ctrl,stat_info_sub)
        kernel_type        = &
             control_get_kernel_type(mcf_ctrl,stat_info_sub)
        symmetry           = &
             control_get_symmetry(mcf_ctrl,stat_info_sub)
        rhs_density_type   = &
             control_get_rhs_density_type(mcf_ctrl,stat_info_sub)
        stateEquation_type = &
             control_get_stateEquation_type(mcf_ctrl,stat_info_sub)
        Newtonian          = &
             control_get_Newtonian(mcf_ctrl,stat_info_sub)
        Brownian           = &
             control_get_Brownian(mcf_ctrl,stat_info_sub)
        random_seed        = &
             control_get_random_seed(mcf_ctrl,stat_info_sub)
        rhs_force_type     = &
             control_get_rhs_force_type(mcf_ctrl,stat_info_sub)
        p_energy           = &
             control_get_p_energy(mcf_ctrl,stat_info_sub)      
        integrate_type     = &
             control_get_integrate_type(mcf_ctrl,stat_info_sub)
        write_restart      = &
             control_get_write_restart(mcf_ctrl,stat_info_sub)
        
        !----------------------------------------------------
        ! Set the debug flag for global Debug object.
        !
        ! Note that each Class obejct uses the same 
        ! global debug object.
        !----------------------------------------------------
        
        CALL debug_new(global_debug,debug_flag,stat_info_sub)

#if __DEBUG

        !----------------------------------------------------
        ! Sychronize every process.
        ! Otherwise, the the time printed above
        ! will be a messed up with other printing afterwards
        ! among different processes.
        !----------------------------------------------------
        
        CALL MPI_Barrier(comm,stat_info_sub)
        
        !----------------------------------------------------
        ! For debug purpose,
        ! printing starting time of simulation.
        !----------------------------------------------------
        IF ( debug_flag == 2 ) THEN
           
           CALL debug_substart(global_debug,rank,&
                "mcf",time_start,stat_info_sub)
           
        END IF

        !----------------------------------------------------
        ! Sychronize every process.
        ! Otherwise, the the time printed above
        ! will be a messed up with other printing afterwards
        ! among different processes.
        !----------------------------------------------------
        
        CALL MPI_Barrier(comm,stat_info_sub)

#endif           
        
        !----------------------------------------------------
        ! For Brownian solvent, create a random 
        ! number generator with Gaussian distribution.
        !
        ! Note that the same random object will be
        ! used for generating random number for random force
        ! calculation, and also for random acceleration of 
        ! eigenvalues in case of non-Newtonian viscoelastic
        ! fluid. 
        !----------------------------------------------------
        
        CALL random_new(mcf_random,2,random_seed,stat_info_sub)
        IF(stat_info_sub /= 0 ) THEN
           PRINT *, "mcf : ", &
                "Creating random object has problem ! "
           stat_info = -1
           GOTO 9999
        END IF
        
        !----------------------------------------------------
        ! Get some basic physics parameters, 
        ! which are used for initializing kernel and PPM.
        !----------------------------------------------------
        
        num_dim = physics_get_num_dim(mcf_phys,stat_info_sub)
        cut_off = physics_get_cut_off(mcf_phys,stat_info_sub)
        num_colloid = physics_get_num_colloid(mcf_phys,stat_info_sub)

        !----------------------------------------------------
        ! Initialize the kernel object to set dimension, 
        ! type and cut off of kernel function.
        ! The kernel object gets kernel type and cut off,
        ! and calculate the smoothing length h accordingly.
        !
        ! Get the smoothing lenght h from kernel object.
        !----------------------------------------------------

        CALL kernel_new(mcf_kern,num_dim,kernel_type,&
             cut_off,stat_info_sub)
        h = kernel_get_h(mcf_kern,stat_info_sub)
        
        !----------------------------------------------------
        ! Set smoothing length in physics object,
        ! since it will be used for calculating dt.
        !----------------------------------------------------
        
        CAll physics_set_h(mcf_phys,h,stat_info_sub)
        
        !----------------------------------------------------
        ! Adjust and calculate physics parameters, such as dt.
        !----------------------------------------------------
        
        CALL physics_adjust_parameters(mcf_phys,stat_info_sub)
        
        IF( stat_info_sub /= 0 ) THEN
           PRINT *,&
                "mcf : Adjusting physics parameters has problem !"
           stat_info = -1
           GOTO 9999
        END IF
        
        CALL physics_get_min_phys_t(mcf_phys,min_phys_t,stat_info_sub)
        CALL physics_get_max_phys_t(mcf_phys,max_phys_t,stat_info_sub)
        CALL physics_get_bcdef(mcf_phys,bcdef,stat_info_sub)

        !----------------------------------------------------
        ! Adjusting IO parameters.
        !----------------------------------------------------
        
        CALL io_adjust_parameters(mcf_io,stat_info_sub)
        
        IF( stat_info_sub /= 0 ) THEN
           PRINT *, &
                "mcf : Adjusting io parameters has problem !"
           stat_info = -1
           GOTO 9999
        END IF

        !----------------------------------------------------
        ! Initialize PPM library, 
        ! which initializes topology also.
        !
        ! Note that min_phys_t and max_phys_t must be given
        ! instead of min_phys and max_phys.
        ! The former includes the solid wall in case of having
        ! solid wall, otherwise it is the same as latter one.
        !----------------------------------------------------

        CALL technique_init_ppm(mcf_tech,num_dim,&
             min_phys_t,max_phys_t,cut_off*1.0001_MK,&
             bcdef,stat_info_sub)
        
        !----------------------------------------------------
        ! Initialize the Rhs object
        ! decide which right hand side formulation,
        ! i.e. density,force, will be used.
        !----------------------------------------------------
        
        CALL rhs_new(mcf_rhs,mcf_ctrl,mcf_phys,&
             mcf_random,stat_info_sub)
        
        !----------------------------------------------------
        ! Initialize the StateEquation object to decide 
        ! which type of state equation will be used.
        !----------------------------------------------------
        
        rho     = physics_get_rho(mcf_phys,stat_info_sub)
        c       = physics_get_c(mcf_phys,stat_info_sub)
        rho_ref = physics_get_rho_ref(mcf_phys,stat_info_sub)
        gamma   = physics_get_gamma(mcf_phys,stat_info_sub)
        
        CALL stateEquation_new(mcf_stateEquation,&
             stateEquation_type,&
             c,rho,rho_ref,gamma,stat_info_sub)
        
        !----------------------------------------------------
        ! Initialize the Particles object,
        ! by importing various objects which
        ! will be used by particles object.
        !----------------------------------------------------

        CALL particles_new(mcf_particles,mcf_ctrl,mcf_phys,&
             mcf_rhs,mcf_stateEquation,mcf_kern,mcf_tech,&
             stat_info_sub)
        
        
        !----------------------------------------------------
        ! Intialize the Marching object,
        ! set the integration scheme.
        !----------------------------------------------------
        
        CALL marching_new(mcf_marching,&
             mcf_io,mcf_particles,stat_info_sub)
        
        
        !----------------------------------------------------
        ! Display parameters on root.
        !----------------------------------------------------
        
        IF ( rank == 0 ) THEN
           
           CALL control_display_parameters(mcf_ctrl,stat_info_sub)
           CALL io_display_parameters(mcf_io,stat_info_sub)
           !CALL kernel_display_parameters(mcf_kern,stat_info_sub)
           !CALL physics_display_parameters(mcf_phys,stat_info_sub)
           !CALL debug_display_parameters(global_debug,stat_info_sub)
           !CALL rhs_display_parameters(mcf_rhs,stat_info_sub)
           CALL stateEquation_display_parameters(mcf_stateEquation,&
           stat_info_sub)
           CALL technique_display_parameters(mcf_tech,stat_info_sub)
           !CALL marching_display_parameters(mcf_marching,&
           !     stat_info_sub)
           
        END IF
        
        !----------------------------------------------------
        ! Sychronize every process.
        !----------------------------------------------------
        
        CALL MPI_Barrier(comm,stat_info_sub)
      
        IF ( relax_run .OR. Brownian ) THEN
           CALL random_display_parameters(mcf_random,stat_info_sub)
        END IF
        !----------------------------------------------------
        ! Root process reads the particles position, 
        ! velocity, density, mass, particle id, 
        ! species ID from file externally,
        ! or generates them on certain lattice internally.
        !----------------------------------------------------
        
        IF ( rank == 0 ) THEN
           
           IF( read_external ) THEN
              
              !----------------------------------------------
              ! Read particles externally,
              ! either from preprocessing or 
              ! restart file of last run.
              !
              ! x,y(, z), vx,vy(, vz), 
              ! rho, mass, p_id, s_id.
              !
              ! Note that rho can be physical density or
              ! number density, depending on rhs_density_type
              ! chosen.
              !----------------------------------------------
              
              CALL io_read_particles(mcf_io,rank,&
                   mcf_particles,stat_info_sub)
              
              !----------------------------------------------
              ! In case of non-Newtonian viscoelastic
              ! oldroyd-B model fluid, read conformation 
              ! tensor:
              !
              !  cxx, cxy, (cxz),
              !  cyx, cyy, (cyz),
              ! (czx, czy,  czz)
              !----------------------------------------------
              
              IF ( .NOT. Newtonian ) THEN
                 
                 CALL io_read_conformation(mcf_io,rank,&
                      mcf_particles,stat_info_sub)
                 
              END IF
              
           ELSE ! not read_external, generate internally
              
              !----------------------------------------------
              ! Generate particles internally 
              ! using certain lattice.
              !
              ! x, y(, z),
              ! vx,vy(, vz),
              ! p_id, s_id of global domain.
              !----------------------------------------------
              
              CALL particles_init_global_inter(mcf_particles,&
                   rank,stat_info_sub)
              CALL particles_init_global_assign_id(mcf_particles,&
                   rank,stat_info_sub)
              
           END IF ! read_external
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "mcf : ", &
                   "Initializing particles has problem !"
              stat_info = -1
              GOTO 9999
           END IF
           
           !-------------------------------------------------
           ! Display particles; parameters on root
           ! before decompsing.
           !-------------------------------------------------
           
           CALL particles_display_parameters(mcf_particles,&
                0,stat_info_sub)
           
        END IF ! rank == 0
        
        !----------------------------------------------------
        ! Decompose the domain of all the particles,
        ! broadcast them from root to all processes.
        !
        ! If particles are read from files, that means
        ! root process has position, velocity, density,
        ! mass and p_ID, s_ID, we have to broadcast all.
        !
        ! Otherwise particles are created internally,
        ! global decomposition only for position, 
        ! velocity and p_ID, s_ID.
        ! After broadcast, allocate memory for 
        ! other quanties on local process, 
        ! such as density and mass.
        !----------------------------------------------------
        
        
        IF ( read_external ) THEN
           
           !-------------------------------------------------
           ! x,y(, z); vx,vy(, vz);
           ! rho; mass; p_ID,s_ID;
           ! for non-Newtonian viscoelastic Oldroyd-B model
           ! also conformation tensor,
           !  cxx, cxy, (cxz),
           !  cyx, cyy, (cyz),
           ! (czx, czy,  czz),
           ! which are all read globally on root process,
           ! we decompose them onto all processes.
           !-------------------------------------------------
           
           CALL particles_decompose_global(mcf_particles,&
                l_map_x   = .TRUE., l_map_v = .TRUE., &
                l_map_rho = .TRUE., l_map_m = .TRUE., &
                l_map_id  = .TRUE.,  &
                l_map_ct  = (.NOT. Newtonian), &
                stat_info = stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *,  "mcf : ", &
                   "particles_decompose_global all failed !"
              stat_info = -1
              GOTO 9999
           END IF
           
           !-------------------------------------------------
           ! Create the rest quntities, besides,
           ! x, y, (z),
           ! vx,vy,(vz),
           ! rho,  mass,
           ! p_ID, s_ID on ach process.
           !
           ! In case of non-Newtonian fluid,
           ! conformation tensor C also has been
           ! read already, is not needed to create.
           !-------------------------------------------------
           
           CALL particles_init_partial_exter(mcf_particles,&
                stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "mcf : ", &
                   "Generating quntities locally failed !"
              stat_info = -1
              GOTO 9999
           END IF
           
        ELSE ! particle created internall, then distribute.
           
           !-------------------------------------------------
           ! x, y(, z); vx, vy(, vz);  p_ID, s_ID
           ! are generated globally on root process(rank=0),
           ! we decompose them onto all processes.
           !-------------------------------------------------
           
           CALL particles_decompose_global(mcf_particles,&
                l_map_x  = .TRUE., l_map_v=.TRUE., &
                l_map_id = .TRUE., stat_info=stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *,  "mcf : ", &
                   "Decomposing x, v, IDs globally failed !"
              stat_info = -1
              GOTO 9999
           END IF
           
           !-------------------------------------------------
           ! Create the rest of quntities, besides,
           !  x,  y(,  z);
           ! vx, vy(, vz);
           ! p_ID, s_ID on ach process.
           !-------------------------------------------------
           
           CALL particles_init_partial_inter(mcf_particles,&
                stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "mcf : ", &
                   "Generating quantities locally failed !"
              stat_info = -1
              GOTO 9999
           END IF
           
           !-------------------------------------------------
           ! Calculate particles' mass initially,
           ! from density.
           !-------------------------------------------------
           
           CALL particles_compute_mass(mcf_particles,stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "mcf : ", &
                   "Computing mass failed !"
              stat_info = -1
              GOTO 9999
           END IF
           
        END IF ! read_external
        
        !----------------------------------------------------
        ! Display physics parameters on root process after 
        ! decomposing, since for internally generated particles, 
        ! e.g., mass of particles, mass of colloid have to be
        ! calculated, therefore we display them here.
        !----------------------------------------------------
        
        IF ( rank == 0 ) THEN
           
           CALL physics_display_parameters(mcf_phys,stat_info_sub)
           
        END IF
        
        
        !----------------------------------------------------
        ! Sychronize every process, since we want to display
        ! parameters on each process simutaneously.
        !----------------------------------------------------
        
        CALL MPI_Barrier(comm,stat_info_sub)
        
        !----------------------------------------------------
        ! Display particles parameters on each rank.
        !----------------------------------------------------
        
        CALL particles_display_parameters(mcf_particles,&
             1,stat_info_sub)
        
        !----------------------------------------------------
        ! Sychronize every process.
        !----------------------------------------------------
        
        CALL MPI_Barrier(comm,stat_info_sub)
        
        !----------------------------------------------------
        ! Relax particles before really run simulation.
        !----------------------------------------------------
        
        IF ( relax_run ) THEN
           
           CALL marching_relax(mcf_marching,stat_info_sub)
           
           IF( stat_info_sub /=0 ) THEN
              PRINT *, "mcf : ", &
                   "Marching relax failed!"
              stat_info = -1
              GOTO 9999
           END IF

        END IF
        
        !----------------------------------------------------
        ! Marching the simulation, 
        ! i.e., integrating with time.
        !----------------------------------------------------

        CALL marching_marching(mcf_marching,stat_info_sub)
        
        IF( stat_info_sub /=0 ) THEN
           PRINT *, "mcf : ", &
                "Marching integration failed!"
           stat_info = -1
           GOTO 9999
        END IF
        
        
9999    CONTINUE
        
        !----------------------------------------------------
        ! Finalization.
        !----------------------------------------------------
        
        IF(ASSOCIATED(min_phys_t))THEN
           DEALLOCATE(min_phys_t)
        END IF
        
        IF(ASSOCIATED(max_phys_t))THEN
           DEALLOCATE(max_phys_t)
        END IF
        
        IF(ASSOCIATED(bcdef))THEN
           DEALLOCATE(bcdef)
        END IF
        
        CALL particles_finalize(mcf_particles,stat_info_sub)
        CALL control_finalize(mcf_ctrl,stat_info_sub)
        CALL physics_finalize(mcf_phys,stat_info_sub)
        IF ( rank == 0 ) THEN
           CALL io_finalize(mcf_io,rank,num_colloid,stat_info_sub)
        END IF
        CALL rhs_finalize(mcf_rhs,stat_info_sub)
        CALL stateEquation_finalize(mcf_stateEquation,stat_info_sub)
        CALL kernel_finalize(mcf_kern,stat_info_sub)
        CALL marching_finalize(mcf_marching,stat_info_sub)
        
        IF ( stat_info == 0 ) THEN
           CALL MPI_Barrier(comm,stat_info_sub)
        END IF

        CALL technique_finalize(mcf_tech,stat_info,stat_info_sub)
        
        
#if __DEBUG
        
        !----------------------------------------------------
        ! For debug purpose,
        ! printing the fininish time of simulation.
        !----------------------------------------------------
        IF ( debug_flag == 2 ) THEN
           
           CALL debug_substop(global_debug,rank,&
                "mcf",time_start,stat_info_sub)
           
        END IF
        
        CALL debug_finalize(global_debug,stat_info_sub)
        
#endif         
        
        
      END PROGRAM mcf
      
#if 0
      !------------------------------------------------------
      ! Terminology. 
      !------------------------------------------------------
      !
      ! host particles : the ones whose positions are used for
      !                  generating ghosts particles.
      !                  the original particle is called the
      !                  host particle of the ghost particle.
      !
      !------------------------------------------------------
#endif
      

      
