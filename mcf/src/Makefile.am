AM_FCFLAGS   = -O2 -D__MPI -I ${top_srcdir}/src
bin_PROGRAMS = mcf
mcf_SOURCES  = mcf.F90\
	       mcf_header.F90\
	       Class_Tool.F90\
	       Class_Random.F90\
	       Class_Debug.F90\
	       Class_Control.F90\
               Class_Boundary.F90\
               Class_Colloid.F90\
	       Class_Physics.F90\
               Class_Rhs.F90\
	       Class_IO.F90\
	       Class_StateEquation.F90\
	       Class_Kernel.F90\
               Class_Technique.F90\
	       Class_Particles.F90\
	       Class_Statistic.F90\
	       Class_Marching.F90

EXTRA_DIST   = \
	boundary_check_wall_solid_particle.F90 boundary_collect_colloid_interaction.F90 boundary_collect_particles_interaction.F90 \
	boundary_finalize.F90 boundary_get.F90 boundary_new.F90 \
	boundary_noslip_mirror.F90 boundary_noslip_solid.F90 boundary_set.F90 \
	boundary_update_boundary.F90 colloid_adjust_colloid.F90 colloid_adjust_parameters.F90 \
	colloid_apply_body_force.F90 colloid_cartesian_ellipse.F90 colloid_cartesian_ellipsoid.F90 \
	colloid_check_boundary_particle.F90 colloid_check_parameters.F90 colloid_check_penetration.F90 \
	colloid_collect_particles_interaction.F90 colloid_compute_acceleration.F90 colloid_compute_accumulation_matrix.F90 \
	colloid_compute_accumulation_vector.F90 colloid_compute_dt_f.F90 colloid_compute_image.F90 \
	colloid_compute_interaction.F90 colloid_compute_lubrication_cc.F90 colloid_compute_lubrication_cw.F90 \
	colloid_compute_repulsion_cc.F90 colloid_compute_repulsion_cw.F90 colloid_compute_rotation_matrix.F90 \
	colloid_compute_rotation_vector.F90 colloid_compute_statistic.F90 colloid_create_boundary_particle_2D.F90 \
	colloid_create_boundary_particle_3D.F90 colloid_create_boundary_particle_3D_dicolloid.F90 colloid_create_boundary_particle_3D_sphere.F90 \
	colloid_finalize.F90 colloid_find_force_extreme.F90 colloid_get.F90 \
	colloid_in_nearest_image.F90 colloid_in_relative_position.F90 colloid_init_accumulation_matrix.F90 \
	colloid_initialize_image.F90 colloid_integrate_position.F90 colloid_integrate_velocity.F90 \
	colloid_nearest_image.F90 colloid_new.F90 colloid_noslip.F90 \
	colloid_noslip_Morris.F90 colloid_noslip_Morris_cylinder_2D.F90 colloid_noslip_Morris_cylinder_3D.F90 \
	colloid_noslip_Morris_dicolloid_3D.F90 colloid_noslip_Morris_disk.F90 colloid_noslip_Morris_ellipse.F90 \
	colloid_noslip_Morris_ellipsoid.F90 colloid_noslip_Morris_sphere.F90 colloid_noslip_Morris_star.F90 \
	colloid_noslip_Morris_star_2D.F90 colloid_noslip_Zhu.F90 colloid_noslip_frozen.F90 \
	colloid_particle_velocity.F90 colloid_polar_ellipse.F90 colloid_polar_star.F90 \
	colloid_set.F90 colloid_set_flow_developed.F90 colloid_spherical_ellipsoid.F90 \
	control_check_parameters.F90 control_finalize.F90 control_get.F90 \
	control_new.F90 control_set.F90 debug_close.F90 \
	debug_finalize.F90 debug_get.F90 debug_new.F90 \
	debug_open.F90 debug_print_msg.F90 debug_set.F90 \
	debug_substart.F90 debug_substop.F90 debug_validate_motion.F90 \
	debug_write_output.F90 debug_write_time.F90 io_adjust_parameters.F90 \
	io_check_parameters.F90 io_close.F90 io_finalize.F90 \
	io_get.F90 io_new.F90 io_open.F90 \
	io_read.F90 io_read_conformation.F90 io_read_ctrl.F90 \
	io_read_io.F90 io_read_particles.F90 io_read_physics.F90 \
	io_write.F90 io_write_boundary.F90 io_write_colloid.F90 \
	io_write_colloid_separate.F90 io_write_condition.F90 io_write_conformation.F90 \
	io_write_particles.F90 io_write_particles_relax.F90 io_write_restart_conformation.F90 \
	io_write_restart_particles.F90 io_write_restart_particles_relax.F90 io_write_restart_physics.F90 \
	io_write_statistic.F90 io_write_statistic_relax.F90 kernel_finalize.F90 \
	kernel_get.F90 kernel_kernel.F90 kernel_kernel_Lucy.F90 \
	kernel_kernel_quintic_spline.F90 kernel_new.F90 marching_adjust_flow_v.F90 \
	marching_finalize.F90 marching_integrate.F90 marching_integrate_Euler.F90 \
	marching_integrate_VV.F90 marching_marching.F90 marching_new.F90 \
	marching_relax.F90 mcf.F90 mcf_header.F90 \
	particles_adjust_particles.F90 particles_apply_body_force.F90 particles_collect_boundary_interaction.F90 \
	particles_collect_colloid_interaction.F90 particles_compute_act.F90 particles_compute_aeval.F90 \
	particles_compute_aevec.F90 particles_compute_colloid_absolute_position.F90 particles_compute_colloid_relative_position.F90 \
	particles_compute_ct.F90 particles_compute_density.F90 particles_compute_dt_f.F90 \
	particles_compute_evgt.F90 particles_compute_interaction.F90 particles_compute_mass.F90 \
	particles_compute_pressure.F90 particles_compute_pressure_tensor.F90 particles_decompose_global.F90 \
	particles_decompose_partial.F90 particles_finalize.F90 particles_find_density_extreme.F90 \
	particles_find_force_extreme.F90 particles_get.F90 particles_init_global_assign_id.F90 \
	particles_init_global_exter.F90 particles_init_global_inter.F90 particles_init_global_inter_cubic.F90 \
	particles_init_global_inter_hexagonal.F90 particles_init_global_inter_square.F90 particles_init_global_inter_staggered.F90 \
	particles_init_partial_exter.F90 particles_init_partial_inter.F90 particles_integrate_boundary_position.F90 \
	particles_integrate_ct.F90 particles_integrate_eval.F90 particles_integrate_evec.F90 \
	particles_integrate_position.F90 particles_integrate_potential_energy.F90 particles_integrate_velocity.F90 \
	particles_map_ghost_get.F90 particles_map_ghost_put.F90 particles_new.F90 \
	particles_normalize_density.F90 particles_reset.F90 particles_reset_boundary_ghost_interaction.F90 \
	particles_reset_boundary_interaction.F90 particles_reset_boundary_velocity.F90 particles_reset_colloid_interaction.F90 \
	particles_reset_colloid_velocity.F90 particles_set.F90 particles_set_boundary_ghost_id.F90 \
	particles_set_boundary_ghost_velocity.F90 particles_set_boundary_velocity.F90 particles_set_colloid_on_lattice.F90 \
	particles_set_colloid_velocity.F90 particles_set_flow_developed.F90 physics_adapt_dt.F90 \
	physics_adjust_parameters.F90 physics_check_parameters.F90 physics_finalize.F90 \
	physics_get.F90 physics_initialize_dt.F90 physics_new.F90 \
	physics_set.F90 random_new.F90 random_random.F90 \
	random_random_Gaussian.F90 random_random_uniform.F90 rhs_density.F90 \
	rhs_finalize.F90 rhs_force_cc_Newtonian.F90 rhs_force_ff_Newtonian.F90 \
	rhs_force_ff_Newtonian_Espanol.F90 rhs_force_ff_Newtonian_HuAdams.F90 rhs_force_ff_Newtonian_HuAdams_angular.F90 \
	rhs_force_ff_Newtonian_Morris.F90 rhs_force_ff_non_Newtonian.F90 rhs_force_ff_non_Newtonian_Espanol.F90 \
	rhs_force_ff_non_Newtonian_HuAdams.F90 rhs_force_ff_non_Newtonian_HuAdams_angular.F90 rhs_get.F90 \
	rhs_new.F90 rhs_set.F90 stateEquation_compute_pressure.F90 \
	stateEquation_finalize.F90 stateEquation_get.F90 stateEquation_new.F90 \
	stateEquation_set.F90 statistic_compute_disorder.F90 statistic_compute_statistic.F90 \
	statistic_compute_v_average.F90 statistic_finalize.F90 statistic_get.F90 \
	statistic_new.F90 statistic_set.F90 technique_build_list.F90 \
	technique_finalize.F90 technique_get.F90 technique_new.F90 \
	tool_cross_product.F90 tool_new.F90 tool_rotation_matrix.F90 \
	tool_rotation_vector.F90 tool_uppercase.F90\
	pp_interaction.inc pp_interaction_cc.inc pp_interaction_cf.inc \
	pp_interaction_cw.inc pp_interaction_fc.inc pp_interaction_ff.inc \
	pp_interaction_fw.inc pp_interaction_wc.inc pp_interaction_wf.inc \
	pp_vgt_ip.inc pp_vgt_jp.inc ppm_param.inc

Makefile.dep:
	$(MAKEDEPF90) -I ${top_srcdir}/src $(mcf_SOURCES) > $@

clean-local:
	-rm *.mod Makefile.dep

-include Makefile.dep



