if HAVE_MPI
lib_LIBRARIES = libppm.a
libppm_a_SOURCES=
else
lib_LIBRARIES = libppm.a  libfakempi.a
libfakempi_a_SOURCES= m_fakempi.F90
libppm_a_LIBADD=libfakempi.a
libppm_a_SOURCES=fakempi.h
endif

EXTRA_DIST=fakempi.h m_fakempi.F90

libppm_a_SOURCES += \
	m_vizing.F90  \
	ppm_module_alloc.F90 ppm_module_bem_basis.F90 ppm_module_bem.F90\
	ppm_module_bem_quadrule.F90 ppm_module_bem_quadrule_npoints.F90 ppm_module_check_meshid.F90\
	ppm_module_check_topoid.F90 ppm_module_clist_destroy.F90 ppm_module_comp_mesh.F90\
	ppm_module_comp_mesh_fd.F90 ppm_module_comp_mesh_fft.F90 ppm_module_comp_part.F90\
	ppm_module_comp_pp_cell.F90 ppm_module_comp_pp_correct.F90 ppm_module_comp_pp_doring.F90\
	ppm_module_comp_pp_mk_table.F90 ppm_module_comp_pp_ring.F90 ppm_module_comp_pp_verlet.F90\
	ppm_module_copy_image_subs.F90 ppm_module_data.F90 ppm_module_data_fieldsolver.F90\
	ppm_module_data_fmm.F90 ppm_module_data_io.F90 ppm_module_data_loadbal.F90\
	ppm_module_data_mesh.F90 ppm_module_data_mg.F90 ppm_module_data_neighlist.F90\
	ppm_module_data_ode.F90 ppm_module_data_part.F90 ppm_module_data_rmsh.F90\
	ppm_module_data_tree.F90 ppm_module_decomp_boxsplit.F90 ppm_module_decomp_cartesian.F90\
	ppm_module_decomp_null.F90 ppm_module_decomp_pruned_cell.F90 ppm_module_decomp_tree.F90\
	ppm_module_define_subs_bc.F90 ppm_module_doio.F90 ppm_module_error.F90\
	ppm_module_estimate_proc_speed.F90 ppm_module_fdsolver_fft_bd.F90 ppm_module_fdsolver_fft_fd.F90\
	ppm_module_fdsolver_finalize.F90 ppm_module_fdsolver_init.F90 ppm_module_fdsolver_map.F90\
	ppm_module_fdsolver_poisson.F90 ppm_module_fdsolver_solve.F90 ppm_module_fieldsolver.F90\
	ppm_module_finalize.F90 ppm_module_find_duplicates.F90 ppm_module_find_neigh.F90\
	ppm_module_fmm_expansion.F90 ppm_module_fmm_expchange.F90 ppm_module_fmm.F90\
	ppm_module_fmm_finalize.F90 ppm_module_fmm_init.F90 ppm_module_fmm_potential.F90\
	ppm_module_fmm_pretraverse.F90 ppm_module_fmm_traverse.F90 ppm_module_get_cost.F90\
	ppm_module_get_revision.F90 ppm_module_impose_part_bc.F90 ppm_module_init.F90\
	ppm_module_interp_m2p.F90 ppm_module_interp_p2m.F90 ppm_module_io_close.F90\
	ppm_module_io_delete.F90 ppm_module_io.F90 ppm_module_io_inquire.F90\
	ppm_module_io_open.F90 ppm_module_io_read_ascii.F90 ppm_module_io_read_binary.F90\
	ppm_module_io_set_unit.F90 ppm_module_io_unused_unit.F90 ppm_module_io_write_ascii.F90\
	ppm_module_io_write_binary.F90 ppm_module_is_initialized.F90 ppm_module_loadbal.F90\
	ppm_module_loadbal_inquire.F90 ppm_module_log.F90 ppm_module_map_connect_distrib.F90\
	ppm_module_map_connect.F90 ppm_module_map_connect_prune.F90 ppm_module_map_connect_send.F90\
	ppm_module_map.F90 ppm_module_map_field.F90 ppm_module_map_field_ghost.F90\
	ppm_module_map_field_ghost_get.F90 ppm_module_map_field_ghost_init.F90 ppm_module_map_field_ghost_put.F90\
	ppm_module_map_field_global.F90 ppm_module_map_field_pop_2d.F90 ppm_module_map_field_pop_3d.F90\
	ppm_module_map_field_push_2d.F90 ppm_module_map_field_push_3d.F90 ppm_module_map_field_send.F90\
	ppm_module_map_part.F90 ppm_module_map_part_get_sub.F90 ppm_module_map_part_ghost.F90\
	ppm_module_map_part_ghost_get.F90 ppm_module_map_part_ghost_pop.F90 ppm_module_map_part_ghost_put.F90\
	ppm_module_map_part_global.F90 ppm_module_map_part_partial.F90 ppm_module_map_part_pop.F90\
	ppm_module_map_part_push.F90 ppm_module_map_part_remap.F90 ppm_module_map_part_ring.F90\
	ppm_module_map_part_ring_shift.F90 ppm_module_map_part_send.F90 ppm_module_mesh_alloc.F90\
	ppm_module_mesh_block_intersect.F90 ppm_module_mesh_define.F90 ppm_module_mesh_derive.F90\
	ppm_module_mesh_finalize.F90 ppm_module_mesh_on_subs.F90 ppm_module_mesh_store.F90\
	ppm_module_mg_alloc.F90 ppm_module_mg_core.F90 ppm_module_mg.F90\
	ppm_module_mg_finalize.F90 ppm_module_mg_init.F90 ppm_module_mg_prolong.F90\
	ppm_module_mg_res.F90 ppm_module_mg_restrict.F90 ppm_module_mg_smooth.F90\
	ppm_module_mg_solv.F90 ppm_module_mktopo.F90 ppm_module_neighlist_clist.F90\
	ppm_module_neighlist.F90 ppm_module_neighlist_MkNeighIdx.F90 ppm_module_neighlist_vlist.F90\
	ppm_module_ode_alldone.F90 ppm_module_ode_create_ode.F90 ppm_module_ode.F90\
	ppm_module_ode_finalize.F90 ppm_module_ode_init.F90 ppm_module_ode_map_pop.F90\
	ppm_module_ode_map_push.F90 ppm_module_ode_start.F90 ppm_module_ode_step.F90\
	ppm_module_print_defines.F90 ppm_module_rmsh_comp_weights.F90 ppm_module_rmsh_create_part.F90\
	ppm_module_rmsh.F90 ppm_module_rmsh_remesh.F90 ppm_module_scale_domain.F90\
	ppm_module_set_decomp_cost.F90 ppm_module_set_proc_speed.F90 ppm_module_substart.F90\
	ppm_module_substop.F90 ppm_module_time.F90 ppm_module_topo_box2subs.F90\
	ppm_module_topo_check.F90 ppm_module_topo_cost.F90 ppm_module_topo.F90\
	ppm_module_topo_inquire.F90 ppm_module_topo_metis_s2p.F90 ppm_module_topo_store.F90\
	ppm_module_topo_subs2proc.F90 ppm_module_tree_alloc.F90 ppm_module_tree_boxcost.F90\
	ppm_module_tree_boxcut.F90 ppm_module_tree_cutdir.F90 ppm_module_tree_cutpos.F90\
	ppm_module_tree_divcheck.F90 ppm_module_tree_done.F90 ppm_module_tree.F90\
	ppm_module_user.F90 ppm_module_user_io.F90 ppm_module_user_util.F90\
	ppm_module_util_cart2sph.F90 ppm_module_util_collapse_list.F90 ppm_module_util_commopt.F90\
	ppm_module_util_cubeq_real.F90 ppm_module_util_eigen_2sym.F90 ppm_module_util_eigen_3sym.F90\
	ppm_module_util_fft_backward.F90 ppm_module_util_fft_forward.F90 ppm_module_util_gmres.F90\
	ppm_module_util_invert_list.F90 ppm_module_util_matinv_2x2.F90 ppm_module_util_matinv_3x3.F90\
	ppm_module_util_qsort.F90 ppm_module_util_quadeq_real.F90 ppm_module_util_rank.F90\
	ppm_module_util_sort.F90 ppm_module_util_time.F90 ppm_module_write.F90

dist_noinst_DATA=\
	ppm_define.h  ppm_error.h  ppm_ode_modalloc.h  ppm_ode_rhsfunc_macro.h  ppm_param.h\
	ppm_alloc_1d.F90 ppm_alloc_1dl.F90 ppm_alloc_2d.F90\
	ppm_alloc_2dl.F90 ppm_alloc_3d.F90 ppm_alloc_3dl.F90\
	ppm_alloc_4d.F90 ppm_alloc_4dl.F90 ppm_alloc_5d.F90\
	ppm_alloc_5dl.F90 ppm_bem_basis.F90 ppm_bem_quadrule.F90\
	ppm_bem_quadrule_npoints.F90 ppm_check_meshid.F90 ppm_check_topoid.F90\
	ppm_clist_destroy.F90 ppm_comp_pp_cell.F90 ppm_comp_pp_correct.F90\
	ppm_comp_pp_doring.F90 ppm_comp_pp_mk_table.F90 ppm_comp_pp_ring.F90\
	ppm_comp_pp_verlet.F90 ppm_copy_image_subs.F90 ppm_decomp_boxsplit.F90\
	ppm_decomp_cartesian.F90 ppm_decomp_null.F90 ppm_decomp_pruned_cell.F90\
	ppm_decomp_tree.F90 ppm_define_subs_bc.F90 ppm_doio.F90\
	ppm_error.F90 ppm_estimate_proc_speed.F90 ppm_fdsolver_fft_bd_2d.F90\
	ppm_fdsolver_fft_bd_3d.F90 ppm_fdsolver_fft_fd_2d.F90 ppm_fdsolver_fft_fd_3d.F90\
	ppm_fdsolver_finalize.F90 ppm_fdsolver_init.F90 ppm_fdsolver_map_2d.F90\
	ppm_fdsolver_map_3d.F90 ppm_fdsolver_poisson_2d.F90 ppm_fdsolver_poisson_3d.F90\
	ppm_fdsolver_solve_2d.F90 ppm_fdsolver_solve_3d.F90 ppm_finalize.F90\
	ppm_find_duplicates.F90 ppm_find_neigh.F90 ppm_fmm_expansion.F90\
	ppm_fmm_expchange.F90 ppm_fmm_finalize.F90 ppm_fmm_init.F90\
	ppm_fmm_potential.F90 ppm_fmm_pretraverse.F90 ppm_fmm_traverse.F90\
	ppm_get_cost.F90 ppm_get_revision.F90 ppm_impose_part_bc.F90\
	ppm_init.F90 ppm_interp_m2p.F90 ppm_interp_p2m.F90\
	ppm_io.F90 ppm_io_close.F90 ppm_io_delete.F90\
	ppm_io_inquire.F90 ppm_io_open.F90 ppm_io_read_ascii.F90\
	ppm_io_read_binary.F90 ppm_io_set_unit.F90 ppm_io_unused_unit.F90\
	ppm_io_write_ascii.F90 ppm_io_write_binary.F90 ppm_is_initialized.F90\
	ppm_loadbal_inquire.F90 ppm_log.F90 ppm_map_connect.F90\
	ppm_map_connect_distrib.F90 ppm_map_connect_prune.F90 ppm_map_connect_send.F90\
	ppm_map_field.F90 ppm_map_field_ghost.F90 ppm_map_field_ghost_get.F90\
	ppm_map_field_ghost_init.F90 ppm_map_field_ghost_put.F90 ppm_map_field_global.F90\
	ppm_map_field_pop_2d.F90 ppm_map_field_pop_3d.F90 ppm_map_field_push_2d.F90\
	ppm_map_field_push_3d.F90 ppm_map_field_send.F90 ppm_map_part.F90\
	ppm_map_part_get_sub.F90 ppm_map_part_ghost.F90 ppm_map_part_ghost_get.F90\
	ppm_map_part_ghost_pop.F90 ppm_map_part_ghost_put.F90 ppm_map_part_global.F90\
	ppm_map_part_partial.F90 ppm_map_part_pop.F90 ppm_map_part_push.F90\
	ppm_map_part_remap.F90 ppm_map_part_ring.F90 ppm_map_part_ring_shift.F90\
	ppm_map_part_send.F90 ppm_mesh_alloc_equi.F90 ppm_mesh_alloc_list.F90\
	ppm_mesh_block_intersect.F90 ppm_mesh_define.F90 ppm_mesh_derive.F90\
	ppm_mesh_finalize.F90 ppm_mesh_on_subs.F90 ppm_mesh_store.F90\
	ppm_mg_alloc_bc.F90 ppm_mg_alloc_field.F90 ppm_mg_core.F90\
	ppm_mg_finalize.F90 ppm_mg_init.F90 ppm_mg_prolong.F90\
	ppm_mg_res_coarse.F90 ppm_mg_res_fine.F90 ppm_mg_restrict.F90\
	ppm_mg_smooth_coarse.F90 ppm_mg_smooth_fine.F90 ppm_mg_solv.F90\
	ppm_neighlist_MkNeighIdx.F90 ppm_neighlist_clist.F90 ppm_neighlist_vlist.F90\
	ppm_ode_alldone.F90 ppm_ode_create_ode.F90 ppm_ode_finalize.F90\
	ppm_ode_init.F90 ppm_ode_map_pop.F90 ppm_ode_map_push.F90\
	ppm_ode_start.F90 ppm_ode_step.F90 ppm_print_defines.F90\
	ppm_rmsh_comp_weights.F90 ppm_rmsh_create_part.F90 ppm_rmsh_remesh.F90\
	ppm_scale_domain.F90 ppm_set_decomp_cost.F90 ppm_set_proc_speed.F90\
	ppm_template_comp_boundary.F90 ppm_template_comp_pp_cell.F90 ppm_template_comp_pp_ring.F90\
	ppm_template_comp_pp_verlet.F90 ppm_template_connect_interaction.F90 ppm_template_ring_interaction.F90\
	ppm_time.F90 ppm_topo_box2subs.F90 ppm_topo_check.F90\
	ppm_topo_cost.F90 ppm_topo_inquire.F90 ppm_topo_metis_s2p.F90\
	ppm_topo_mkfield.F90 ppm_topo_mkgeom.F90 ppm_topo_mkpart.F90\
	ppm_topo_mktree.F90 ppm_topo_store.F90 ppm_topo_subs2proc.F90\
	ppm_tree.F90 ppm_tree_alloc.F90 ppm_tree_boxcost.F90\
	ppm_tree_boxcut.F90 ppm_tree_cutdir.F90 ppm_tree_cutpos.F90\
	ppm_tree_divcheck.F90 ppm_tree_done.F90 ppm_util_cart2sph.F90\
	ppm_util_collapse_list.F90 ppm_util_commopt.F90 ppm_util_cubeq_real.F90\
	ppm_util_eigen_2sym.F90 ppm_util_eigen_3sym.F90 ppm_util_fft_backward_2d.F90\
	ppm_util_fft_backward_3d.F90 ppm_util_fft_forward_2d.F90 ppm_util_fft_forward_3d.F90\
	ppm_util_gmres.F90 ppm_util_gmres_solveupper.F90 ppm_util_invert_list.F90\
	ppm_util_matinv_2x2.F90 ppm_util_matinv_3x3.F90 ppm_util_qsort.F90\
	ppm_util_quadeq_real.F90 ppm_util_rank2d.F90 ppm_util_rank3d.F90\
	ppm_util_sort2d.F90 ppm_util_sort3d.F90 ppm_util_time.F90\
	ppm_write.F90 substart.F90 substop.F90\
	ppm_bcdef_neumann_lower_x.inc ppm_bcdef_neumann_lower_y.inc ppm_bcdef_neumann_lower_z.inc \
	ppm_bcdef_neumann_upper_x.inc ppm_bcdef_neumann_upper_y.inc ppm_bcdef_neumann_upper_z.inc \
	ppm_bcdef_periodic_x.inc ppm_bcdef_periodic_y.inc \
	ppm_bcdef_periodic_z.inc ppm_bcdef_symmetry_lower_x.inc ppm_bcdef_symmetry_lower_y.inc \
	ppm_bcdef_symmetry_lower_z.inc ppm_bcdef_symmetry_upper_x.inc ppm_bcdef_symmetry_upper_y.inc \
	ppm_bcdef_symmetry_upper_z.inc ppm_comp_pp_kernels.inc ppm_find_neigh_subs_2d.inc \
	ppm_find_neigh_subs_3d.inc \
	ppm_bcdef_le_x.inc ppm_bcdef_le_y.inc ppm_bcdef_le_z.inc\
	Makefile.dep

install-data-local: 	
	@for p in *; do \
	  case "$$p" in \
	    *.mod) \
		echo " $(INSTALL_HEADER) '$$p' '$(includedir)'"; \
		$(INSTALL_HEADER) "$$p" $(includedir)/;; \
	  esac \
	done


-include Makefile.dep

MAKEDEPF90_FLAGS=-DHAVE_CONFIG_H
Makefile.dep:
	$(MAKEDEPF90) $(MAKEDEPF90_FLAGS) -I ${top_builddir} $(libppm_a_SOURCES) | \
	sed "s,config.h,../config.h,g" >  $@


clean-local:
	-rm  *.mod