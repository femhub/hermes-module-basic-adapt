from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.pair cimport pair

from hermes2d.hermes2d_defs cimport Solution

# FIXME; Space needs wrapping
# from hermes2d.hermes2d_defs cimport Space

cdef extern from "basicadapt.h":

    cdef cppclass ModuleBasicAdapt:

        # This part is identical to ModuleBasic:
        bool set_mesh_str(char *mesh_str)
        void set_initial_mesh_refinement(int init_ref_num)
        void set_matrix_solver(char *solver_name)
        void set_initial_poly_degree(int p)
        void set_material_markers(vector[int] &mat_markers)
        void set_c1_array(vector[double] &c1_array)
        void set_c2_array(vector[double] &c2_array)
        void set_c3_array(vector[double] &c3_array)
        void set_c4_array(vector[double] &c4_array)
        void set_c5_array(vector[double] &c5_array)
        void set_dirichlet_markers(vector[int] &bdy_markers_dirichlet)
        void set_dirichlet_values(vector[int] &bdy_markers_dirichlet, vector[double] &bdy_values_dirichlet)
        void set_neumann_markers(vector[int] &bdy_markers_neumann)
        void set_neumann_values(vector[double] &bdy_values_neumann)
        void set_newton_markers(vector[int] &bdy_markers_newton)
        void set_newton_values(vector[pair[double, double]] &bdy_values_newton)
        double get_assembly_time()
        double get_solver_time()
        void get_solution(Solution* s)
        # FIXME: Space needs wrapping.
        # void get_space(Space* s)
        bool calculate()
        int get_num_base_elements()

        # New things which are not in ModuleBasic:
        void set_adaptivity_threshold(double threshold)  
        void set_adaptivity_strategy(int strategy)
        void set_cand_list(char *cand_name)
        void set_adaptivity_error_weights(double w1, double w2, double w3)
        void set_mesh_regularity(int mesh_regularity)
        void set_err_stop(double err_stop)
        void set_conv_exp(double conv_exp)
        void set_ndof_stop(int ndof_stop)
        void set_max_num_adapt_steps(int max_num_adapt_steps)
        bool prepare_for_adaptivity()
        void get_ref_solution(Solution* s)
        int get_num_adapt_steps_done()
        double get_assembly_time_total()
        double get_solver_time_total()
        double get_adaptivity_time_last()
        double get_adaptivity_time_total()
        # FIXME: Space needs wrapping.
        # void get_ref_space(H1Space* rs)
        bool adaptivity_step(bool &finished, double &err_est_rel)
        bool adapt(double &err_est_rel)
        int get_ndof_coarse()
        int get_ndof_fine()






























