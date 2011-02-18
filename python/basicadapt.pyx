from numpy cimport ndarray
from libcpp cimport bool
from libcpp.vector cimport vector
from libcpp.pair cimport pair

cimport basicadapt_defs
from hermes2d.hermes2d cimport Solution
#from hermes2d.hermes2d cimport Space

cdef vector[int] array2vector_int(a):
    cdef vector[int] v
    for i in range(len(a)):
        v.push_back(a[i])
    return v

cdef vector[double] array2vector_double(a):
    cdef vector[double] v
    for i in range(len(a)):
        v.push_back(a[i])
    return v

cdef vector[pair[double, double]] array2vector_double_pair(a):
    cdef vector[pair[double, double]] v
    for i in range(len(a)):
        p1, p2 = a[i]
        v.push_back(pair[double, double](float(p1), float(p2)))
    return v

cdef class ModuleBasicAdapt:
    cdef basicadapt_defs.ModuleBasicAdapt *thisptr

    def __init__(self):
        self.thisptr = new basicadapt_defs.ModuleBasicAdapt()

    def __dealloc__(self):
        del self.thisptr

    # ---------------------------------------
    # THIS PART IS THE SAME AS IN ModuleBasic
    # ---------------------------------------

    def set_mesh_str(self, mesh_str):
        return self.thisptr.set_mesh_str(mesh_str)

    def set_initial_mesh_refinement(self, int init_ref_num):
        self.thisptr.set_initial_mesh_refinement(init_ref_num)

    def set_initial_poly_degree(self, int p):
        self.thisptr.set_initial_poly_degree(p)

    def set_matrix_solver(self, solver_name):
        self.thisptr.set_matrix_solver(solver_name)

    def set_material_markers(self, mat_markers):
        self.thisptr.set_material_markers(array2vector_int(mat_markers))

    def set_c1_array(self, c1_array):
        self.thisptr.set_c1_array(array2vector_double(c1_array))

    def set_c2_array(self, c2_array):
        self.thisptr.set_c2_array(array2vector_double(c2_array))

    def set_c3_array(self, c3_array):
        self.thisptr.set_c3_array(array2vector_double(c3_array))

    def set_c4_array(self, c4_array):
        self.thisptr.set_c4_array(array2vector_double(c4_array))

    def set_c5_array(self, c5_array):
        self.thisptr.set_c5_array(array2vector_double(c5_array))

    def set_dirichlet_markers(self, bdy_markers_dirichlet):
        self.thisptr.set_dirichlet_markers(array2vector_int(bdy_markers_dirichlet))

    def set_dirichlet_values(self,  bdy_markers_dirichlet, bdy_values_dirichlet):
        self.thisptr.set_dirichlet_values(array2vector_int(bdy_markers_dirichlet), array2vector_double(bdy_values_dirichlet))

    def set_neumann_markers(self, bdy_markers_neumann):
        self.thisptr.set_neumann_markers(array2vector_int(bdy_markers_neumann))

    def set_neumann_values(self, bdy_values_neumann):
        self.thisptr.set_neumann_values(array2vector_double(bdy_values_neumann))

    def set_newton_markers(self, bdy_markers_newton):
        self.thisptr.set_newton_markers(array2vector_int(bdy_markers_newton))

    def set_newton_values(self, bdy_values_newton):
        vvv = array2vector_double_pair(bdy_values_newton)
        self.thisptr.set_newton_values(vvv)

    def get_assembly_time(self):
        return self.thisptr.get_assembly_time()

    def get_assembly_time_total(self):
        return self.thisptr.get_assembly_time_total()

    def get_solver_time(self):
        return self.thisptr.get_solver_time()

    def get_solver_time_total(self):
        return self.thisptr.get_solver_time_total()

    def get_solution(self):
        s = Solution()
        self.thisptr.get_solution(s.getptr())
        return s

    #def get_space(self):
    #    s = Space()
    #    self.thisptr.get_space(s.getptr())
    #    return s

    def calculate(self):
        success = self.thisptr.calculate()
        return success

    # -----------------------------------------
    # THIS PART IS NOT CONTAINED IN ModuleBasic    
    # -----------------------------------------

    def set_adaptivity_threshold(self, double threshold):
        self.thisptr.set_adaptivity_threshold(threshold)  

    def set_adaptivity_strategy(self, int strategy):
        self.thisptr.set_adaptivity_strategy(strategy)  

    def set_cand_list(self, cand_name):
        self.thisptr.set_cand_list(cand_name)

    def set_adaptivity_error_weights(self, double w1, double w2, double w3):
        self.thisptr.set_adaptivity_error_weights(w1, w2, w3)  

    def set_mesh_regularity(self, int mesh_regularity):
        self.thisptr.set_mesh_regularity(mesh_regularity)  

    def set_err_stop(self, double err_stop):
        self.thisptr.set_err_stop(err_stop)  

    def set_conv_exp(self, double conv_exp):
        self.thisptr.set_conv_exp(conv_exp)  

    def set_ndof_stop(self, int ndof_stop):
        self.thisptr.set_ndof_stop(ndof_stop)  

    def set_max_num_adapt_steps(self, int max_num_adapt_steps):
        self.thisptr.set_max_num_adapt_steps(max_num_adapt_steps)  

    def prepare_for_adaptivity(self):
        cdef bool success
        success = self.thisptr.prepare_for_adaptivity()
        return success

    def get_ref_solution(self):
        s = Solution()
        self.thisptr.get_ref_solution(s.getptr())
        return s

    def get_num_adapt_steps_done(self):
        return self.thisptr.get_num_adapt_steps_done()

    def get_adaptivity_time_last(self):
        return self.thisptr.get_adaptivity_time_last()

    def get_adaptivity_time_total(self):
        return self.thisptr.get_adaptivity_time_total()

    def adaptivity_step(self):
        cdef bool _finished = False
        cdef double _err_est_rel
        success = self.thisptr.adaptivity_step(_finished, _err_est_rel)
        return (success, _finished, _err_est_rel)

    def adapt(self):
        cdef double _err_est_rel
        success = self.thisptr.adapt(_err_est_rel)
        return (success, _err_est_rel)

    def get_ndof_coarse(self):
        return self.thisptr.get_ndof_coarse()

    def get_ndof_fine(self):
        return self.thisptr.get_ndof_fine()

    def get_num_base_elements(self):
        return self.thisptr.get_num_base_elements()










