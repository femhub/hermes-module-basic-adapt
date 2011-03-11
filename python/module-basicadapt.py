# If you have installed the Hermes wrappers and the module wrapper into
# directories which are not on default Python path, you must point to them here:
import sys
sys.path.append("/home/milan/myhermes/lib/python")

from hermes2d.hermes2d import Linearizer
from hermes2d.plot import sln2png

# If you ran 'make install' and appended 'sys.path' appropriately (see above),
# you may import the module from anywhere by the following line:
#from hermes_modules.basicadapt import ModuleBasicAdapt

# In current module's directory, you may import the wrapper also locally:
from basicadapt import ModuleBasicAdapt

def main():
    e = ModuleBasicAdapt()

    # Set problem-dependent data.
    mesh_ok = e.set_mesh_str("\na = 1.0  # size of the mesh\nb = sqrt(2)/2\n\nvertices =\n{\n  { 0, -a },    # vertex 0\n  { a, -a },    # vertex 1\n  { -a, 0 },    # vertex 2\n  { 0, 0 },     # vertex 3\n  { a, 0 },     # vertex 4\n  { -a, a },    # vertex 5\n  { 0, a },     # vertex 6\n  { a*b, a*b }  # vertex 7\n}\n\nelements =\n{\n  { 0, 1, 4, 3, 0 },  # quad 0\n  { 3, 4, 7, 0 },     # tri 1\n  { 3, 7, 6, 0 },     # tri 2\n  { 2, 3, 6, 5, 0 }   # quad 3\n}\n\nboundaries =\n{\n  { 0, 1, 1 },\n  { 1, 4, 2 },\n  { 3, 0, 4 },\n  { 4, 7, 2 },\n  { 7, 6, 2 },\n  { 2, 3, 4 },\n  { 6, 5, 2 },\n  { 5, 2, 3 }\n}\n\ncurves =\n{\n  { 4, 7, 45 },  # +45 degree circular arcs\n  { 7, 6, 45 }\n}\n");
    assert mesh_ok is True
    e.set_initial_mesh_refinement(2)
    e.set_initial_poly_degree(4)
    solver_name = "umfpack"   # choose from amesos, aztecoo, mumps, petsc, superlu, umfpack
    e.set_matrix_solver(solver_name)
    print "Matrix solver:", solver_name
    e.set_material_markers([0])
    e.set_c1_array([1])
    e.set_c2_array([0])
    e.set_c3_array([0])
    e.set_c4_array([0])
    e.set_c5_array([1])
    e.set_dirichlet_markers([4])
    e.set_dirichlet_values([4], [0])
    e.set_neumann_markers([1, 3])
    e.set_neumann_values([0, 0])
    e.set_newton_markers([2])
    e.set_newton_values( [(1, 1)] )

    # Set adaptivity data.
    e.set_adaptivity_threshold(0.3)
    e.set_adaptivity_strategy(0)
    e.set_cand_list("hp_aniso_h")
    e.set_adaptivity_error_weights(2.0, 1.0, 1.4142136)
    e.set_mesh_regularity(-1)
    e.set_err_stop(0.1)
    e.set_conv_exp(1.0)
    e.set_ndof_stop(60000)
    e.set_max_num_adapt_steps(5)

    # Solve the problem adaptively.
    success, err_rel_est = e.adapt()
    assert success is True
    sln = e.get_ref_solution()
    print "Final relative error:", err_rel_est
    print "Final ndof_coarse, ndof_fine:", (e.get_ndof_coarse(), e.get_ndof_fine())
    print "Assembly time (last step):", e.get_assembly_time()
    print "Assembly time (total):", e.get_assembly_time_total()
    print "Solver time (last step):", e.get_solver_time()
    print "Solver time (total):", e.get_solver_time_total()
    print "Adaptivity time (last step):", e.get_adaptivity_time_last()
    print "Adaptivity time (total):", e.get_adaptivity_time_total()
    print "Saving reference solution to 'solution.png'"
    sln2png(sln, "ref_solution.png")
    raw_input()


if __name__ == "__main__":
    main()
