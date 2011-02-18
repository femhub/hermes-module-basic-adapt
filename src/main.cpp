// WARNING: Hermes needs to be compiled with Glut 
// for this module to work. Check the CMake.vars file
// in ~/repos/hermes/hermes2d/. 

// How to run it:   module-basicadapt model.cfg

//#define HERMES_REPORT_INFO
#define HERMES_REPORT_ALL
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "disc.h"
#include "basicadapt.h"

using namespace RefinementSelectors;

int main(int argc, char* argv[])
{
  // Check the number of command line arguments.
  if(argc != 2) error("Configuration file missing.");

  // Open configuration file.
  FILE* f = fopen(argv[1], "r");
  if(f == NULL) error("Cannot open file %s.", argv[1]);

  // Initialize the ModuleBasicAdapt class.
  ModuleBasicAdapt BA;

  /*** READ MESH (TEMPORARILY AS A STRING )***/

  // Convert mesh string into mesh, and check its correctness.
  bool mesh_ok = BA.set_mesh_str("\na = 1.0  # size of the mesh\nb = sqrt(2)/2\n\nvertices =\n{\n  { 0, -a },    # vertex 0\n  { a, -a },    # vertex 1\n  { -a, 0 },    # vertex 2\n  { 0, 0 },     # vertex 3\n  { a, 0 },     # vertex 4\n  { -a, a },    # vertex 5\n  { 0, a },     # vertex 6\n  { a*b, a*b }  # vertex 7\n}\n\nelements =\n{\n  { 0, 1, 4, 3, 0 },  # quad 0\n  { 3, 4, 7, 0 },     # tri 1\n  { 3, 7, 6, 0 },     # tri 2\n  { 2, 3, 6, 5, 0 }   # quad 3\n}\n\nboundaries =\n{\n  { 0, 1, 1 },\n  { 1, 4, 2 },\n  { 3, 0, 4 },\n  { 4, 7, 2 },\n  { 7, 6, 2 },\n  { 2, 3, 4 },\n  { 6, 5, 2 },\n  { 5, 2, 3 }\n}\n\ncurves =\n{\n  { 4, 7, 45 },  # +45 degree circular arcs\n  { 7, 6, 45 }\n}\n");
  if(!mesh_ok) error("No elements found in mesh.");

  /*** READ DATA - GENERAL ***/

  // Read number of initial uniform mesh refinements.
  int init_ref_num;
  if(!Get(f, &init_ref_num)) error("Could not read number of initial mesh refinements.");
  info("init_ref_num: %d", init_ref_num);
  BA.set_initial_mesh_refinement(init_ref_num);

  // Read initial polynomial degree of elements.
  int init_p;
  if(!Get(f, &init_p)) error("Could not read initial polynomial degree.");
  info("init_p: %d", init_p);
  BA.set_initial_poly_degree(init_p);

  // Read matrix solver.
  char* str = new char[255];
  if(!Get(f, str)) error("Could not read matrix solver name.");
  info("matrix solver: %s", str);
  BA.set_matrix_solver(std::string(str));
  delete str;

  /*** READ DATA - MODEL-SPECIFIC ***/

#include "read_model_data.cpp"

  /*** READ DATA - SPECIFIC FOR ADAPTIVITY ***/

#include "read_adaptivity_data.cpp"


  /*** SOLVE THE PROBLEM ***/

  // Solve the problem using automatic adaptivity.
  double err_rel_est;
  bool success = BA.adapt(err_rel_est);
  info("=========================");
  info("Adaptivity step %d finished.", BA.get_num_adapt_steps_done());
  info("ndof_coarse: %d, ndof_fine: %d", BA.get_ndof_coarse(), BA.get_ndof_fine());
  info("assembly time (last step): %g s", BA.get_assembly_time());
  info("assembly time (total): %g s", BA.get_assembly_time_total());
  info("solver time (last step): %g s", BA.get_solver_time());
  info("solver time (total): %g s", BA.get_solver_time_total());
  info("adaptivity time (last step): %g s", BA.get_adaptivity_time_last());
  info("adaptivity time (total): %g s", BA.get_adaptivity_time_total());
  if (!success) error("Computation failed.");

  // Show reference solution.
  char* title = new char[100];
  sprintf(title, "Reference solution, step %d", BA.get_num_adapt_steps_done());
  Solution sln_ref;
  BA.get_ref_solution(&sln_ref);
  ScalarView view(title, new WinGeom(0, 0, 440, 350));
  view.show_mesh(false);
  view.show(&sln_ref);

  // Show gradient of reference solution.
  sprintf(title, "Gradient, step %d", BA.get_num_adapt_steps_done());
  ScalarView gradview(title, new WinGeom(445, 0, 440, 350));
  MagFilter grad(Hermes::Tuple<MeshFunction *>(&sln_ref, &sln_ref), 
                 Hermes::Tuple<int>(H2D_FN_DX, H2D_FN_DY));
  gradview.show_mesh(false);
  gradview.show(&grad);

  // Show reference space.
  BCTypes bctypes;
  BCValues bcvalues;
  H1Space space_ref(BA.get_ref_mesh(), &bctypes, &bcvalues, 1); // FIXME: this is a hack since constructor 
                                                                  // to Space needs some mesh. The mesh is not used.
  BA.get_ref_space(&space_ref);
  sprintf(title, "Reference space, step %d", BA.get_num_adapt_steps_done());
  OrderView oview(title, new WinGeom(890, 0, 440, 350));
  oview.show(&space_ref);

  // Wait.
  View::wait();
  
  return 1;
}
