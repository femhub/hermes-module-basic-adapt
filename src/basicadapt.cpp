#include "basicadapt.h"

// Constructor.
ModuleBasicAdapt::ModuleBasicAdapt() : ModuleBasic()
{ 
  this->space_ref = NULL;
  this->sln_ref = new Solution();
  this->ref_selector = NULL;
  this->num_adapt_steps_done = 0;
  this->adaptivity_time_last = 0;
  this->adaptivity_time_total = 0;
  this->assembly_time_total = 0;
  this->solver_time_total = 0;
  this->ndof_coarse = this->ndof_fine = 0;
}

// Destructor.
ModuleBasicAdapt::~ModuleBasicAdapt()
{
  if (sln_ref != NULL) delete sln_ref;
}

// Set adaptivity threshold.
void ModuleBasicAdapt::set_adaptivity_threshold(double threshold) 
{
  this->threshold = threshold;
}

// Set adaptivity strategy.
void ModuleBasicAdapt::set_adaptivity_strategy(int strategy) 
{
  this->strategy = strategy;
}

// Set adaptivity candidate list.
void ModuleBasicAdapt::set_cand_list(std::string cand_name)
{
  bool found = false;
  if (cand_name == "h_iso"     ) {this->cand_list = H2D_H_ISO;      found = true;}
  if (cand_name == "h_aniso"   ) {this->cand_list = H2D_H_ANISO;    found = true;}
  if (cand_name == "p_iso"     ) {this->cand_list = H2D_P_ISO;      found = true;}
  if (cand_name == "p_aniso"   ) {this->cand_list = H2D_P_ANISO;    found = true;}
  if (cand_name == "hp_iso"    ) {this->cand_list = H2D_HP_ISO;     found = true;}
  if (cand_name == "hp_aniso_h") {this->cand_list = H2D_HP_ANISO_H; found = true;}
  if (cand_name == "hp_aniso_p") {this->cand_list = H2D_HP_ANISO_P; found = true;}
  if (cand_name == "hp_aniso"  ) {this->cand_list = H2D_HP_ANISO;   found = true;}
  if (!found) {
    warn("Possible adaptivity candidate lists: h_iso, h_aniso, p_iso, p_aniso, hp_iso, hp_aniso_h, hp_aniso_p, hp_aniso.");
    error("Unknown adaptivity candidate list %s.", cand_name.c_str());
  }
}

// Set adaptivity error weights.
void ModuleBasicAdapt::set_adaptivity_error_weights(double w1, double w2, double w3) 
{
  this->adaptivity_weight_1 = w1;
  this->adaptivity_weight_2 = w2;
  this->adaptivity_weight_3 = w3;
}

// Set mesh regularity.
void ModuleBasicAdapt::set_mesh_regularity(int mesh_regularity)
{
  this->mesh_regularity = mesh_regularity;
}

// Set stopping relative error for adaptivity.
void ModuleBasicAdapt::set_err_stop(double err_stop)
{
  this->err_stop = err_stop;
}

// Set convergence exponent for adaptivity.
void ModuleBasicAdapt::set_conv_exp(double conv_exp) 
{
  this->conv_exp = conv_exp;  
}

// Set maximum number of degrees of freedom for adaptivity.
void ModuleBasicAdapt::set_ndof_stop(int ndof_stop) 
{
  this->ndof_stop = ndof_stop;
}

// Set maximum number of adaptivity steps.
void ModuleBasicAdapt::set_max_num_adapt_steps(int max_num_adapt_steps) 
{
  this->max_num_adapt_steps = max_num_adapt_steps;
}

double ModuleBasicAdapt::get_adaptivity_time_last()
{
  return this->adaptivity_time_last;
}

double ModuleBasicAdapt::get_adaptivity_time_total()
{
  return this->adaptivity_time_total;
}

double ModuleBasicAdapt::get_assembly_time_total()
{
  return this->assembly_time_total;
}

double ModuleBasicAdapt::get_solver_time_total()
{
  return this->solver_time_total;
}

// Get reference space.
void ModuleBasicAdapt::get_ref_space(H1Space* rs) 
{
  rs->dup(this->space_ref->get_mesh());
  rs->copy_orders(this->space_ref);
}

// Get number of DOF on coarse mesh.
int ModuleBasicAdapt::get_ndof_coarse() {
  return this->ndof_coarse;
}

// Get number of DOF on fine mesh.
int ModuleBasicAdapt::get_ndof_fine() {
  return this->ndof_fine;
}

// Get reference mesh.
Mesh* ModuleBasicAdapt::get_ref_mesh() 
{
  return this->space_ref->get_mesh();
}

// Sets some constants, performs uniform mesh refinement
// and calculates reference solution. This needs to get 
// done prior to adaptivity.
bool ModuleBasicAdapt::prepare_for_adaptivity() 
{
  // Perform basic sanity checks, create mesh, perform 
  // uniform refinements, create space, register weak forms.
  bool mesh_ok = this->create_space_and_forms();
  if (!mesh_ok) return false;
  this->ndof_coarse = Space::get_num_dofs(this->space); 
  if (this->ndof_coarse <= 0) return false;

  // Initialize refinement selector.
  this->ref_selector = new H1ProjBasedSelector(this->cand_list, this->conv_exp, H2DRS_DEFAULT_ORDER);
  this->ref_selector->set_error_weights(this->adaptivity_weight_1, this->adaptivity_weight_2, 
                                        this->adaptivity_weight_3);

  // Construct globally refined reference mesh and setup reference space.
  this->space_ref = (H1Space*)construct_refined_space(this->space);
  this->ndof_fine = Space::get_num_dofs(this->space_ref);

  // Initialize the FE problem on reference mesh.
  bool is_linear = true;
  DiscreteProblem dp(this->wf, this->space_ref, is_linear);

  // Set up the solver, matrix, and rhs according to the solver selection.
  info("Initializing matrix solver, matrix, and rhs vector.");
  SparseMatrix* matrix = create_matrix(this->matrix_solver);
  Vector* rhs = create_vector(this->matrix_solver);
  Solver* solver = create_linear_solver(this->matrix_solver, matrix, rhs);

  // Begin assembly time measurement.
  TimePeriod cpu_time_assembly;
  cpu_time_assembly.tick();

  // Assemble the stiffness matrix and right-hand side vector.
  info("Assembling matrix and vector on reference mesh.");
  dp.assemble(matrix, rhs);

  // End assembly time measurement.
  this->assembly_time = cpu_time_assembly.accumulated();
  this->assembly_time_total += this->assembly_time;

  // Begin solver time measurement.
  TimePeriod cpu_time_solver;
  cpu_time_solver.tick();

  // Solve the linear system and if successful, obtain the solution.
  info("Solving on reference mesh.");
  if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), this->space_ref, this->sln_ref);
  else {
    info("Matrix solver failed.\n");
    return false;
  }

  // End solver time measurement.
  cpu_time_solver.tick();
  this->solver_time = cpu_time_solver.accumulated();
  this->solver_time_total += this->solver_time;

  // Clean up.
  info("Deleting matrix solver, matrix and rhs vector.");
  delete solver;
  delete matrix;
  delete rhs;

  return true;
}

void ModuleBasicAdapt::get_ref_solution(Solution* s) 
{
  s->copy(this->sln_ref);
}


// Perform one adaptivity step (assumes that reference solution was computed).
bool ModuleBasicAdapt::adaptivity_step(bool &finished, double &err_est_rel)
{ 
  // Initialize matrix solver.
  info("Initializing matrix solver, matrix, and rhs vector.");
  SparseMatrix* matrix = create_matrix(this->matrix_solver);
  Vector* rhs = create_vector(this->matrix_solver);
  Solver* solver = create_linear_solver(this->matrix_solver, matrix, rhs);

  // Begin solver time measurement.
  TimePeriod cpu_time_solver_1;
  cpu_time_solver_1.tick();

  // Project the reference mesh solution onto the coarse mesh.
  info("Projecting reference solution on coarse mesh (ndof: %d).", this->ndof_coarse);
  OGProjection::project_global(this->space, this->sln_ref, 
                               this->sln, this->matrix_solver); 

  // End solver time measurement.
  cpu_time_solver_1.tick();
  this->solver_time = cpu_time_solver_1.accumulated();
  this->solver_time_total += this->solver_time;

  // Clean up.
  info("Deleting matrix solver, matrix and rhs vector.");
  delete rhs;
  delete solver;
  delete matrix;

  // Begin adaptivity time measurement.
  TimePeriod cpu_time_adaptivity;
  cpu_time_adaptivity.tick();

  // Calculate element errors and total error estimate.
  info("Calculating error estimate."); 
  Adapt* adaptivity = new Adapt(this->space, HERMES_H1_NORM);
  bool solutions_for_adapt = true;
  err_est_rel = adaptivity->calc_err_est(this->sln, this->sln_ref, solutions_for_adapt, 
                                         HERMES_TOTAL_ERROR_REL | HERMES_ELEMENT_ERROR_REL) * 100;

  // Report results.
  info("ndof_coarse: %d, ndof_fine: %d, err_est_rel: %g%%", 
       Space::get_num_dofs(this->space), Space::get_num_dofs(this->space_ref), err_est_rel);

  // If err_est too large, adapt the mesh.
  if (err_est_rel < this->err_stop) {
    info("Error estimate ($g) \% below tolerance (%g) \%, not adapting.", err_est_rel, this->err_stop);
    finished = true;
    return true;
  }
  else {
    info("Adapting coarse mesh (step %d).", this->num_adapt_steps_done + 1);
    finished = adaptivity->adapt(this->ref_selector, this->threshold, this->strategy, this->mesh_regularity);
    this->ndof_coarse = Space::get_num_dofs(this->space); 
  }

  // Clean up.
  delete adaptivity;

  // End adaptivity time measurement.
  cpu_time_adaptivity.tick();
  this->adaptivity_time_last = cpu_time_adaptivity.accumulated();
  this->adaptivity_time_total += this->adaptivity_time_last;

  // Delete last reference mesh and space.
  info("Deleting reference mesh and space.");
  delete this->space_ref->get_mesh();
  delete this->space_ref;

  // Begin assembly time measurement.
  TimePeriod cpu_time_assembly;
  cpu_time_assembly.tick();

  // Construct new globally refined reference mesh and setup reference space.
  info("Constructing new reference mesh and space.");
  this->space_ref = (H1Space*)construct_refined_space(this->space);

  // Checking whether the number of degrees of freedom did not exceed the limit.
  this->ndof_fine = Space::get_num_dofs(this->space_ref); 
  if (this->ndof_fine >= this->ndof_stop) finished = true;

  // Initialize the FE problem on reference mesh.
  bool is_linear = true;
  DiscreteProblem dp(this->wf, this->space_ref, is_linear);

  // Set up the solver, matrix, and rhs according to the solver selection.
  info("Initializing matrix solver, matrix, and rhs vector.");
  matrix = create_matrix(this->matrix_solver);
  rhs = create_vector(this->matrix_solver);
  solver = create_linear_solver(this->matrix_solver, matrix, rhs);

  // Assemble the stiffness matrix and right-hand side vector.
  info("Assembling matrix and vector on reference mesh.");
  dp.assemble(matrix, rhs);

  // End assembly time measurement.
  cpu_time_assembly.tick();
  this->assembly_time = cpu_time_assembly.accumulated();
  this->assembly_time_total += this->assembly_time;

  // Begin solver time measurement.
  TimePeriod cpu_time_solver_2;
  cpu_time_solver_2.tick();

  // Solve the linear system and if successful, obtain the solution.
  info("Solving on reference mesh (ndof: %d).", this->ndof_fine);
  if(solver->solve()) Solution::vector_to_solution(solver->get_solution(), this->space_ref, this->sln_ref);
  else {
    info("Matrix solver failed.");
    return false;
  }

  // End solver time measurement.
  cpu_time_solver_2.tick();
  this->solver_time += cpu_time_solver_2.accumulated();
  this->solver_time_total += cpu_time_solver_2.accumulated();;

  // Clean up.
  info("Deleting matrix solver, matrix and rhs vector.");
  delete solver;
  delete matrix;
  delete rhs;

  // Increase counter of done adaptivity steps.
  this->num_adapt_steps_done++;

  // If computationa did not fail, return true.
  return true;
}

bool ModuleBasicAdapt::adapt(double &rel_est_rel) {
  bool finished = false;
  bool success;
  int n_step = 1;
  this->prepare_for_adaptivity();
  while (!finished) {
    success = adaptivity_step(finished, rel_est_rel);
    if (!success) {
      info("Adaptivity step %d not successful, stopping.", this->num_adapt_steps_done);
      return success;
    }
    if (this->num_adapt_steps_done >= this->max_num_adapt_steps) {
      info("Maximum number of adaptivity steps %d reached, stopping.", this->max_num_adapt_steps);
      return success;
    }
  }
  return success;
}
