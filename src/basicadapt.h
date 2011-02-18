#define HERMES_REPORT_INFO
#define HERMES_REPORT_FILE "application.log"
#include "hermes2d.h"
#include "../basic/basic.h"

// This class adds adaptivity to the module Basic.
//
// PDE: -div(c1 \nabla u) + (c2, c3) \cdot \nabla u + c4 u = c5 
//      c1 ... equation parameter, element-wise constant
//      c2 ... equation parameter, element-wise constant
//      c3 ... equation parameter, element-wise constant
//      c4 ... equation parameter, element-wise constant
//      c2 ... equation parameter, element-wise constant
//
// Possible BC: Dirichlet with piecewise-constant values (u = const)
//              Neumann with piecewise-constant normal derivatives (du/dn = const)
//              Newton with piecewise-constant parameters (const_1 u + du/dn = const_2)
//

typedef std::pair<double, double> double_pair;

using namespace RefinementSelectors;

class HERMES_API ModuleBasicAdapt : public ModuleBasic {
public:
  ModuleBasicAdapt();

  ~ModuleBasicAdapt();

  // Set adaptivity threshold.
  void set_adaptivity_threshold(double threshold);

  // Set adaptivity strategy.
  void set_adaptivity_strategy(int strategy);

  // Set adaptivity candidate list.
  void set_cand_list(std::string cand_name);

  // Set adaptivity error weights.
  void set_adaptivity_error_weights(double w1, double w2, double w3);

  // Set mesh regularity.
  void set_mesh_regularity(int mesh_regularity);

  // Set stopping relative error for adaptivity.
  void set_err_stop(double err_stop);

  // Set convergence exponent for adaptivity.
  void set_conv_exp(double conv_exp);

  // Set maximum number of degrees of freedom for adaptivity.
  void set_ndof_stop(int ndof_stop);

  // Set maximum number of adaptivity steps.
  void set_max_num_adapt_steps(int max_num_adapt_steps);
 
  // Sets some constants, performs uniform mesh refinement
  // and calculates reference solution. This needs to get 
  // done prior to adaptivity.
  bool prepare_for_adaptivity();

  // Get solution (on reference mesh).
  void get_ref_solution(Solution* s);

  // Get number of adaptivity steps done.
  int get_num_adapt_steps_done() {return this->num_adapt_steps_done;};

  // Get adaptivity time for the last step.
  double get_adaptivity_time_last();

  // Get adaptivity time (total).
  double get_adaptivity_time_total();

  // Get matrix solver time (total).
  double get_solver_time_total();

  // Get assembly time (total).
  double get_assembly_time_total();

  // Get reference space.
  void get_ref_space(H1Space* rs);

  // Get number of DOF on coarse mesh.
  int get_ndof_coarse();

  // Get number of DOF on reference mesh.
  int get_ndof_fine();

  // Get reference mesh.
  Mesh* get_ref_mesh();

  // Assumes that reference solution was computed.
  // One adaptivity step means: 
  // (1) Project ref. solution to coarse mesh.
  // (2) Calculate error estimates. 
  //     If (5a) is used, remove old reference mesh and space. 
  // (3) Refine coarse mesh.
  // (4) Create new reference mesh and space. 
  // (5) Solve on new reference mesh. There are two ways:
  //     (a) Assemble the stiffness matrix and right-hand side, and solve as usual.
  //     (b) TODO: Solve iteratively via JFNK, using last ref. solution as starting 
  //         guess. This only means to assemble the right-hand side a couple of times, 
  //         not the stiffness matrix, so it should be pretty fast.
  // (6) If (5b) is used, remove old reference mesh and space. 
  bool adaptivity_step(bool &finished, double &rel_est_rel);

  // Prepares for adaptivity and performs severla adaptivity steps.
  bool adapt(double &rel_est_rel);


protected:

  // Reference space.
  H1Space* space_ref;

  // Reference solution.
  Solution* sln_ref;

  // Reference selector.
  H1ProjBasedSelector* ref_selector;

  // Adaptivity threshold.
  double threshold;

  // Adaptivity strategy.
  int strategy;

  // Adaptivity candidate list.
  CandList cand_list;

  // Adaptivity error weights.
  double adaptivity_weight_1, adaptivity_weight_2, adaptivity_weight_3;

  // Mesh regularity.
  int mesh_regularity;

  // Stopping relative error for adaptivity (in %).
  double err_stop;

  // Convergence exponent for adaptivity.
  double conv_exp;

  // Maximum number of degrees of freedom for adaptivity.
  int ndof_stop;

  // Maximum number of adaptivity steps.
  int max_num_adapt_steps;

  // Number of adaptivity steps already done.
  int num_adapt_steps_done;

  // Time measurement of adaptivity.
  double adaptivity_time_last;
  double adaptivity_time_total;

  // Time measurement of matrix solver.
  double solver_time_total;

  // Timemeasurement of matrix assembly.
  double assembly_time_total;

  // NDOF on coarse and fine meshes.
  int ndof_coarse, ndof_fine;
};

