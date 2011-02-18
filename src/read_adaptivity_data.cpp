  /*** READ DATA - SPECIFIC FOR ADAPTIVITY ***/

  // Read adaptivity threshold.
  double threshold;
  if(!Get(f, &threshold)) error("Could not read adaptivity threshold.");
  info("threshold: %g", threshold);
  BA.set_adaptivity_threshold(threshold);

  // Read adaptivity strategy.
  int strategy;
  if(!Get(f, &strategy)) error("Could not read adaptivity strategy.");
  info("strategy: %d", strategy);
  BA.set_adaptivity_strategy(strategy);
  
  // Read candidate list.
  str = new char[255];
  if(!Get(f, str)) error("Could not read adaptivity candidate list.");
  info("cand_list: %s", str);
  BA.set_cand_list(std::string(str));
  delete str;

  // Read error weights for adaptivity.
  double w1, w2, w3;
  if(!Get(f, &w1)) error("Could not read an adaptivity error weight.");
  if(!Get(f, &w2)) error("Could not read an adaptivity error weight.");
  if(!Get(f, &w3)) error("Could not read an adaptivity error weight.");
  info("error weights: %g %g %g", w1, w2, w3);
  BA.set_adaptivity_error_weights(w1, w2, w3);

  // Read mesh regularity.
  int mesh_regularity;
  if(!Get(f, &mesh_regularity)) error("Could not read mesh regularity.");
  info("mesh regularity: %d", mesh_regularity);
  BA.set_mesh_regularity(mesh_regularity);

  // Read stopping relative error for adaptivity.
  double err_stop;
  if(!Get(f, &err_stop)) error("Could not read stopping error value for adaptivity.");
  info("err_stop: %g", err_stop);
  BA.set_err_stop(err_stop);

  // Read convergence exponent for adaptivity.
  double conv_exp;
  if(!Get(f, &conv_exp)) error("Could not read convergence exponent for adaptivity.");
  info("conv_exp: %g", conv_exp);
  BA.set_conv_exp(conv_exp);

  // Read maximum number of degrees of freedom for adaptivity.
  int ndof_stop;
  if(!Get(f, &ndof_stop)) error("Could not read maximum number of degrees of freedom for adaptivity.");
  info("ndof_stop: %d", ndof_stop);
  BA.set_ndof_stop(ndof_stop);

  // Read maximum number of adaptivity steps.
  int max_num_adapt_steps;
  if(!Get(f, &max_num_adapt_steps)) error("Could not read maximum number of adaptivity steps.");
  info("max_num_adapt_steps: %d", max_num_adapt_steps);
  BA.set_max_num_adapt_steps(max_num_adapt_steps);
