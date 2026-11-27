// Read-only output helpers shared by the existing SLOT and Newton entry points.
nlohmann::json writeWeakRowLedger(const std::filesystem::path& path,
    const NewtonProblem& problem,const vela::NewtonCarrierTermDiagnosticsEvaluation& terms) {
 using namespace vela;
 std::vector<NewtonCarrierRowConvergenceViolation> rows;
 const auto evaluation=evaluateCarrierRowConvergence(terms.rows,problem.newton.carrierRowConvergence,&rows);
 std::vector<std::set<std::string>> materials(problem.mesh.numNodes());
 for(const auto& cell:problem.mesh.cells())for(auto node:cell.node_ids)
  materials[node].insert(problem.mesh.getRegion(cell.region_id).material);
 std::ofstream f(path);f.exceptions(std::ios::badbit|std::ios::failbit);
 const auto& units=problem.newton.inputScaling.unitSystem();
 const double currentScale=terms.continuityResidualToAPerUm;
 if(!(currentScale>0.) || !std::isfinite(currentScale))
  throw std::runtime_error("Weak row diagnostics require the assembler's continuity current conversion");
 f<<std::setprecision(17)<<"node_id,carrier,materials,residual_internal,absolute_residual_A_per_um,physical_scale_internal,scale_internal,ratio,density_gate_units,flux_internal,recombination_internal,impact_internal,source_qualified,density_qualified,flux_qualified,qualified,floor_limited,blocks_acceptance\n";
 for(const auto& row:rows){std::string material;
  for(const auto& x:materials[row.nodeId]){if(!material.empty())material+='|';material+=x;}
  // CSV quoting also permits future material names containing commas/quotes.
  std::string quoted;for(char c:material){quoted+=c;if(c=='"')quoted+='"';}
  f<<row.nodeId<<','<<row.carrier<<",\""<<quoted<<"\","<<row.residual<<','<<std::abs(row.residual)*currentScale<<','
   <<row.physicalScale<<','<<row.scale<<','<<row.ratio<<','<<row.carrierDensity_m3<<','<<row.flux<<','
   <<row.recombination<<','<<row.impact<<','<<row.sourceQualified<<','<<row.densityQualified<<','<<row.fluxQualified<<','
   <<row.qualified<<','<<row.floorLimited<<','<<(evaluation.enforced&&row.blocksAcceptance)<<'\n';
 }
 auto report=carrierRowConvergenceJson(evaluation);
 report["ledger_file"]=path.string();report["current_to_A_per_um"]=currentScale;
 report["current_conversion"]="q * continuityResidualScale * physical line factor * 1e-6";
 report["scope"]="all active continuity rows if configured row gate is enabled; original production gate, unchanged";
 report["materials"]="all incident cell materials; interface membership is retained";
 report["density_gate_units_to_m3"]=problem.newton.inputScaling.isUnitScaling()?units.concentrationM3PerInternal():1.;
 double maxFlux[2]={0.,0.},maxSource[2]={0.,0.};
 for(const auto& r:terms.rows)for(int b=0;b<2;++b){
  maxFlux[b]=std::max(maxFlux[b],std::abs(b==0?r.electronFluxAbsSum:r.holeFluxAbsSum));
  maxSource[b]=std::max({maxSource[b],std::abs(b==0?r.electronRecombination:r.holeRecombination),std::abs(b==0?r.electronImpact:r.holeImpact)});
 }
 const auto& c=problem.newton.carrierRowConvergence;
 for(int b=0;b<2;++b)report["qualification_scales"][b==0?"electron":"hole"]={
  {"max_flux_internal",maxFlux[b]},{"max_source_internal",maxSource[b]},
  {"flux_threshold_internal",std::max(c.minFluxScale,c.minFluxScaleFraction*maxFlux[b])},
  {"source_global_threshold_internal",c.minSourceGlobalFraction*maxSource[b]},
  {"source_absolute_threshold_internal",c.minSourceScale},{"source_row_fraction",c.minSourceScaleFraction},
  {"density_threshold_gate_units",c.minCarrierDensity_m3}};
 return report;
}
