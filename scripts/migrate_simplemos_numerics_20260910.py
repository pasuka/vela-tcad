"""Apply the reviewed numerical migration to the current checkout, once."""
from pathlib import Path
import validate_simplemos_poisson_precision_20260910 as candidate
s=candidate.s;REPO=s.REPO
ROOT=REPO/'build-release/phumob_numerics_production_20260910'

def replace(text,old,new,count=1):
    assert text.count(old)==count,(old[:100],text.count(old),count)
    return text.replace(old,new)

def write(rel,text):
    (REPO/rel).write_text(text,encoding='utf-8')

def main():
    for rel in ('src/solver/NewtonSolver.cpp','include/vela/solver/NewtonSolver.h','src/equation/CoupledDDAssembler.cpp','include/vela/equation/CoupledDDAssembler.h'):
        assert s.a.sha(REPO/rel)==s.a.sha(ROOT/'before'/rel),rel
    norm=(REPO/'scripts/diagnostics/simplemos_stable_merit.hpp').read_text()
    norm=norm.replace('simplemos_stable_merit','vela::stable_merit')
    write('include/vela/numerics/StableMeritComparison.h',norm)
    poisson=candidate.HEADER.read_text().replace('simplemos_poisson_precision','vela::poisson_precision').replace('// Diagnostic only: preserve the frozen double geometry/material coefficients.','// Extended arithmetic for Poisson; geometry and material inputs retain their stored precision.')
    write('include/vela/equation/ExtendedPoissonResidual.h',poisson)
    rel='include/vela/equation/CoupledDDAssembler.h';text=(REPO/rel).read_text()
    marker='    VectorXd pack(const CoupledDDState& state) const;'
    text=replace(text,marker,'    /// Opt-in binary128 Poisson evaluation from packed Boltzmann coordinates.\n    void setExtendedPoissonResidual(bool enabled);\n    bool usesExtendedPoissonResidual() const { return extendedPoissonResidual_; }\n\n'+marker)
    text=replace(text,'    bool usesFermiDirac_ = false;','    bool usesFermiDirac_ = false;\n    bool extendedPoissonResidual_ = false;');write(rel,text)
    rel='src/equation/CoupledDDAssembler.cpp';text=(REPO/rel).read_text()
    text='#include "vela/equation/ExtendedPoissonResidual.h"\n'+text
    marker='VectorXd CoupledDDAssembler::residual('
    setter='''void CoupledDDAssembler::setExtendedPoissonResidual(bool enabled)
{
    if (enabled && (usesFermiDirac_ || electronQuantumPotentialConfig_.enabled ||
                    (electronQuantumPotential_V_.size() && !electronQuantumPotential_V_.isZero(0.0))))
        throw std::invalid_argument("Extended Poisson residual requires classical Boltzmann carriers without a quantum potential.");
    extendedPoissonResidual_ = enabled;
}

'''
    text=replace(text,marker,setter+marker)
    start=text.index('VectorXd CoupledDDAssembler::residualImpl(');end=text.index('std::vector<CoupledDDCarrierTermDiagnostic>',start);frag=text[start:end]
    hook=candidate.HOOK
    hook=replace(hook,'    if (const char* precisionText=std::getenv("VELA_POISSON_PRECISION")) {','    if (extendedPoissonResidual_) {')
    hook=replace(hook,'        const std::string precisionMode(precisionText);\n        if(precisionMode!="kernel" && precisionMode!="packed")throw std::runtime_error("Invalid Poisson precision mode");\n','')
    hook=hook.replace('simplemos_poisson_precision','poisson_precision').replace('precisionMode=="packed"','true').replace('Unsupported Poisson precision diagnostic branch','Extended Poisson residual does not support feedback substitutions or thermionic contacts')
    frag=replace(frag,'    for (const auto& [node, value] : bcs.psi)',hook+'\n    for (const auto& [node, value] : bcs.psi)')
    text=text[:start]+frag+text[end:];write(rel,text)
    rel='include/vela/solver/NewtonSolver.h';text=(REPO/rel).read_text()
    text=replace(text,'    int linearRefinementIterations = 0;','''    int linearRefinementIterations = 0;
    std::string poissonResidualPrecision = "double"; ///< "double" or opt-in "binary128" packed Boltzmann evaluation.
    bool stableMeritComparison = false; ///< Exact-sign squared-norm comparison for unweighted L2 merit searches.
    bool exactDirichletUpdates = false; ///< Remove direct-solve roundoff from exact Dirichlet identity-row updates.''')
    text=text.replace('configureQuasiFermiReferences','configureAssemblerNumerics');write(rel,text)
    rel='src/solver/NewtonSolver.cpp';text=(REPO/rel).read_text();text='#include "vela/numerics/StableMeritComparison.h"\n'+text
    marker='NewtonConfig newtonConfigFromJson('
    validation='''namespace {
void validateNumericalPrecisionOptions(const NewtonConfig& cfg)
{
    if (cfg.poissonResidualPrecision != "double" && cfg.poissonResidualPrecision != "binary128")
        throw std::invalid_argument("poisson_residual_precision must be double or binary128.");
    if (cfg.poissonResidualPrecision == "binary128" &&
        (cfg.carrierStatistics.model != "boltzmann" || cfg.electronQuantumPotential.enabled))
        throw std::invalid_argument("binary128 Poisson residual requires classical Boltzmann statistics without a quantum potential.");
    if (cfg.stableMeritComparison && (cfg.residualNorm != "l2" ||
        cfg.lineSearchMode != "merit" || cfg.globalContinuityClosure.mode != "off"))
        throw std::invalid_argument("stable_merit_comparison requires unweighted l2, merit line search, and global continuity merit off.");
    if (cfg.exactDirichletUpdates && cfg.carrierRegularizationScale != 0.0)
        throw std::invalid_argument("exact_dirichlet_updates requires unregularized Dirichlet identity rows.");
}
} // namespace

'''
    text=replace(text,marker,validation+marker)
    marker='    cfg.linearRefinementIterations = json.value("linear_refinement_iterations", 0);'
    text=replace(text,marker,marker+'''
    cfg.poissonResidualPrecision = json.value("poisson_residual_precision", cfg.poissonResidualPrecision);
    cfg.stableMeritComparison = json.value("stable_merit_comparison", cfg.stableMeritComparison);
    cfg.exactDirichletUpdates = json.value("exact_dirichlet_updates", cfg.exactDirichletUpdates);''')
    start=text.index('NewtonConfig newtonConfigFromJson(');end=text.index('NewtonSolver::NewtonSolver(',start);frag=text[start:end]
    frag=replace(frag,'    return cfg;','    validateNumericalPrecisionOptions(cfg);\n    return cfg;');text=text[:start]+frag+text[end:]
    marker='    if (cfg_.linearRefinementIterations < 0 || cfg_.linearRefinementIterations > 10)'
    text=replace(text,marker,'    validateNumericalPrecisionOptions(cfg_);\n'+marker)
    text=text.replace('configureQuasiFermiReferences','configureAssemblerNumerics')
    marker='void NewtonSolver::configureAssemblerNumerics(\n    CoupledDDAssembler& assembler) const\n{'
    text=replace(text,marker,marker+'\n    assembler.setExtendedPoissonResidual(cfg_.poissonResidualPrecision == "binary128");')
    for car,short in (('electron','phin'),('hole','phip')):
        old=f'''static_cast<long double>(
                        initial.{car}QuasiFermiReferenceAt(i)) +
                    static_cast<long double>(initial.{short}Increment(i)) -
                    static_cast<long double>(
                        assembler.{car}QuasiFermiReferenceAt(node))'''
        new=f'''(static_cast<long double>(
                        initial.{car}QuasiFermiReferenceAt(i)) -
                     static_cast<long double>(
                        assembler.{car}QuasiFermiReferenceAt(node))) +
                    static_cast<long double>(initial.{short}Increment(i))'''
        text=replace(text,old,new)
    marker='        const VectorXd rawStep = step;'
    text=replace(text,marker,'''
        // These rows are exact identity equations. A sparse direct solve can
        // pollute their zero update with tiny roundoff, creating artificial
        // merit decreases after all physical rows have reached their floor.
        if (cfg_.exactDirichletUpdates) {
            for (const auto& [node, value] : bcs.psi)
                step(static_cast<int>(node)) = -r(static_cast<int>(node));
            for (const auto& [node, value] : bcs.phin)
                step(N + static_cast<int>(node)) = -r(N + static_cast<int>(node));
            for (const auto& [node, value] : bcs.phip)
                step(2 * N + static_cast<int>(node)) = -r(2 * N + static_cast<int>(node));
        }
'''+marker)
    marker='        const auto runLineSearch = [&](const VectorXd& trialStep) {'
    text=replace(text,marker,'''        const auto stableDecrease = [&](const VectorXd& trialResidual, Real) {
            return stable_merit::compare(r, trialResidual).accepted;
        };
'''+marker)
    start=text.index(marker);end=text.index('        auto ls = runLineSearch(step);',start);frag=text[start:end]
    frag=replace(frag,'                    : BacktrackingLineSearch::DecreaseAcceptFunction{});','                    : (cfg_.stableMeritComparison\n                        ? BacktrackingLineSearch::DecreaseAcceptFunction(stableDecrease)\n                        : BacktrackingLineSearch::DecreaseAcceptFunction{}));')
    text=text[:start]+frag+text[end:];write(rel,text)
    print('Migrated explicit numerical options and cancellation-free warm restarts.')

if __name__=='__main__':main()
