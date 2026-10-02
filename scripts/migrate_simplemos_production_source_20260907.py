"""One-shot, checked source migration; preimages are preserved separately."""
from pathlib import Path
R=Path(__file__).resolve().parents[1]
def edit(name,old,new):
    p=R/name;s=p.read_text();assert s.count(old)==1,(name,old[:100],s.count(old));p.write_text(s.replace(old,new),newline='\n')

edit('include/vela/equation/DDAssembler.h','    bool poissonEdgeCoupling = false;', '''    // Charge-only policy: continuity/recombination volumes are independent.
    std::string poissonChargeNodeVolume = "inherit";
    // Explicit mesh-edge multipliers. Empty retains the historical geometry.
    std::vector<Real> transportEdgeCouplingRatios;
    bool poissonEdgeCoupling = false;''')
edit('include/vela/equation/DDAssembler.h','        return poissonEdgeCoupling || transportEdgeCoupling ||','''        return poissonChargeNodeVolume != "inherit" ||
               !transportEdgeCouplingRatios.empty() ||
               poissonEdgeCoupling || transportEdgeCoupling ||''')
edit('include/vela/equation/AssemblerUtils.h','#include "vela/equation/ChargeSpec.h"','#include "vela/equation/ChargeSpec.h"\n#include "vela/equation/DDAssembler.h"')
name='include/vela/equation/AssemblerUtils.h'
s=(R/name).read_text();start=s.index('inline std::vector<Real> computeTransportNodeVolumes')
helper='''// One effective geometry for assembly, terminal integration and diagnostics.
inline std::vector<Real> computeEffectiveTransportEdgeCouplings(
    const DeviceMesh& mesh, const std::vector<std::vector<Index>>& edgeCells,
    const std::vector<Material>& cellMaterials,
    const RegionResolvedInterfaceAssemblyConfig& config)
{
    const auto& ratios = config.transportEdgeCouplingRatios;
    if (!ratios.empty() && (config.transportEdgeCoupling ||
                           ratios.size() != mesh.edges().size()))
        throw std::invalid_argument("Transport edge ratios require one value per mesh edge and cannot combine with transport_edge_coupling.");
    auto result = config.transportEdgeCoupling
        ? computeTransportEdgeCouplings(mesh, edgeCells, cellMaterials)
        : computeEdgeCouplings(mesh);
    for (std::size_t e = 0; e < ratios.size(); ++e) {
        if (!std::isfinite(ratios[e]) || ratios[e] < 0.0)
            throw std::invalid_argument("Transport edge ratios must be finite and nonnegative.");
        result[e] *= ratios[e];
    }
    return result;
}

'''
(R/name).write_text(s[:start]+helper+s[start:],newline='\n')
name='src/equation/CoupledDDAssembler.cpp'
for car in ('Electron','Hole','Dopant'):
    old=f'''    , poisson{car}Vol_(
          scaling.regionResolvedInterfaceAssembly'''
    new=f'''    , poisson{car}Vol_(
          scaling.regionResolvedInterfaceAssembly.poissonChargeNodeVolume == "signed_transport"
              ? detail::computeTransportSignedAverageBoxNodeVolumes(mesh, cellMaterials_)
              : scaling.regionResolvedInterfaceAssembly'''
    edit(name,old,new)
edit(name,'''    , couple_(scaling.regionResolvedInterfaceAssembly.transportEdgeCoupling
          ? detail::computeTransportEdgeCouplings(
                mesh, edgeCells_, cellMaterials_)
          : detail::computeEdgeCouplings(mesh))''','''    , couple_(detail::computeEffectiveTransportEdgeCouplings(
          mesh, edgeCells_, cellMaterials_, scaling.regionResolvedInterfaceAssembly))''')
edit('src/post/ContactCurrent.cpp','    , couple_(detail::computeEdgeCouplings(mesh))','''    , couple_(detail::computeEffectiveTransportEdgeCouplings(
          mesh, edgeCells_, detail::buildCellMaterials(mesh, matdb, temperature_K),
          scaling.regionResolvedInterfaceAssembly))''')
edit('src/tools/vela_example_runner.cpp','''    vela::MobilityModelConfig mobilityConfig = newton.mobility;
    vela::detail::updateSurfaceMobilityCellGeometry(''','''    const auto effectiveCoupling = vela::detail::computeEffectiveTransportEdgeCouplings(
        mesh, edgeCells, cellMaterials, newton.regionResolvedInterfaceAssembly);
    vela::MobilityModelConfig mobilityConfig = newton.mobility;
    vela::detail::updateSurfaceMobilityCellGeometry(''')
edit('src/tools/vela_example_runner.cpp','units.internalLengthToMeters(edge.couple)','units.internalLengthToMeters(effectiveCoupling[edgeId])')
edit('src/solver/NewtonSolver.cpp','''    if (state.hasReferencedElectronQuasiFermi()) {
        for (int i = 0; i < N; ++i) {''','''    // Subtract physical references before nondimensionalizing even for
    // six-column external states; scaled subtraction loses tiny QF drops.
    for (int i = 0; i < N; ++i) {
        packed(N + i) = static_cast<Real>((static_cast<long double>(state.phin(i)) -
            assembler.electronQuasiFermiReferenceAt(static_cast<Index>(i))) / potentialScale);
        packed(2 * N + i) = static_cast<Real>((static_cast<long double>(state.phip(i)) -
            assembler.holeQuasiFermiReferenceAt(static_cast<Index>(i))) / potentialScale);
    }
    if (state.hasReferencedElectronQuasiFermi()) {
        for (int i = 0; i < N; ++i) {''')
edit('src/solver/NewtonSolver.cpp','''            assembly.poissonEdgeCoupling = value.value(''','''            assembly.poissonChargeNodeVolume = value.value("poisson_charge_node_volume", "inherit");
            assembly.transportEdgeCouplingRatios = value.value(
                "transport_edge_coupling_ratios", std::vector<Real>{});
            assembly.poissonEdgeCoupling = value.value(''')
edit('src/solver/NewtonSolver.cpp','''            if (assembly.transportSignedAverageBoxNodeVolumeScope !=''','''            if (assembly.poissonChargeNodeVolume != "inherit" &&
                assembly.poissonChargeNodeVolume != "signed_transport")
                throw std::invalid_argument("poisson_charge_node_volume must be inherit or signed_transport.");
            if (assembly.poissonChargeNodeVolume != "inherit" &&
                (assembly.poissonElectronTransportNodeVolume || assembly.poissonHoleTransportNodeVolume ||
                 assembly.poissonDopantTransportNodeVolume))
                throw std::invalid_argument("Charge volume policy cannot combine with per-term volume flags.");
            if (!assembly.transportEdgeCouplingRatios.empty() && assembly.transportEdgeCoupling)
                throw std::invalid_argument("Explicit transport ratios cannot combine with transport_edge_coupling.");
            for (Real ratio : assembly.transportEdgeCouplingRatios)
                if (!std::isfinite(ratio) || ratio < 0)
                    throw std::invalid_argument("Transport ratios must be finite and nonnegative.");
            if (assembly.transportSignedAverageBoxNodeVolumeScope !=''')

# Numerically stable analytic psi partials, restricted to unclipped equal-ni kernels.
helper=(R/'scripts/diagnostics/simplemos_stable_sg_psi_derivative.hpp').read_text().replace('simplemos_stable_sg','vela::stable_sg').replace('// Isolated','// Production')
(R/'include/vela/discretization/StableSGDerivative.h').write_text(helper,newline='\n')
edit(name,'#include "vela/equation/CoupledDDAssembler.h"','#include "vela/equation/CoupledDDAssembler.h"\n#include "vela/discretization/StableSGDerivative.h"')
for carrier,qf,offset in [('Electron','phin','electron'),('Hole','phip','hole')]:
    marker=f'                add({qf}Offset() + i, psiOffset() + i, dF_dpsi_i);'
    psi0='psi_i - electronQuantumPotential_V_(i)' if carrier=='Electron' else 'psi_i'
    psi1='psi_j - electronQuantumPotential_V_(j)' if carrier=='Electron' else 'psi_j'
    patch=f'''                // Fixed-mobility partial; mobility chain terms are assembled separately.
                const Real stablePsi0 = {psi0} - {offset}QuasiFermiReferenceAt(idxI);
                const Real stablePsi1 = {psi1} - {offset}QuasiFermiReferenceAt(idxI);
                if (compensatedEqualNiFlux_ && !bgnEnabled_ && niI == niJ &&
                    std::abs((stablePsi0 - {qf}_i) / Vt_) < 500.0 &&
                    std::abs((stablePsi1 - {qf}_j_from_i) / Vt_) < 500.0 &&
                    {'Bminus' if carrier=='Electron' else 'Bplus'} > 0.0) {{
                    const Real stableFlux = sg{carrier}BoltzmannContinuityFlux(
                        niI, niJ, stablePsi0, stablePsi1, {qf}_i, {qf}_j_from_i,
                        Vt_, coef, SGBoltzmannFluxPolicy{{false,true}});
                    const auto derivative = stable_sg::psiDerivative(stableFlux,
                        {'Bminus, dBminusArg' if carrier=='Electron' else 'Bplus, dBplus'}, Vt_, {'true' if carrier=='Electron' else 'false'});
                    dF_dpsi_i = derivative[0]; dF_dpsi_j = derivative[1];
                }}
'''
    edit(name,marker,patch+marker)

header=(R/'scripts/diagnostics/simplemos_linear_refinement.hpp').read_text().replace('simplemos_diagnostic','vela::linear_refinement').replace('// Isolated validation helper. This header is not part of the production solver.','// Optional high-precision residual accumulation for sparse iterative refinement.').replace('isolated refinement','linear refinement')
header=header.replace('    std::vector<MP> r(F.size()), denominator(F.size());','''    if (J.rows() != F.size() || J.cols() != step.size() || !F.allFinite() || !step.allFinite())
        throw std::invalid_argument("Invalid refinement defect dimensions or input");
    std::vector<MP> r(F.size()), denominator(F.size());''')
header=header.replace('    vela::SparseMatrixd scaled=J;','''    if (!rowWeights.allFinite() || (rowWeights.array() <= 0).any())
        throw std::invalid_argument("Refinement row weights must be finite and positive");
    vela::SparseMatrixd scaled=J;''')
(R/'include/vela/solver/LinearRefinement.h').write_text(header,newline='\n')
edit('include/vela/solver/NewtonSolver.h','    Real stallResidualFloor = 1.0e-9;','''    // Zero preserves the direct-solve path; each correction recomputes F+J dx in MP100.
    int linearRefinementIterations = 0;
    Real stallResidualFloor = 1.0e-9;''')
edit('src/solver/NewtonSolver.cpp','#include "vela/solver/NewtonSolver.h"','#include "vela/solver/NewtonSolver.h"\n#include "vela/solver/LinearRefinement.h"')
edit('src/solver/NewtonSolver.cpp','    cfg.stallResidualFloor = json.value("stall_residual_floor", cfg.stallResidualFloor);','''    cfg.linearRefinementIterations = json.value("linear_refinement_iterations", 0);
    if (cfg.linearRefinementIterations < 0 || cfg.linearRefinementIterations > 10)
        throw std::invalid_argument("linear_refinement_iterations must be between 0 and 10.");
    cfg.stallResidualFloor = json.value("stall_residual_floor", cfg.stallResidualFloor);''')
edit('src/solver/NewtonSolver.cpp','''                step = linearSolver.solve(J, -r);
            }
        } catch (const std::runtime_error&) {''','''                step = linearSolver.solve(J, -r);
            }
            if (cfg_.linearRefinementIterations > 0) {
                ScopedPerformanceTimer timer("newton.linear_refinement");
                const VectorXd weights = cfg_.continuityRowScaling.enabled
                    ? activeRowWeights : VectorXd::Ones(r.size());
                step = linear_refinement::refine(J, r, weights, std::move(step),
                    linearSolver, cfg_.linearRefinementIterations,
                    [](int, const VectorXd&, const linear_refinement::Defect&) {});
            }
        } catch (const std::runtime_error&) {''')
