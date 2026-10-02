"""Exclude adjoint-amplified nonlinear residual as the channel FD mismatch."""
import math
import validate_simplemos_vela_channel_charge as run


def audit():
    run.audit.verify(run.FREEZE)
    adj = {int(r['node_id']):r for r in run.audit.rows(run.native.upstream.LOCAL/'m65_n23_vd_0p050000_endpoint/adjoint.csv')}
    residuals = {}
    for name, _ in run.CASES:
        dest = run.LOCAL/name
        cfg = run.audit.read(dest/'config.json'); cfg.pop('output_state_file')
        cfg.update(state_file=str(dest/'state.csv'),simulation_type='newton_residual_probe',output_csv=str(dest/'final_residual.csv'))
        status = run.execute(cfg,dest/'final_residual.json',dest/'charge.csv')
        if status['exit_code'] != 0: raise ValueError('Residual probe failed')
        residuals[name] = {int(r['node_id']):r for r in run.audit.rows(dest/'final_residual.csv')}
    grad = run.audit.rows(run.OUT/'actual_direction_gradient.csv'); out = []
    for suffix, amp in (('1e13',1e13),('5e12',5e12)):
        p,m = residuals['plus_'+suffix],residuals['minus_'+suffix]
        parts = {}
        for field,weight in (('psi_residual','lambda_poisson'),('phin_residual','lambda_electron'),('phip_residual','lambda_hole')):
            parts[weight] = math.fsum(float(a[weight])*(float(p[i][field])-float(m[i][field]))/2 for i,a in adj.items())
        remainder = math.fsum(parts.values())
        g = next(r for r in grad if float(r['amplitude_cm_3'])==amp)
        discrepancy = float(g['gradient_minus_adjoint_A_per_um'])
        out.append({'amplitude_cm_3':amp,'gradient_minus_adjoint_A_per_um':discrepancy,
            'lambda_dot_charged_residual_delta_A_per_um':remainder,
            'residual_fraction_of_prediction_gap':abs(remainder/discrepancy),
            'inferred_lambda_dot_J_delta_x_minus_F_delta_A_per_um':discrepancy-remainder,**parts})
    run.audit.write_csv(run.OUT/'stationarity_projection.csv',out)
    for r in out: print(r)


if __name__=='__main__':
    audit()
