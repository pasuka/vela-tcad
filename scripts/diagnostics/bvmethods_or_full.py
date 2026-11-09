"""Independent OR qualification on the committed ablation supervisor.

Uses the same experimental C++ executable. This controller explicitly qualifies OR accepted states
and never changes an existing run, source tree, or production default.
"""
import json
import math
from pathlib import Path

from bvmethods_full_ablation import Supervisor, accepted, put


def finite(value):
    try:
        value = float(value)
        return value if math.isfinite(value) else None
    except (TypeError, ValueError):
        return None


def or_gate(row):
    g = json.loads(row['validation_diagnostics'])
    f = finite(g['coupled_device_residual_inf'])
    load = finite(g['actual_load_line_residual_V'])
    rms = finite(g['density_update_rms'])
    residual = f is not None and f >= 0 and f <= 1e-8
    update = g['density_update_valid'] and rms is not None and 0 <= rms < 1
    return bool(g['density_update_mode'] == 'enforce'
                and g['residual_stopping_scope'] == 'device_and_circuit'
                and g['load_line_enforced']
                and g['load_line_convergence_mode'] == 'enforce'
                and g['physical_residuals_finite']
                and g['terminal_kcl_satisfied']
                and g['actual_load_line_satisfied']
                and load is not None and abs(load) <= 1e-6
                and f is not None and f >= 0
                and (residual or update)
                and (not g['carrier_row_enforced'] or g['carrier_row_satisfied'])
                and (not g['global_continuity_enforced'] or g['global_continuity_satisfied'])
                and all(finite(row[k]) is not None for k in
                        ['outer_voltage_V', 'inner_voltage_V', 'current_total_A_per_um']))


def gate_summary(data):
    result = {'residual_only': 0, 'update_only': 0, 'both': 0,
              'zero_update_points': 0, 'failed_points': 0,
              'carrier_report_unsatisfied_points': 0,
              'continuity_report_unsatisfied_points': 0,
              'max_device_residual': 0., 'max_load_residual_V': 0.,
              'max_valid_update_rms': None}
    for row in data:
        g = json.loads(row['validation_diagnostics'])
        f, load, rms = (finite(g[k]) for k in
                        ['coupled_device_residual_inf', 'actual_load_line_residual_V', 'density_update_rms'])
        residual = f is not None and f <= 1e-8 and load is not None and abs(load) <= 1e-6
        update = g['density_update_valid'] and rms is not None and 0 <= rms < 1
        if residual and update:
            result['both'] += 1
        elif residual:
            result['residual_only'] += 1
        elif update:
            result['update_only'] += 1
        result['failed_points'] += not or_gate(row)
        result['zero_update_points'] += int(row['iterations']) == 0
        result['carrier_report_unsatisfied_points'] += not g['carrier_row_enforced'] and not g['carrier_row_satisfied']
        result['continuity_report_unsatisfied_points'] += not g['global_continuity_enforced'] and not g['global_continuity_satisfied']
        if f is not None:
            result['max_device_residual'] = max(result['max_device_residual'], f)
        if load is not None:
            result['max_load_residual_V'] = max(result['max_load_residual_V'], abs(load))
        if g['density_update_valid'] and rms is not None:
            result['max_valid_update_rms'] = max(result['max_valid_update_rms'] or 0., rms)
    return result


class ORSupervisor(Supervisor):
    def config(self, template, folder, group):
        if group not in ('G1', 'G2'):
            raise ValueError('expected G1 Full_OR or G2 Automatic_guard_OR')
        cfg = super().config(template, folder, group)
        cfg['sweep']['external_circuit']['coupled_update_convergence']['mode'] = 'enforce'
        return cfg

    def sweep(self, name, template, group, **kwargs):
        result = super().sweep(name, template, group, **kwargs)
        folder = kwargs.get('folder') or self.out/name/group
        data = accepted(folder/'iv.csv')
        result['all_accepted_gates'] = all(or_gate(row) for row in data)
        result['OR_gate_summary'] = gate_summary(data)
        result['passed'] = bool(result['solver_exit_code'] == 0 and not result['timed_out']
                                and result['all_accepted_gates'] and result['requested_reached'])
        result['stopping_combination'] = 'OR'
        put(folder/'qualification.json', result)
        return result

    def p1(self, group):
        result = super().p1(group)
        result['short_window_target_test'] = result['passed']
        result['passed'] = result['full_bv_completed']
        result['qualification_scope'] = 'complete_BV_OR_with_original_independent_gates'
        folder = self.out/'P1'/group
        put(folder/'status.json', result)
        put(folder/'qualification.json', result)
        self.save()
        return result
