#!/usr/bin/env python3
"""Audit an installed-library resistor acquisition, independently of E2 WP14.

The profile identifies one native run, not Vela model equivalence, branch
stability, repeatability or threshold-field qualification. Commercial reference
inputs are supplied explicitly and are never copied into repository fixtures.
"""
from __future__ import annotations
import argparse
from decimal import Decimal
import hashlib
import json
from pathlib import Path
import re

from analyze_sentaurus_bvmethods import first_upward_crossing
from vela_sdevice.cmd_parser.lexer import Source, strip_comments
from vela_sdevice.cmd_parser.syntax import group_blocks
from vela_sdevice.plt.document import read_plt_document
from validate_wp14_zero_qs import block, trajectory_checks

STOP = 1.3e-3
RESISTANCE = 1e7
SCHEMA = 'vela.bvmethods.installed_library_resistor_reference.v1'
# Authored declaration contract, not a preprocessed commercial command file.
PHYSICS = '''EffectiveIntrinsicDensity(OldSlotboom)
Mobility(DopingDep HighFieldsaturation(GradQuasiFermi) Enormal)
Recombination(SRH(DopingDep) Band2Band(Model=NonLocalPath) Avalanche(GradQuasiFermi))
Fermi'''
ELECTRODES = '''{Name="drain" Voltage=0 Resistor=1e7}
{Name="source" Voltage=0} {Name="gate" Voltage=0 Barrier=-.55}
{Name="substrate" Voltage=0}'''
MATH = '''Extrapolate Iterations=20 Notdamped=100 RelErrControl AvalDerivatives
ErrRef(Electron)=1e10 ErrRef(Hole)=1e10
BreakCriteria{Current(Contact="drain" AbsVal=1.3e-3)} Transient=BE'''
SOLVE = '''Coupled(Iterations=100){Poisson} Coupled{Poisson Electron Hole}
Quasistationary(InitialStep=1e-7 Increment=1.41 MinStep=1e-10 MaxStep=.025
Goal{Name="drain" Voltage=100000}){Coupled{Poisson Electron Hole}}'''
TOKEN = re.compile(r'"[^"\n]*"|[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?|[A-Za-z_]\w*|\S')
NUMBER = re.compile(r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?\Z')


def canonical(text: str) -> tuple:
    """Retain every non-whitespace token; fold numeric spelling only."""
    return tuple(Decimal(t) if NUMBER.fullmatch(t) else t
                 for t in TOKEN.findall(strip_comments(text)))


def source_checks(text: str) -> dict:
    source = Source(text, 'reference.cmd')
    sections, issue = group_blocks(source.tokens(), source)
    names = [s.name.value.lower() for s in sections]
    checks = {'single_device_sections': issue is None and sorted(names) ==
              sorted(['file', 'electrode', 'physics', 'plot', 'math', 'solve']),
              'no_preprocessor': not source.preprocessor_findings}
    clean = strip_comments(text)
    for name, expected in [('Physics', PHYSICS), ('Electrode', ELECTRODES),
                           ('Math', MATH), ('Solve', SOLVE)]:
        checks['original_' + name.lower()] = canonical(block(clean, name)) == canonical(expected)
    return checks


def stop_record_checks(log: str, actual_last_current: float) -> dict:
    records = re.findall(r'Contact current exceeded limit:\s*electrode:\s*(\S+)\s*'
                         r'value\s*:\s*([-+\d.eE]+)\s*limit\s*:\s*([-+\d.eE]+)', log)
    if len(records) != 1:
        return {'single_drain_stop_record': False, 'declared_stop_matches': False,
                'reported_value_matches_last': False}
    name, value, limit = records[0]
    value, limit = float(value), float(limit)
    # The native stop report uses five significant digits; PLT has more.
    return {'single_drain_stop_record': name == 'drain',
            'declared_stop_matches': limit == STOP,
            'reported_value_matches_last': abs(abs(actual_last_current)-abs(value)) <=
            max(1e-20, abs(value)*5e-5)}


def validate(cmd: Path, par: Path, mesh: Path, plt: Path, log: Path,
             reference_cmd: Path, reference_par: Path, expected_mesh_sha256: str) -> dict:
    text = cmd.read_text(encoding='utf-8');reference = reference_cmd.read_text(encoding='utf-8')
    checks = source_checks(text)
    checks.update({'reference_' + k: v for k, v in source_checks(reference).items()})
    checks['parameters_match_frozen_reference'] = canonical(par.read_text(encoding='utf-8')) == canonical(reference_par.read_text(encoding='utf-8'))
    clean_file = block(strip_comments(text), 'File')
    for field, path in [('Grid', mesh), ('Parameter', par), ('Current', plt), ('Output', log)]:
        values = re.findall(r'\b'+field+r'\s*=\s*"([^"\n]+)"', clean_file)
        checks['file_' + field.lower() + '_bound'] = values == [path.name]
    sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    checks['mesh_matches_frozen_hash'] = sha(mesh) == expected_mesh_sha256
    names, rows = read_plt_document(plt.read_text(encoding='utf-8'))
    if set(names) != {'drain', 'source', 'gate', 'substrate'}:
        raise ValueError('all four contacts, and no additional contacts, required')
    log_text = log.read_text(encoding='utf-8')
    checks['native_version_matches'] = bool(re.search(r'\bVersion\s+T-2022\.03-SP2\b', log_text))
    trajectory = trajectory_checks(rows, log_text, stop=STOP, resistance=RESISTANCE)
    checks.update(stop_record_checks(log_text, rows[-1]['drain']['total_current_A_per_um']))
    currents = [abs(r['drain']['total_current_A_per_um']) for r in rows]
    inner = [r['drain']['inner_voltage_V'] for r in rows]
    bv = first_upward_crossing(inner, currents, 1e-4)
    checks['BV_criterion_bracketed'] = bv is not None
    return {'schema': SCHEMA, 'profile': 'installed_T2022_resistor',
            'status': 'pass' if all(checks.values()) and trajectory['status']=='pass' else 'fail',
            'source_and_binding_checks': checks, 'trajectory': trajectory,
            'BV': {'current_A_per_um': 1e-4, 'voltage_V': bv,
                   'method': 'first upward crossing, linear interpolation in current',
                   'native_visualizer_executed': False},
            'inputs': {k: {'name': p.name, 'sha256': sha(p)} for k,p in
                       [('cmd',cmd),('par',par),('mesh',mesh),('plt',plt),('log',log),
                        ('reference_cmd',reference_cmd),('reference_par',reference_par)]},
            'qualified_scope': 'one original-model native from-zero QS acquisition only',
            'repeated_runs_qualified': False, 'Vela_equivalence_qualified': False,
            'threshold_field_state_qualified': False}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ['cmd','par','mesh','plt','log','reference-cmd','reference-par','output']:
        parser.add_argument('--'+name, type=Path, required=True)
    parser.add_argument('--expected-mesh-sha256', required=True)
    args = vars(parser.parse_args());output=args.pop('output')
    try:
        result = validate(**args)
    except (OSError, ValueError, KeyError, IndexError) as exc:
        result = {'schema':SCHEMA, 'status':'fail', 'input_error':str(exc)}
    output.parent.mkdir(parents=True,exist_ok=True)
    output.write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print('original-library resistor reference: '+result['status'])
    return 0 if result['status']=='pass' else 1


if __name__ == '__main__':
    raise SystemExit(main())
