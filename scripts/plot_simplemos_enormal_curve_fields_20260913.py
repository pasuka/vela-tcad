"""Scientific field-error curves from the qualified 204-point reconstruction."""
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import simplemos_enormal_curves_20260912 as e

a, d = e.a, e.d
O = e.O / 'fields'


def main():
    a.verify(O / 'evidence.json')
    assert a.read(O / 'summary.json')['all_reconstruction_qualified']
    assert a.read(e.O / 'summary.json')['comparison_qualified'] == 204
    fields = a.rows(O / 'fields.csv')
    srh = a.rows(O / 'srh.csv')
    specs = [('ElectrostaticPotential', 'Potential: max |Vela - native|', 'V'),
             ('eQuasiFermiPotential', 'Electron quasi-Fermi: max difference', 'V'),
             ('hQuasiFermiPotential', 'Hole quasi-Fermi: max difference', 'V'),
             ('eDensity', 'Electron density: max |log10(Vela/native)|', 'dex'),
             ('hDensity', 'Hole density: max |log10(Vela/native)|', 'dex')]
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=True)
    maxima = []
    for device, color in [('n19', '#1768ac'), ('n23', '#d04b35')]:
        for vd, style in [(.05, '-'), (1., '--')]:
            label = f'{device}, Vd={vd:g} V'
            for ax, (field, title, unit) in zip(axes.flat, specs):
                rows = sorted([r for r in fields if r['device'] == device and
                               float(r['vd']) == vd and r['scope'] == 'all_Si' and r['field'] == field],
                              key=lambda r: float(r['vg']))
                assert len(rows) == 51 and len({r['vg'] for r in rows}) == 51
                assert all(r['units'] == unit for r in rows)
                ax.semilogy([float(r['vg']) for r in rows], [float(r['max_abs']) for r in rows],
                            color=color, linestyle=style, label=label)
                worst = max(rows, key=lambda r: float(r['max_abs']))
                maxima.append(dict(device=device, vd=vd, field=field, units=unit,
                                   max_abs=float(worst['max_abs']), vg=float(worst['vg']),
                                   node=int(worst['max_node'])))
                ax.set_title(title, fontsize=10)
                ax.set_ylabel(unit)
            rows = sorted([r for r in srh if r['device'] == device and float(r['vd']) == vd],
                          key=lambda r: float(r['vg']))
            assert len(rows) == 51 and len({r['vg'] for r in rows}) == 51
            axes[1, 2].semilogy([float(r['vg']) for r in rows], [float(r['normalized_L1']) for r in rows],
                                color=color, linestyle=style, label=label)
            worst = max(rows, key=lambda r: float(r['normalized_L1']))
            maxima.append(dict(device=device, vd=vd, field='SRH_common_Si_volume_L1', units='relative',
                               max_abs=float(worst['normalized_L1']), vg=float(worst['vg']), node=''))
    axes[1, 2].set_title('SRH rate: normalized L1 (common Si volume)', fontsize=10)
    axes[1, 2].set_ylabel('relative')
    for ax in axes.flat:
        ax.grid(alpha=.25, which='both')
        ax.set_xlim(0, 1)
    for ax in axes[1]:
        ax.set_xlabel('Vg (V)')
    axes[0, 0].legend(fontsize=8)
    fig.suptitle('PhuMob + Enormal: Vela versus Sentaurus, all Si nodes / 204 qualified points')
    fig.text(.5, .012, 'Field diagnostics, not dual-initialization errors. SRH uses common positive Si volumes; it is not a drain-current error.',
             ha='center', fontsize=9)
    fig.tight_layout(rect=(0, .04, 1, .95))
    for extension in ('png', 'pdf'):
        fig.savefig(O / ('field_error_curves.' + extension), dpi=180)
    plt.close(fig)
    a.write_csv(O / 'curve_maxima.csv', maxima)
    d.matrix.freeze(O / 'visual_evidence.json',
                    [Path(__file__).resolve(), O / 'evidence.json', O / 'field_error_curves.png',
                     O / 'field_error_curves.pdf', O / 'curve_maxima.csv'])
    print(maxima, flush=True)


if __name__ == '__main__':
    main()
