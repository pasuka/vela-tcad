"""Offline SDE package entry point; Gmsh is an explicit external runtime."""
import hashlib
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
PACKAGE_INTEGRAL_TARGET = 0.00025
PACKAGE_INTEGRAL_MAX_ITERATIONS = 16


def verify():
    manifest = json.loads((ROOT / 'PACKAGE_MANIFEST.json').read_text(encoding='utf-8'))
    for item in manifest['files']:
        path = (ROOT / item['path']).resolve()
        if not path.is_relative_to(ROOT) or not path.is_file():
            raise ValueError('missing or escaping package file: '+item['path'])
        if hashlib.sha256(path.read_bytes()).hexdigest() != item['sha256']:
            raise ValueError('package hash mismatch: '+item['path'])
    return manifest


def main():
    sys.dont_write_bytecode = True
    try:
        manifest = verify()
        if sys.argv[1:] == ['--version']:
            print(json.dumps({k:manifest[k] for k in ('schema','platform','python_version','gmsh_profile','source_commit')}))
            return 0
        # The optional external CLI can itself be a Python script. Do not
        # pass host Python injection variables into that child process.
        os.environ.pop('PYTHONPATH',None)
        os.environ.pop('PYTHONHOME',None)
        sys.path[:0] = [str(ROOT/'vendor'),str(ROOT/'scripts')]
        from sentaurus_gmsh_pipeline import build_parser, main as pipeline_main
        arguments = sys.argv[1:]
        args = build_parser().parse_args(arguments)
        provided_options = {arg.split('=', 1)[0] for arg in arguments}
        components = [args.gmsh_python or os.environ.get('VELA_GMSH_PYTHON'),
                      args.gmsh_library or os.environ.get('VELA_GMSH_LIBRARY'),
                      args.gmsh_executable or os.environ.get('VELA_GMSH_EXECUTABLE')]
        if not all(components):
            raise ValueError('configure explicit Gmsh 4.15.2 wrapper, shared library and executable paths')
        handles = []
        for component in components[1:]:
            directory = str(Path(component).resolve().parent)
            os.environ['PATH'] = directory+os.pathsep+os.environ.get('PATH','')
            if os.name == 'nt':
                handles.append(os.add_dll_directory(directory))
        if '--runtime-dir' not in provided_options:
            arguments += ['--runtime-dir',str(args.output_dir.resolve()/'.gmsh-runtime')]
        if args.capability_level in {'M3','M5'}:
            integral_target = args.integral_relative_target
            if '--integral-relative-target' not in provided_options:
                integral_target = PACKAGE_INTEGRAL_TARGET
                arguments += ['--integral-relative-target',str(integral_target)]
            # Integral splits can introduce new edges that violate an authored
            # size rule even after the integral error has met its target.
            if integral_target > 0:
                arguments = ['--maxtransdiff-max-iterations', str(PACKAGE_INTEGRAL_MAX_ITERATIONS)] + arguments
        return pipeline_main(arguments)
    except (ValueError,OSError,ImportError) as error:
        print('vela-sde: '+str(error),file=sys.stderr)
        return 2


if __name__ == '__main__':
    raise SystemExit(main())
