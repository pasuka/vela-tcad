"""Resolve copied evidence without rewriting immutable run JSON files."""
from pathlib import Path, PurePosixPath


def candidate_path(value, mapping=None):
    if mapping is None:
        return Path(value)
    source, destination = mapping
    source_path, value_path = PurePosixPath(source), PurePosixPath(value)
    if not source_path.is_absolute() or '..' in source_path.parts or '..' in value_path.parts:
        raise ValueError('Evidence mapping requires an absolute source and no parent traversal')
    try:
        relative = value_path.relative_to(source_path)
    except ValueError as error:
        raise ValueError('Candidate evidence is outside the declared source root') from error
    return Path(destination).joinpath(*relative.parts)
