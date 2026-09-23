"""Experimental public DD restart schemas; VDS1 supplies common validation.

HDF5 uses named contiguous datasets; FlatBuffers uses the checked-in schema;
NPY uses a structured little-endian float64 array, never pickle/object arrays.
All formats preserve the independent referenced-QF coordinates.
"""
from pathlib import Path
import struct
import sys

import numpy as np
import h5py
import flatbuffers
sys.path.insert(0, str(Path(__file__).resolve().parent/'generated'))
from vela_state import DDState as fb
import dd_state_binary as binary


def _flags(names):
    for flags in (0,1,3,4,5,7):
        if names == binary.fields(flags):
            return flags
    raise ValueError('Unsupported field order or incomplete restart coordinates')


def _decode(names, count, values):
    flags = _flags(names)
    if count < 1 or count > 2**31-1 or values.size != count*len(names):
        raise ValueError('Invalid state dimensions')
    raw = struct.pack('<8sIIQ',binary.MAGIC,1,flags,count)
    return binary.decode(raw+np.asarray(values,dtype='<f8').tobytes(order='C'))


def write(path, rows):
    path = Path(path)
    raw = binary.encode(rows)
    _,_,flags,count = struct.unpack_from('<8sIIQ',raw)
    names = binary.fields(flags)
    values = np.frombuffer(raw,dtype='<f8',offset=24).reshape(len(names),count)
    if path.suffix == '.h5':
        with h5py.File(path,'x') as f:
            f.attrs['schema'] = 'vela.ddstate/1'
            f.attrs['units'] = 'V;m^-3'
            f.attrs['node_count'] = np.uint64(count)
            f.create_dataset('field_names',data=names,dtype=h5py.string_dtype('utf-8'))
            for name, column in zip(names,values):
                f.create_dataset(name,data=column,dtype='<f8')
    elif path.suffix == '.npy':
        records = np.empty(count,dtype=[(name,'<f8') for name in names])
        for name, column in zip(names,values): records[name] = column
        with path.open('xb') as f: np.save(f,records,allow_pickle=False)
    elif path.suffix == '.vfb':
        builder = flatbuffers.Builder(1024)
        strings = [builder.CreateString(n) for n in names]
        fb.DDStateStartFieldNamesVector(builder,len(names))
        for s in reversed(strings): builder.PrependUOffsetTRelative(s)
        fields = builder.EndVector()
        data = builder.CreateNumpyVector(values.ravel())
        units = builder.CreateString('V;m^-3')
        fb.DDStateStart(builder)
        fb.DDStateAddVersion(builder,1)
        fb.DDStateAddNodeCount(builder,count)
        fb.DDStateAddFieldNames(builder,fields)
        fb.DDStateAddValues(builder,data)
        fb.DDStateAddUnits(builder,units)
        root = fb.DDStateEnd(builder)
        builder.Finish(root,file_identifier=b'VDSF')
        with path.open('xb') as f: f.write(builder.Output())
    else: raise ValueError('Unsupported public state extension')


def read(path):
    path = Path(path)
    if path.stat().st_size > 2**30: raise ValueError('State exceeds size limit')
    if path.suffix == '.h5':
        with h5py.File(path,'r') as f:
            def attr(name):
                v = f.attrs[name]
                return v.decode() if isinstance(v,bytes) else v
            if attr('schema') != 'vela.ddstate/1' or attr('units') != 'V;m^-3':
                raise ValueError('Unsupported HDF5 schema or units')
            names = list(f['field_names'].asstr()[...]);_flags(names)
            count = int(f.attrs['node_count'])
            if count < 1 or count*len(names)*8 > 2**30: raise ValueError('Invalid node count')
            columns = []
            for name in names:
                d = f[name]
                if d.shape != (count,) or d.dtype.kind != 'f' or d.dtype.itemsize != 8:
                    raise ValueError('Invalid HDF5 field dtype or shape')
                columns.append(d[...])
            return _decode(names,count,np.asarray(columns))
    if path.suffix == '.npy':
        records = np.load(path,allow_pickle=False)
        names = list(records.dtype.names or []);_flags(names)
        if records.ndim != 1 or records.dtype != np.dtype([(n,'<f8') for n in names]):
            raise ValueError('Unsupported NPY dtype or shape')
        return _decode(names,len(records),np.asarray([records[n] for n in names]))
    if path.suffix == '.vfb':
        raw = path.read_bytes()
        if len(raw)<8 or raw[4:8]!=b'VDSF': raise ValueError('Invalid FlatBuffers identifier')
        try:
            s = fb.DDState.GetRootAs(raw)
            if s.Version()!=1 or s.Units()!=b'V;m^-3': raise ValueError('Unsupported FlatBuffers schema')
            if not 5<=s.FieldNamesLength()<=11: raise ValueError('Invalid fields')
            names = [s.FieldNames(i).decode() for i in range(s.FieldNamesLength())]
            _flags(names)
            count = s.NodeCount()
            if count<1 or s.ValuesLength()!=count*len(names): raise ValueError('Invalid values vector')
            return _decode(names,count,s.ValuesAsNumpy())
        except (IndexError,TypeError,struct.error,OverflowError) as exc:
            raise ValueError('Malformed FlatBuffers state') from exc
    raise ValueError('Unsupported public state extension')
