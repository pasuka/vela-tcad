"""VDS1 binary64 restart interchange: SI densities, volts, implicit node IDs."""
from array import array
import math
from pathlib import Path
import struct
import sys

MAGIC=b'VELADS01'
BASE=['psi','phin','phip','electrons_m3','holes_m3']
QF=['electron_qf_increment_V','hole_qf_increment_V','electron_qf_reference_V','hole_qf_reference_V']

def fields(flags):
    if flags&~7 or flags&2 and not flags&1:raise ValueError('Unsupported VDS1 flags')
    return BASE+(['electron_quantum_potential_V'] if flags&1 else [])+(['electron_quantum_potential_like_V'] if flags&2 else [])+(QF if flags&4 else [])

def decode(raw):
    if len(raw)<24:raise ValueError('Truncated VDS1 header')
    magic,version,flags,count=struct.unpack_from('<8sIIQ',raw)
    if magic!=MAGIC or version!=1:raise ValueError('Unsupported VDS1 header')
    names=fields(flags)
    if count>2**31-1 or len(raw)!=24+count*len(names)*8:raise ValueError('VDS1 length mismatch')
    values=array('d');values.frombytes(raw[24:])
    if sys.byteorder!='little':values.byteswap()
    if not all(map(math.isfinite,values)):raise ValueError('Non-finite state')
    columns=[values[i*count:(i+1)*count] for i in range(len(names))]
    result=[dict(zip(['node_id']+names,[str(i)]+list(row))) for i,row in enumerate(zip(*columns))]
    if flags&4:
        for r in result:
            for physical,ref,inc in [('phin',QF[2],QF[0]),('phip',QF[3],QF[1])]:
                value=r[ref]+r[inc]
                if not math.isfinite(value) or abs(r[physical]-value)>32*sys.float_info.epsilon*max(1.,abs(r[physical]),abs(value)):
                    raise ValueError('Inconsistent referenced coordinates')
    return result

def encode(rows):
    if not rows:raise ValueError('Empty state')
    keys=set(rows[0]);flags=(1 if 'electron_quantum_potential_V' in keys else 0)|(2 if 'electron_quantum_potential_like_V' in keys else 0)
    if keys.intersection(QF):
        if not set(QF)<=keys:raise ValueError('Partial referenced coordinates')
        flags|=4
    names=fields(flags)
    if keys!={'node_id',*names}:raise ValueError('Unsupported state fields')
    for i,row in enumerate(rows):
        if set(row)!=keys or int(row['node_id'])!=i:raise ValueError('Noncanonical nodes or fields')
    data=array('d',(float(r[n]) for n in names for r in rows))
    for i,v in enumerate(data):
        if not math.isfinite(v):raise ValueError('Non-finite state')
        if v!=0 and abs(v)<sys.float_info.min:data[i]=0.
    if sys.byteorder!='little':data.byteswap()
    raw=struct.pack('<8sIIQ',MAGIC,1,flags,len(rows))+data.tobytes()
    # Structural checks above are cheap; coordinate consistency is checked on read.
    return raw

def read(path):return decode(Path(path).read_bytes())
def write(path,rows):
    raw=encode(rows)
    with Path(path).open('xb') as f:f.write(raw)
