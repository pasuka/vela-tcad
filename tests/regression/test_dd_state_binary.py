import math
from pathlib import Path
import struct
import sys
import tempfile
import unittest
sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'scripts'))
import dd_state_binary as b
from ldmos_binary_adapter import translated
from translate_dd_state import translate_state_csv
import csv

class BinaryStateTest(unittest.TestCase):
    def state(self):
        return [dict(node_id=str(i),psi=.1*i,phin=28.,phip=-28.,electrons_m3=1e23,holes_m3=2e9,
            electron_qf_reference_V=28.,hole_qf_reference_V=-28.,
            electron_qf_increment_V=1e-17,hole_qf_increment_V=-1e-17) for i in range(3)]
    def test_normal_values_and_low_coordinate_bits_roundtrip(self):
        s=self.state();self.assertEqual(b.decode(b.encode(s)),s)
        s[0]['electron_qf_increment_V']=4*math.ulp(0.)
        self.assertEqual(b.decode(b.encode(s))[0]['electron_qf_increment_V'],0.)
    def test_independent_wire_encoding_and_node_order(self):
        raw=struct.pack('<8sIIQ5d',b.MAGIC,1,0,1,.1,.2,.3,1e23,2e9)
        self.assertEqual(b.decode(raw)[0]['phin'],.2)
        self.assertEqual(b.encode(b.decode(raw)),raw)
        s=self.state();s[1]['node_id']='0'
        with self.assertRaises(ValueError):b.encode(s)
    def test_invalid_headers_lengths_and_nonfinite_fields(self):
        raw=b.encode(self.state())
        for bad in [raw[:10],raw[:-1],raw+b'x',raw[:8]+struct.pack('<I',2)+raw[12:],raw[:12]+struct.pack('<I',8)+raw[16:],raw[:24]+struct.pack('<d',float('nan'))+raw[32:]]:
            with self.assertRaises(ValueError):b.decode(bad)
        s=self.state();s[0]['electron_qf_reference_V']=1.
        with self.assertRaises(ValueError):b.decode(b.encode(s))
    def test_translation_matches_csv_and_keeps_inactive_placeholders(self):
        with tempfile.TemporaryDirectory() as d:
            d=Path(d);s=self.state();p=d/'in.csv'
            with p.open('w',newline='') as f:
                w=csv.DictWriter(f,fieldnames=list(s[0]));w.writeheader();w.writerows(s)
            translate_state_csv(p,d/'out.csv',28.,[0,1])
            with (d/'out.csv').open() as f:csv_rows=list(csv.DictReader(f))
            result=b.decode(b.encode(translated(s,28.,{0,1})))
            for a,c in zip(result,csv_rows):
                for key in a:self.assertEqual(float(a[key]),float(c[key]))

if __name__=='__main__':unittest.main()
