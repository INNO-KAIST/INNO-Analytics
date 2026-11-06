"""Independently reparse the JS writer's ZIP/CSV using Python stdlib."""
import csv
import io
import json
import math
import pathlib
import subprocess
import zipfile

root = pathlib.Path(__file__).resolve().parents[3]
node = r"""const io=require('./modules/image/src/img_io.js');
const values=new Float64Array([0,-0,1+2**-52,-3.25,NaN,Infinity,-Infinity,1e-300,Number.MAX_VALUE]);
const p={name:'../경로/one',w:3,h:3,data:values,cal:{pw:0.2,ph:0.3,unit:'µm'},valueUnit:'nm',c:0,z:2,t:1};
process.stdout.write(Buffer.from(io.exportAfmCsvBundle({format:process.argv[1],planes:[p,p]}).bytes).toString('base64'));"""
import base64
expected=[0.0,-0.0,1+2**-52,-3.25,math.nan,math.inf,-math.inf,1e-300,float.fromhex('0x1.fffffffffffffp+1023')]
for fmt in ('matrix','coordinates'):
    encoded=subprocess.run(['node','-e',node,fmt],cwd=root,check=True,capture_output=True,text=True).stdout
    with zipfile.ZipFile(io.BytesIO(base64.b64decode(encoded))) as archive:
        assert archive.testzip() is None
        assert len(archive.namelist())==len(set(archive.namelist()))==3
        metadata=json.loads(archive.read('metadata.json'))
        for plane in metadata['planes']:
            assert '/' not in plane['file'] and '\\' not in plane['file']
            assert plane['xy']['unit']=='µm'
            rows=list(csv.reader(io.StringIO(archive.read(plane['file']).decode('utf-8'))))
            if fmt=='matrix': cells=[x for row in rows for x in row]
            else:
                assert rows.pop(0)==['x','y','value']
                assert [float(x) for x in rows[4][:2]]==[0.2,0.3]
                cells=[r[2] for r in rows]
            assert len(cells)==len(expected)
            for got,want in zip(map(float,cells),expected):
                assert (math.isnan(got) and math.isnan(want)) or (got==want and (want!=0 or math.copysign(1,got)==math.copysign(1,want)))
print('Python zipfile CRC/central directory and CSV numeric roundtrip passed')
