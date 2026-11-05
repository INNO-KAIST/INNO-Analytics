import struct, pathlib
D=pathlib.Path(__file__).resolve().parent/'curve_fixtures';D.mkdir(exist_ok=True)
def wave(filename,*,curve=False,labels=('Bias','Current','Aux'),keys=('ForceDist','ForceScanRate','NumPtsPerWave'),cal=False,three=False,truncate=False):
 nx,ny,nl=4,3,(2 if three else 0); n=nx*ny*(nl or 1); note=''.join(k+': 1\r' for k in keys).encode(); label=(b'\0'*32+b''.join(v.encode().ljust(32,b'\0') for v in labels)) if curve or labels else b''; dim=2 if three else 1
 data=bytearray(384+4*n+len(note)+len(label)); struct.pack_into('<H',data,0,5); struct.pack_into('<i',data,12,len(note));struct.pack_into('<i',data,36+4*dim,len(label));struct.pack_into('<i',data,76,n);struct.pack_into('<H',data,80,2);struct.pack_into('<4i',data,132,nx,ny,nl,0)
 if cal:struct.pack_into('<2d',data,148,1e-9,2e-9);data[216:218]=b'm\0';data[220:222]=b'm\0'
 for i in range(n):struct.pack_into('<f',data,384+4*i,i+0.5)
 data[384+4*n:384+4*n+len(note)]=note;data[384+4*n+len(note):]=label
 if truncate:data=data[:-3]
 (D/filename).write_bytes(data)
wave('curve.ibw',curve=True)
wave('casefold.ibw',curve=True,keys=('forcedist','FORCESCANRATE','NumPtsPerWave'))
for key in ('ForceDist','ForceScanRate','NumPtsPerWave'):wave('missing_'+key+'.ibw',curve=True,keys=tuple(x for x in ('ForceDist','ForceScanRate','NumPtsPerWave') if x!=key))
wave('emptylabel.ibw',curve=True,labels=('Bias','','Aux'))
wave('broken.ibw',curve=True,truncate=True)
wave('calibrated.ibw',labels=(),keys=(),cal=True)
wave('three.ibw',labels=('Height','Phase'),keys=('ForceDist','ForceScanRate','NumPtsPerWave'),three=True,cal=True)

wave('calibrated_strong.ibw',curve=True,cal=True)
wave('calibrated_partial.ibw',curve=True,cal=True,keys=('ForceDist',))
wave('malformed_block.ibw',curve=True)
p=D/'malformed_block.ibw'; data=bytearray(p.read_bytes());struct.pack_into('<i',data,40,127);p.write_bytes(data)
wave('calibrated_key_only.ibw',labels=(),keys=('ForceDist',),cal=True)
wave('calibrated_labels_only.ibw',curve=True,keys=(),cal=True)
