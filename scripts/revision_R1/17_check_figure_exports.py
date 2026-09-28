"""Normalize TIFF technical metadata without resampling or changing RGB pixels."""
from pathlib import Path
from PIL import Image
import csv
root=Path(__file__).resolve().parents[2];out=root/'results/revision_R1/figures';rows=[]
for p in sorted(out.glob('*.tif')):
 with Image.open(p) as im:
  im.load();rgb=im.convert('RGB');before=rgb.tobytes()
  if im.mode!='RGB' or im.info.get('dpi')!=(600.0,600.0) or im.info.get('compression')!='tiff_lzw':
   if im.mode=='RGBA':assert im.getchannel('A').getextrema()==(255,255),'Unexpected transparency requires explicit handling'
   temp=p.with_suffix('.normalized.tif');rgb.save(temp,dpi=(600,600),compression='tiff_lzw')
   with Image.open(temp) as chk:assert chk.tobytes()==before
   temp.replace(p)
 with Image.open(p) as im:
  rows.append(dict(file=p.name,width_px=im.width,height_px=im.height,dpi=str(im.info['dpi']),mode=im.mode,bytes=p.stat().st_size))
with (root/'results/revision_R1/tables/figure_parameters.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
print('TIFF metadata verified:',len(rows),'files; RGB pixels preserved.')
