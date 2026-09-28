from pathlib import Path
import urllib.request,json,hashlib,concurrent.futures,time
root=Path(__file__).resolve().parents[2]; out=root/'data/revision_R1'; out.mkdir(parents=True,exist_ok=True)
tasks=[]
for study,expr in [('brca_metabric','data_mrna_illumina_microarray.txt'),('brca_tcga_pan_can_atlas_2018','data_mrna_seq_v2_rsem.txt')]:
 for name in ['data_clinical_patient.txt','data_clinical_sample.txt',expr,'meta_'+expr[5:]]:
  base=('https://raw.githubusercontent.com/cBioPortal/datahub/master/public/' if name.startswith('meta_') else 'https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/')
  tasks.append((study+'/'+name,base+study+'/'+name))
tasks += [('MCPcounter.R','https://raw.githubusercontent.com/ebecht/MCPcounter/master/Source/R/MCPcounter.R'),('MCP_genes.txt','https://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt'),('MCP_probesets.txt','https://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/probesets.txt'),('GSE176078_family.soft.gz','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE176nnn/GSE176078/soft/GSE176078_family.soft.gz')]
for name in ['rankAndScoring.R','rankGenesGeneric.R','simpleScoreGeneric.R','singscore.R']:
 tasks.append(('singscore_source/'+name,'https://raw.githubusercontent.com/DavisLaboratory/singscore/master/R/'+name))

def get(t):
 name,url=t;p=out/name;p.parent.mkdir(parents=True,exist_ok=True)
 if p.exists() and p.stat().st_size>200:return dict(file=name,url=url,status='cached',bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
 try:
  with urllib.request.urlopen(url,timeout=45) as r, p.with_suffix(p.suffix+'.partial').open('wb') as f:
   while True:
    b=r.read(1024*1024)
    if not b:break
    f.write(b)
  p.with_suffix(p.suffix+'.partial').rename(p)
  row=dict(file=name,url=url,status='downloaded',bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
 except Exception as e:row=dict(file=name,url=url,status='failed',error=str(e))
 print(json.dumps(row),flush=True);return row
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as ex:rows=list(ex.map(get,tasks))
(out/'download_manifest.json').write_text(json.dumps(rows,indent=2))
