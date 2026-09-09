from pathlib import Path
import urllib.request,urllib.parse,json,time
P=Path(__file__).resolve().parents[1]
queries={
'human_obstruction':'Homo sapiens[Organism] AND (kidney OR renal) AND (obstruction OR obstructive OR hydronephrosis OR ureteropelvic) AND gse[Entry Type]',
'human_distal_injury':'Homo sapiens[Organism] AND (kidney) AND (single nucleus OR single cell) AND (injury OR fibrosis) AND gse[Entry Type]'
}
rows=[]
for name,q in queries.items():
 base='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
 url=base+'esearch.fcgi?'+urllib.parse.urlencode({'db':'gds','term':q,'retmode':'json','retmax':100})
 try:
  raw=urllib.request.urlopen(url,timeout=60).read();d=json.loads(raw);(P/'01_sources'/f'{name}_esearch.json').write_bytes(raw)
  ids=d['esearchresult']['idlist'];time.sleep(.4)
  su=base+'esummary.fcgi?'+urllib.parse.urlencode({'db':'gds','id':','.join(ids),'retmode':'json'})
  raw=urllib.request.urlopen(su,timeout=60).read();s=json.loads(raw);(P/'01_sources'/f'{name}_esummary.json').write_bytes(raw)
  for k in s['result']['uids']:
   v=s['result'][k];rows.append({'search':name,'accession':v.get('accession'),'title':v.get('title'),'summary':v.get('summary'),'samples':v.get('n_samples')})
  print(name,'count',d['esearchresult']['count'],flush=True)
 except Exception as e:print(name,str(e),flush=True)
 time.sleep(.4)
(P/'01_sources/native_dataset_candidates.json').write_text(json.dumps(rows,ensure_ascii=False,indent=2))
for row in rows: print(row['accession'],row['samples'],row['title'],flush=True)

