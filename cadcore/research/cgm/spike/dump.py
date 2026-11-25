import sys,os,hashlib
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
sys.argv=[sys.argv[0]]+sys.argv[1:]
import cgmlib
out=sys.argv[2]
for p in open(sys.argv[1]).read().split('\n'):
  if not p: continue
  try: s=cgmlib.cgm_stream(p)
  except Exception: continue
  if s: open(os.path.join(out,hashlib.md5(p.encode()).hexdigest()[:10]+'.cgm'),'wb').write(s); print(hashlib.md5(p.encode()).hexdigest()[:10],len(s),p,sep='\t')
