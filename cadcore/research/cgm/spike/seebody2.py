import sys,os,struct,collections
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import classes
from concurrent.futures import ThreadPoolExecutor
def lists(s,F):
  k=0; out=[]
  while k<len(s):
    if k+4>len(s): return None
    n=struct.unpack_from(F+'I',s,k)[0]; k+=4
    if n>100000 or k+4*n>len(s): return None
    out.append(list(struct.unpack_from(F+'%dI'%n,s,k))); k+=4*n
  return out
def one(p):
  try: cs=classes.classes(p)
  except Exception: return p,None
  sb=[s for c,n,s in cs if c=='CATSeeBodyCont' and s is not None]
  return p,(sb[0] if sb else None)
if __name__=='__main__':
  parts=[p for p in open('parts.txt').read().split('\n') if p]
  c=collections.Counter(); ex={}
  with ThreadPoolExecutor(4) as e:
    for p,s in e.map(one,parts):
      if s is None: c['none']+=1; continue
      L=lists(s,'<') or lists(s,'>')
      key='ok %d lists'%len(L) if L is not None else 'fail'
      c[key]+=1; ex.setdefault(key,(p.split('/')[-1],s[:60].hex(' ')))
  for k,v in c.most_common(): print(v,k,ex.get(k))
