import sys,os,glob,struct,collections
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import old
def hdr(m):
  hi,lo=m>>4,m&15
  lw={0xb:1,0xa:4}.get(hi)
  if lw is None or lo<lw+1: return None
  return lw,lo-lw
def chain(s):
  k=s.find(b'\xa8\x06\x7a')
  if k<0: return None,None,-1
  lvl,size=struct.unpack_from('<II',s,k+3); k+=11; st=k; recs=[]
  while k<st+size:
    h=hdr(s[k])
    if not h: return recs,(st,size,lvl),k
    lw,iw=h
    ln=int.from_bytes(s[k+3:k+3+lw],'little'); idf=s[k+3+lw:k+3+lw+iw]
    recs.append((k,s[k],s[k+1],s[k+2],ln,idf,k+3+lw+iw)); k+=3+lw+iw+ln
  return recs,(st,size,lvl),k
if __name__=='__main__':
  c=collections.Counter(); ex={}
  for f in sorted(glob.glob('cgm2/*.cgm')):
    s=open(f,'rb').read()
    if old.old_records(s): continue
    r,d,k=chain(s)
    if r is None: c['nodir']+=1; ex.setdefault('nodir',[]).append(f); continue
    if k==d[0]+d[1]: c['ok']+=1; continue
    key=s[k:k+1].hex(); c[key]+=1; ex.setdefault(key,(f,k,s[k-6:k+16].hex(' ')))
  print(c); [print(k,v if k!='nodir' else v[:5]) for k,v in ex.items()]
