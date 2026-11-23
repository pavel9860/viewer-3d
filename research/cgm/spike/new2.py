import sys,os,glob,struct,collections
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import old
def chain(s,H):
  k=s.find(b'\xa8\x06\x7a')
  if k<0: return None,None,-1
  lvl,size=struct.unpack_from('<II',s,k+3); k+=11; st=k; recs=[]
  while k<st+size:
    m=s[k]
    if m not in H: return recs,(st,size),k
    hl,lenof,lenw=H[m]
    ln=int.from_bytes(s[k+lenof:k+lenof+lenw],'little')
    recs.append((k,m,s[k+1],s[k+2],ln)); k+=hl+ln
  return recs,(st,size),k
H={0xb5:(8,3,1),0xb2:(5,3,1)}
if __name__=='__main__':
  c=collections.Counter(); ex={}
  for f in sorted(glob.glob('cgm2/*.cgm')):
    s=open(f,'rb').read()
    if old.old_records(s): continue
    r,d,k=chain(s,H)
    if r is None: c['nodir']+=1; continue
    if k==d[0]+d[1]: c['ok']+=1; continue
    key=s[k:k+1].hex(); c[key]+=1; ex.setdefault(key,(f,k,s[k-6:k+16].hex(' ')))
  print(c); [print(k,v) for k,v in ex.items()]
