import sys,os,glob,collections,struct; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import chain
c=collections.Counter(); ex={}
for f in sorted(glob.glob('cgm2/*.cgm')):
  s=open(f,'rb').read()
  if b'\xe5\x0d' not in s: continue
  r,k=chain.chain(s)
  d=r[0]; size=d[4] if False else struct.unpack_from('<I',s,d[0]+9)[0]
  start=d[0]+13
  key=('dirsize==recbytes' if k-start==size else 'dirsize!=') , s[k:k+1].hex()
  c[key]+=1; ex.setdefault(key,(f,k,len(s),k-start,size,len(r)))
for k,v in c.most_common(): print(v,k,ex[k])
