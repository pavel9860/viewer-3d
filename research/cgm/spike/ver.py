import sys,os,glob,collections,re; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import chain
c=collections.Counter()
for f in sorted(glob.glob('cgm2/*.cgm')):
  s=open(f,'rb').read(); m=re.search(rb'[ -~]{6,}',s[:64]); v=m.group().decode() if m else '?'
  v=re.sub(r'[:@].*','',v).strip()
  k=s.find(b'\xe5\x0d'); c[(v,s[:4].hex(), k>=0)]+=1
for k,n in sorted(c.items(),key=lambda x:-x[1]): print(n,k)
