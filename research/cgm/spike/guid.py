import sys,os,struct,collections
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import store_streams as data
def info(p):
  st=data.streams(p); d=st.get('Data',b'')
  big=d[8:10]==b'\0\1'; F='>' if big else '<'; a=0; out=[]
  while a+16<=len(d):
    end=a+8*struct.unpack_from(F+'I',d,a+4)[0]; q=a+16; cls=g=None
    if end<=a: break
    while q+8<=end:
      tag,fl,l=struct.unpack_from(F+'HHI',d,q)
      if l==0: break
      if tag==2: cls=d[q+8:q+8*l].split(b'\0')[0]
      if tag==3: g=d[q+8:q+8*l]
      q+=8*l
    if cls==b'CGMGeom': out.append(g)
    a=end
  return out,[k for k in st if k.startswith('Cont_')],big
c=collections.Counter()
for p in open('parts.txt').read().split('\n')[:1565:40]:
  if not p: continue
  try: g,names,big=info(p)
  except Exception as e: continue
  print(big,[x.hex() for x in g],names)
