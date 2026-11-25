import sys,os,struct,collections
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)));import store_streams as data
def cgm_stream(p):
  st=data.streams(p); d=st.get('Data',b'')
  if len(d)<16: return None
  big=d[8:10]==b'\0\1'; F='>' if big else '<'; a=0
  while a+16<=len(d):
    end=a+8*struct.unpack_from(F+'I',d,a+4)[0]; q=a+16; cls=g=None
    if end<=a: break
    while q+8<=end:
      tag,fl,l=struct.unpack_from(F+'HHI',d,q)
      if l==0: break
      if tag==2: cls=d[q+8:q+8*l].split(b'\0')[0]
      if tag==3: g=d[q+16:q+32]
      q+=8*l
    if cls==b'CGMGeom' and g:
      for k,v in st.items():
        if k.startswith('Cont_') and k.endswith('_'+g[4:8].hex().lstrip('0')) : return v
      for k,v in st.items():
        if k.startswith('Cont_') and g[0:4].hex().lstrip('0') in k and g[5:8].hex().lstrip('0') in k: return v
    a=end
  return None
c=collections.Counter(); types=collections.Counter(); heads=collections.Counter(); ex={}
for p in open(sys.argv[1]).read().split('\n'):
  if not p: continue
  try: s=cgm_stream(p)
  except Exception as e: c['err']+=1; continue
  if s is None: c['no cgm']+=1; continue
  k=s.find(b'\xe5\x0d')
  heads[s[:k].hex()[:12]]+=1
  ok=True; n=0
  while k<len(s):
    if s[k:k+2]!=b'\xe5\x0d': ok=False; break
    t=s[k+2:k+5].hex(); l=struct.unpack_from('<I',s,k+5)[0]
    types[t]+=1; ex.setdefault(t,p); n+=1
    k+= 13 if n==1 else 9+l
  c['chain ok' if ok and k==len(s) else 'chain broken']+=1
  if not(ok and k==len(s)) and c['chain broken']<4: print('BROKEN',p,k,len(s),s[k:k+16].hex())
print(c); print(heads.most_common(5)); print(len(types),'types'); print(types.most_common(40))
