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
