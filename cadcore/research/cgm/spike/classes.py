import sys,os,struct
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import store_streams as data
def classes(p):
  st=data.streams(p); d=st.get('Data',b'')
  big=d[8:10]==b'\0\1'; F='>' if big else '<'; a=0; out=[]
  while a+16<=len(d):
    end=a+8*struct.unpack_from(F+'I',d,a+4)[0]; q=a+16; cls=g=None
    if end<=a: break
    while q+8<=end:
      tag,fl,l=struct.unpack_from(F+'HHI',d,q)
      if l==0: break
      if tag==2: cls=d[q+8:q+8*l].split(b'\0')[0].decode('latin1')
      if tag==3: g=d[q+8:q+8*l]
      q+=8*l
    if g:
      name='Cont'+''.join('_%x'%struct.unpack_from('>I',g,4*i)[0] for i in range(4))
      out.append((cls,name,st.get(name)))
    a=end
  return out
if __name__=='__main__':
  for cls,name,s in classes(sys.argv[1]): print(cls,name,len(s) if s is not None else None,(s or b'')[:32].hex(' '))
