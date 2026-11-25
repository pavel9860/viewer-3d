import sys,os,glob,struct
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import chain
def old_records(s):
  k=s.find(b'\xe5\x0d')
  if k<0: return None
  r,e=chain.chain(s)
  size=struct.unpack_from('<I',s,k+9)[0]
  if e-(k+13)!=size: return None
  return [(o,fl,t,l,i,s[o+13:o+13+l]) for (o,fl,t,l,i) in r[1:]], e
def files():
  for f in sorted(glob.glob(os.path.join(os.path.dirname(os.path.abspath(__file__)),'cgm2/*.cgm'))):
    s=open(f,'rb').read(); x=old_records(s)
    if x: yield f,s,x[0],x[1]
import re
def version(s):
  m=re.search(rb'[ -~]{6,}',s[:64]); return m.group().decode() if m else '?'
def hdr(s):
  k=s.find(b'\xe5\x0d'); return s[:k]
def levels(s,e):
  # trailer: refs list, then f32?/u32, 6 f64 tolerances (48 bytes), then two varints
  import vi
  c,k,_=vi.rd(s,e)
  for _ in range(c): x,k,_=vi.rd(s,k)
  # skip u32 then 6 doubles (variable: count leading 80s)
  return s[k:k+80]
