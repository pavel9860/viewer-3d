import sys,os,glob,struct,collections
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import vi,old
def dirpos(s):
  k=s.find(b'\xa8\x06')
  return k
def chain(s):
  k=dirpos(s)
  if k<0: return None
  # a8 06 <varint?> u32 count u32 size ?
  return k
if __name__=='__main__':
  c=collections.Counter()
  for f in sorted(glob.glob('cgm2/*.cgm')):
    s=open(f,'rb').read()
    if old.old_records(s): continue
    k=s.find(b'\xa8\x06'); 
    c[(k>=0)]+=1
    if k>=0 and c[True]<12: print(f, k, s[k-12:k+24].hex(' '))
  print(c)
