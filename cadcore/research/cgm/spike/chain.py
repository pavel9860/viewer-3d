import sys,struct
# old framing: e5 0d|0f  flag  u16 type  u32 len  u32 id  [i16 extra if 0f]  body[len]
def chain(s):
  k=s.find(b'\xe5\x0d'); n=0; recs=[]
  while k<len(s):
    if s[k]!=0xe5 or s[k+1] not in (0x0d,0x0f): return recs,k
    fl=s[k+2]; t=struct.unpack_from('<H',s,k+3)[0]; l=struct.unpack_from('<I',s,k+5)[0]
    if n==0: recs.append((k,fl,t,l,None)); k+=13; n+=1; continue
    rid=struct.unpack_from('<I',s,k+9)[0]; hl=13 if s[k+1]==0x0d else 15
    recs.append((k+hl-13,fl,t,l,rid)); n+=1
    k+=hl+l
  return recs,k
