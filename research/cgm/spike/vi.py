def rd(b,k):
  t=b[k]
  if t>=0x80: return t&0x7f,k+1,0
  m=(t>>3)&0xf; v=0; k+=1
  for j in range(4):
    if m>>j&1: v|=b[k]<<(8*j); k+=1
  return v,k,t&7
