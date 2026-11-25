import sys,os,collections; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import cpt
T=collections.defaultdict(lambda: collections.defaultdict(collections.Counter)); N=collections.Counter()
nf=0
for f,s,R in cpt.files():
  nf+=1
  if nf>250: break
  typ={i:t for (i,t,fl,b) in R}
  for (i,t,fl,b) in R:
    N[t]+=1; k=0
    for j in range(7):
      if k>=len(b): break
      try: v,k=cpt.sint(b,k)
      except Exception: break
      T[t][j][typ.get(i-v,'-') if v!=0 else 'zero']+=1
for t in sorted(N, key=lambda x:-N[x])[:30]:
  print(t,N[t],' | '.join('%d:%s'%(j,','.join('%s=%d'%(a,b) for a,b in T[t][j].most_common(3))) for j in range(6)))
