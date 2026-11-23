import sys,os,collections; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import old
c=collections.Counter(); L=collections.defaultdict(collections.Counter); F=collections.defaultdict(collections.Counter); n=0
for f,s,r,e in old.files():
  n+=1
  for (o,fl,t,l,i,b) in r: c[t]+=1; L[t][l]+=1; F[t][fl]+=1
print(n,'files')
for t,v in c.most_common(): print(t,v,dict(F[t]),L[t].most_common(4))
