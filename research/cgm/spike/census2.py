import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import sys,glob,collections,chain
c=collections.Counter(); brk=collections.Counter(); ex={}
for f in sorted(glob.glob('cgm2/*.cgm')):
  s=open(f,'rb').read(); r,k=chain.chain(s)
  if k==len(s): c['ok']+=1
  else:
    c['broken']+=1; key=s[k:k+2].hex(); brk[key]+=1; ex.setdefault(key,(f,k,len(s)))
print(c); print(brk.most_common(20)); 
for k,v in list(ex.items())[:10]: print(k,v)
