import sys,os,collections; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import old
s=open(sys.argv[1],'rb').read(); r,e=old.old_records(s)
c=collections.Counter(t for (_,_,t,_,_,_) in r); print(sorted(c.items()))
want=set(map(int,sys.argv[2].split(','))) if len(sys.argv)>2 else None
for (o,fl,t,l,i,b) in r:
  if want is None or t in want: print(i,t,fl,l,b.hex(' '))
