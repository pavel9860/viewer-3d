import sys,struct
def be(b,i): return struct.unpack_from('>I',b,i)[0]
def streams(fn):
    b=open(fn,'rb').read()
    off,ln=be(b,8),be(b,12);c=b[off:off+ln];p=0;out={}
    while p+4<=len(c):
        nl=be(c,p)
        if 4<=nl<=1024 and p+4+nl+26<=len(c) and nl%2==0 and c[p+4+nl-2:p+4+nl]==b'\0\0':
            try: name=c[p+4:p+4+nl-2].decode('utf-16-be')
            except: name=''
            r=p+4+nl; rl=be(c,r)
            if name and name.isprintable() and rl>=26 and r+rl<=len(c) and c[r+8:r+10] in(b'\0\x01',b'\0\x02'):
                if c[r+9]==2:
                    size=be(c,r+14)
                    for n in range(1,4000):
                        t=r+rl-4-20*n
                        if t<r+26: break
                        if be(c,t)!=n: continue
                        es=[struct.unpack_from('>5I',c,t+4+20*k) for k in range(n)]
                        if sum(e[1] for e in es)==size:
                            es.sort(key=lambda e:e[3]); out[name]=b''.join(b[e[0]:e[0]+e[1]] for e in es); break
                p=r+rl; continue
        p+=1
    return out
