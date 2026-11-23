import sys,os; sys.path.insert(0,os.path.dirname(os.path.abspath(__file__))); import old,new3
def dialect(s):
  a=s.find(b'\xa8\x06\x7a'); e=s.find(b'\xe5\x0d')
  if a>=0 and (e<0 or a<e):
    r,d,k=new3.chain(s); return 'compact' if r is not None and k==d[0]+d[1] else 'compact?'
  if old.old_records(s): return 'varint-be' if s[0]==0 else 'varint'
  if b'\x00\xb3\x00' in s[:80]: return 'fixed'
  return '?'
