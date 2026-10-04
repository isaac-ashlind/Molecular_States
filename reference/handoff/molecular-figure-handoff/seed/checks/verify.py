import json, math
from pathlib import Path
E=(tuple(range(5)),0)
def mul(a,b):return(tuple(a[0][b[0][i]] for i in range(5)),a[1]^b[1])
def perm(*cycles):
 p=list(range(5))
 for c in cycles:
  for i,j in zip(c,c[1:]+c[:1]):p[i-1]=j-1
 return tuple(p)
def close(gs):
 s={E};todo=[E]
 while todo:
  a=todo.pop()
  for b in gs:
   c=mul(a,b)
   if c not in s:s.add(c);todo.append(c)
 return frozenset(s)
b=(perm((2,3),(4,5)),1);t=(perm((1,2,3)),0);u=(perm((4,5)),0)
H=close([b]);G6=close([b,t]);G12=close([b,t,u]);B=close([b,t,u,(E[0],1)])
assert [len(g) for g in [H,G6,G12,B]]==[2,6,12,24]
assert mul(b,mul(t,b))==mul(t,t) and mul(b,t)!=mul(t,b)
assert close([b,t,mul(t,u)])==G12
assert len({frozenset(mul(g,h) for h in H) for g in G12})==6
assert 240//len(G12)==20
interval={H};todo=[H]
while todo:
 k=todo.pop()
 for g in B-k:
  z=close(list(k)+[g])
  if z not in interval:interval.add(z);todo.append(z)
assert len(interval)==10
assert len([(a,z) for a in interval for z in interval if a<z and not any(a<k<z for k in interval)])==17
def countcycles(p):
 seen=set();n=0
 for i in range(5):
  if i not in seen:
   n+=1;j=i
   while j not in seen:seen.add(j);j=p[j]
 return n
weights={}
for parity in [1,-1]:
 vals={'A1':0,'A2':0,'E':0}
 for g in G6:
  cy=countcycles(g[0]);factor=2**cy*(-1)**(5-cy)*parity**g[1]
  chars={'A1':1,'A2':(-1)**g[1],'E':0 if g[1] else(2 if g==E else -1)}
  for s in vals:vals[s]+=chars[s]*factor
 assert all(v%6==0 for v in vals.values())
 weights[str(parity)]={s:v//6 for s,v in vals.items()}
assert weights=={'1':{'A1':12,'A2':4,'E':8},'-1':{'A1':4,'A2':12,'E':8}}
for w in weights.values():assert w['A1']+w['A2']+2*w['E']==32
vs=[(1/math.sqrt(3),)*3,(2/math.sqrt(6),-1/math.sqrt(6),-1/math.sqrt(6)),(0,1/math.sqrt(2),-1/math.sqrt(2))]
for i,a in enumerate(vs):
 for j,z in enumerate(vs):assert abs(sum(x*y for x,y in zip(a,z))-(i==j))<1e-12
for i in range(1,401):
 x=i/200;c=math.exp(-1/(8*x*x));assert abs((1+2*c)/3+2*(1-c)/3-1)<1e-12
report={'passed':['group orders','generator extension','coset counts','bond interval and covers','spin weights and dimension sum','coefficient orthogonality','Gaussian sum rule'],'spin_weights':weights,'not_verified':['reference geometry','global normal topology','cell partitions','covering domain','seam/Wigner conventions','other molecule tables','LaTeX compilation','renders']}
print(json.dumps(report,indent=2))
