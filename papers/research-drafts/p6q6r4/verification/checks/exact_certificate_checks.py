"""Exact finite-moment separation with explicit input preconditions.

Every accepted certificate is evaluated on all supplied rows. The caller must
prove that these rows contain every admissible row; this function cannot infer
completeness. The first coordinates encode nonnegative atom-count/norm budgets.
"""
def validate(vectors,target):
 assert target and all(type(x) is int for x in target)
 assert any(target),'The zero target is represented by the empty combination'
 assert all(len(row)==len(target) and all(type(x) is int for x in row) for row in vectors)

def upper_bounds(vectors,target):
 validate(vectors,target)
 nblocks=next((n for n in range(1,20) if 3+3*n+n*(n-1)//2==len(target)),None)
 assert nblocks is not None
 positive=(0,1,2)+tuple(3+3*i+j for i in range(nblocks) for j in (0,2))
 assert all(target[j]>=0 for j in positive)
 upper=[]
 for row in vectors:
  assert all(row[j]>=0 for j in positive)
  active=[j for j in positive if row[j]>0]
  assert active,'Each row needs a positive nonnegative-budget coordinate'
  bound=min(target[j]//row[j] for j in active)
  assert bound>=0
  upper.append(bound)
 return upper

def dot(a,b):
 assert len(a)==len(b)
 return sum(x*y for x,y in zip(a,b))

def verify_certificate(vectors,target,certificate):
 validate(vectors,target)
 if certificate is None:return False
 kind=certificate['kind']
 if kind=='empty':
  assert not vectors;return True
 if kind=='bounded_empty':
  assert not any(upper_bounds(vectors,target));return True
 z=certificate['coefficients'];assert all(type(x)is int for x in z)
 goal=dot(z,target)
 if kind=='integer_separator':
  values=[dot(z,row) for row in vectors]
  assert values and min(values)>=0 and goal<0
  assert goal==certificate['target_value'] and min(values)==certificate['minimum_row_value']
 elif kind=='integer_box_separator':
  upper=upper_bounds(vectors,target)
  lower=sum(u*min(0,dot(z,row)) for row,u in zip(vectors,upper))
  assert goal<lower
  assert goal==certificate['target_value'] and lower==certificate['lower_bound']
 elif kind=='modular_separator':
  upper=upper_bounds(vectors,target);modulus=certificate['modulus']
  assert type(modulus) is int and modulus>=2
  assert all(dot(z,row)%modulus==0 for row,u in zip(vectors,upper) if u)
  assert goal%modulus==certificate['target_residue']!=0
 else:raise AssertionError(kind)
 return True
