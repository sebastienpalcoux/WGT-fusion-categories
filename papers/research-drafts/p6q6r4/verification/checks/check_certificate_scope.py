#!/usr/bin/env python3
"""Synthetic controls for row-set enlargement using the hardened exact checker."""
import argparse, hashlib, json
from pathlib import Path
from exact_certificate_checks import verify_certificate

def digest(value):
 return hashlib.sha256(json.dumps(value,separators=(',',':')).encode()).hexdigest()

def audit():
 if not __debug__: raise RuntimeError('Assertions must be enabled')
 old=[(1,0,0,0,0,0)]; target=(1,0,1,0,0,0); enlarged=old+[target]
 separator=[0,0,-1,0,0,0]
 certificates=[dict(kind='integer_separator',coefficients=separator,minimum_row_value=0,target_value=-1),
  dict(kind='integer_box_separator',coefficients=separator,lower_bound=0,target_value=-1),
  dict(kind='modular_separator',coefficients=[0,0,1,0,0,0],modulus=2,target_residue=1)]
 outcomes=[]
 for bad_rows,bad_target,bad_cert in [([], (0,0,0,0,0,0), {'kind':'empty'}),
   ([(0,0,0,0,0,0)],target,{'kind':'bounded_empty'}),
   ([(-1,0,0,0,0,0)],target,{'kind':'bounded_empty'})]:
  rejected=False
  try:verify_certificate(bad_rows,bad_target,bad_cert)
  except AssertionError:rejected=True
  assert rejected,'Invalid budget domain or zero target was accepted'
 for cert in certificates:
  assert verify_certificate(old,target,cert)
  rejected=False
  try:verify_certificate(enlarged,target,cert)
  except AssertionError:rejected=True
  assert rejected,'Witness was incorrectly reused after enlargement'
  outcomes.append(dict(kind=cert['kind'],accepted_on_old_set=True,rejected_after_enlargement=True))
 larger=old+[(1,0,0,1,0,1)]
 assert verify_certificate(larger,target,certificates[0])
 assert verify_certificate(old,target,certificates[0])
 return dict(status='PASS_SYNTHETIC_SCOPE_CONTROLS',old_rows=old,added_row=target,target=target,
  old_rows_sha256=digest(old),enlarged_rows_sha256=digest(enlarged),cases=outcomes,
  exact_checker_sha256=hashlib.sha256(Path(__file__).with_name('exact_certificate_checks.py').read_bytes()).hexdigest(),
  valid_superset_implies_valid_subset=True,global_field_catalogue_coverage_verified=False,
  scope='Synthetic regression only. No large certificate is cleared by this test.')

if __name__=='__main__':
 parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--output',type=Path,required=True)
 args=parser.parse_args();args.output.parent.mkdir(parents=True,exist_ok=True)
 args.output.write_text(json.dumps(audit(),indent=2)+'\n')
 print('PASS: linear, bounded and modular witnesses are rechecked on enlarged rows.')
