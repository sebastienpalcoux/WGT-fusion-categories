#!/usr/bin/env python3
"""Run compact bounded checks. This is not a full certificate replay."""
import argparse,ast,hashlib,json,platform,shutil,subprocess,sys,time,zipfile
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def require(condition,message):
 if not condition:raise ValueError(message)

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def integrity():
 names=set()
 for line in (ROOT/'SHA256SUMS.txt').read_text().splitlines():
  expected,name=line.split('  ',1);p=ROOT/name
  require(name not in names and not Path(name).is_absolute() and '..' not in Path(name).parts
   and p.is_file() and not p.is_symlink(),'Unsafe or duplicate manifest path: '+name)
  require(sha(p)==expected,'Hash differs: '+name)
  require(p.stat().st_size<10_000_000,'File exceeds 10 MB: '+name);names.add(name)
 actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file()
  and '__pycache__' not in p.parts and p.suffix!='.pyc'}
 require(actual==names|{'SHA256SUMS.txt'},'Unmanifested or missing file')
 index=json.loads((ROOT/'historical_sources.json').read_text())['files']
 with zipfile.ZipFile(ROOT/'historical_sources.zip') as z:
  require(set(z.namelist())==set(index),'Source archive entries differ')
  for name,entry in index.items():
   data=z.read(name)
   require(len(data)==entry['bytes'] and hashlib.sha256(data).hexdigest()==entry['sha256'],'Historical source differs: '+name)
   ast.parse(data.decode(),filename=name)
 return dict(status='PASS_FILE_INTEGRITY',files=len(names),historical_sources=len(index),numerical_replay_of_historical_sources=False)

def main():
 require(__debug__,'Do not use Python -O or PYTHONOPTIMIZE.')
 parser=argparse.ArgumentParser(description=__doc__)
 parser.add_argument('group',nargs='?',choices=('quick','standard-library','structural','integrity'),default='quick')
 parser.add_argument('--output',type=Path,default=ROOT.parent/'verification_output');args=parser.parse_args()
 output=args.output.resolve();require(not output.is_relative_to(ROOT),'Place generated output outside verification/.')
 output.mkdir(parents=True,exist_ok=True);start=time.monotonic()
 report=dict(integrity=integrity(),checks=[]);versions={'python':platform.python_version()}
 def run(name,argv,result_file=None,expected_file=None):
  process=subprocess.run([sys.executable,'-B',str(ROOT/'checks'/name),*argv],cwd=ROOT,text=True,capture_output=True)
  (output/(name+'.log')).write_text(process.stdout+process.stderr)
  require(process.returncode==0,'Check failed: '+name+'; see its log.')
  result=json.loads((output/result_file).read_text()) if result_file else json.loads(process.stdout)
  if expected_file:
   expected=json.loads((ROOT/'reference'/expected_file).read_text());require(result==expected,'Result differs from reference: '+name)
  report['checks'].append(dict(script=name,status='PASS',result=result,compared_with_reference=bool(expected_file)))
  print('PASS:',name,flush=True)
 if args.group not in ('integrity','structural'):
  run('check_certificate_scope.py',['--output',str(output/'certificate_scope.json')],'certificate_scope.json')
  run('first_odd_family_certificate.py',['--output',str(output/'odd_frontier.json')],'odd_frontier.json')
  run('verify_two_blocks.py',['--output',str(output/'two_blocks.json')],'two_blocks.json')
  run('verify_length10.py',['--output',str(output/'length10.json')],'length10.json')
  run('attachment_test_identities.py',[],expected_file='attachment_identity_results.json')
  run('check_ch_pair_independent.py',['--explorer',str(ROOT/'controls'),'--output',str(output/'coefficient_controls.json')],
   'coefficient_controls.json','ch_pair_independent_results.json')
  if args.group=='quick':
   import numpy,sympy
   versions.update(numpy=numpy.__version__,sympy=sympy.__version__)
   run('verify_p6q6r4_review.py',['--output',str(output/'independent_arithmetic.json')],'independent_arithmetic.json','p6q6r4_independent_checks.json')
   run('verify_dyadic_compression.py',['--output',str(output/'dyadic.json')],'dyadic.json','p6q6r4_dyadic_checks.json')
 if args.group=='structural':
  runtime=output/'structural_runtime';runtime.mkdir(exist_ok=True)
  for source in (ROOT/'checks/structural').glob('*.py'):shutil.copy2(source,runtime/source.name)
  for name in ['orthogonal_six_catalogue.py','audit_orthogonal_six.py','cubic_four_e3.py','audit_cubic_four_e3.py','cubic_four_e4.py','audit_cubic_four_e4.py','symplectic_profiles.py','verify_symplectic_profiles.py']:
   process=subprocess.run([sys.executable,'-B',str(runtime/name)],cwd=runtime,text=True,capture_output=True)
   (output/(name+'.log')).write_text(process.stdout+process.stderr)
   require(process.returncode==0,'Structural check failed: '+name)
   report['checks'].append(dict(script='structural/'+name,status='PASS'))
   print('PASS:',name,flush=True)
  for name in ['orthogonal_six_independent_audit.json','cubic_four_e3_independent_audit.json','cubic_four_e4_independent_audit.json','symplectic_profile_verification.json']:
   result=json.loads((runtime/name).read_text());expected=json.loads((ROOT/'reference'/name).read_text())
   # Timings are diagnostic; all mathematical data and source hashes must match.
   result.pop('seconds',None);expected.pop('seconds',None)
   require(result==expected,'Structural audit differs: '+name)
 report.update(status='PASS_COMPACT_BOUNDED_SCOPE',group=args.group,versions=versions,seconds=time.monotonic()-start,
  global_large_exclusions_replayed=False,hypothesis_V_verified=False,formal_proof_assistant_verification=False)
 (output/'verification_summary.json').write_text(json.dumps(report,indent=2,sort_keys=True)+'\n')
 print('PASS: compact bounded scope; global certificate coverage remains unverified.')
 print('Report:',output/'verification_summary.json')
if __name__=='__main__':main()
