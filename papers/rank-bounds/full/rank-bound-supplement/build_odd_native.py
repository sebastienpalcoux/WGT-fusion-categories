#!/usr/bin/env python3
"""Build the exact odd Gauss continuation (C++17 and GMP development files)."""
from pathlib import Path
import os,shlex,shutil,subprocess
root=Path(__file__).resolve().parent
compiler=shlex.split(os.environ.get('CXX','c++'))
flags=['-lgmpxx','-lgmp']
if shutil.which('pkg-config'):
    found=subprocess.run(['pkg-config','--cflags','--libs','gmpxx'],capture_output=True,text=True)
    if found.returncode==0:flags=shlex.split(found.stdout)
command=compiler+['-O3','-std=c++17',str(root/'arithmetic/gauss_native.cpp'),*flags,'-o',str(root/'gauss_native')]
subprocess.run(command,check=True)
print('Built the exact odd Gauss continuation.')
