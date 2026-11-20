#!/usr/bin/env python3
"""Bounded native-width handoff and completion-accounting checks.

No census is started. Two near-terminal positive inputs exercise both native
backends, with production and artificially small lcm ceilings. Each native
run and Python continuation has a five-second limit. The compiler has a
30-second limit. C++ undefined-behavior checks are enabled. These finite
fixtures do not constitute a proof of every native arithmetic operation.
"""
import argparse
from contextlib import contextmanager
import gzip
import hashlib
import importlib.util
import json
from math import isqrt, lcm
from pathlib import Path
import re
import signal
import subprocess
import sys
import tempfile
from zipfile import ZipFile


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


@contextmanager
def deadline(seconds):
    def stop(*_):
        raise TimeoutError('bounded exact continuation exceeded its time limit')
    previous = signal.signal(signal.SIGALRM, stop)
    signal.alarm(seconds)
    try:
        yield
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, previous)


def positive_input(archive, rank, task_index):
    base = f'ranks/{rank:02d}/'
    tasks = json.loads(archive.read(base+'manifest.json'))['tasks']
    task = tasks[task_index]
    rows = [json.loads(line)[-1] for line in
            gzip.decompress(archive.read(base+'raw-census.jsonl.gz')).splitlines()
            if json.loads(line)[1] == task_index]
    assert rows
    # At rank 41 the selected task's row is known to survive the stronger
    # partner-capacity prefilter; the main terminal conditions still exclude it.
    original = min(rows)
    M = task['B']*original[-1]
    prefix = sorted(M//d for d in original[1:])[:-2]
    return dict(task, prefix=prefix), original


def main():
    if not __debug__:
        raise SystemExit('Run without -O or -OO.')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write', action='store_true')
    args = parser.parse_args()
    here = Path(__file__).resolve().parent
    root = here.parents[1]
    census = root/'census'
    sys.path[:0] = [str(census), str(census/'engine_src'), str(census/'engine')]
    import square
    import nonsquare
    from terminal_arithmetic import obstruction
    archive_path = census/'reference-census.zip'
    with ZipFile(archive_path) as archive:
        fixtures = [('square',39,*positive_input(archive,39,3)),
                    ('nonsquare',41,*positive_input(archive,41,395))]
    records = []
    with tempfile.TemporaryDirectory(prefix='rank-native-boundary-') as temporary:
        build = Path(temporary)
        for kind, rank, task, original in fixtures:
            adapter = square if kind == 'square' else nonsquare
            source = census/'engine_src'/('square44.cpp' if kind == 'square' else 'nonsquare.cpp')
            base_options = ['g++','-O1','-std=c++17','-fsanitize=undefined',
                            '-fno-sanitize-recover=all','-I',str(census/'engine')]
            observations = []
            rowsets = []
            for label, define in [('production', None),
                                  ('forced-handoff','SQ_LCM_BITS=3' if kind == 'square' else 'NS_LCM_BITS=1')]:
                executable = build/(kind+'-'+label)
                command = [*base_options]
                if define:
                    command.append('-D'+define)
                subprocess.run([*command,str(source),'-o',str(executable)],
                               check=True, capture_output=True, text=True, timeout=30)
                if kind == 'square':
                    tail = ['--rank',str(rank),'--prefix',','.join(map(str,task['prefix']))]
                else:
                    tail = [str(rank),str(task['q']),str(task['B']),str(task['m']),
                            '--prefix',','.join(map(str,task['prefix']))]
                run = subprocess.run([str(executable),*tail],check=True,
                                     capture_output=True,text=True,timeout=5)
                assert not run.stdout or run.stdout.endswith('\n')
                items = [json.loads(line) for line in run.stdout.splitlines()]
                native = [x for x in items if isinstance(x,list)]
                handoffs = [x for x in items if isinstance(x,dict)]
                assert len(native)+len(handoffs) == len(items)
                exact_rows = []
                for item in handoffs:
                    assert item['exact'][:len(task['prefix'])] == task['prefix']
                    with deadline(5):
                        exact_rows.extend(adapter.exact_rows(item,task,rank))
                values = {tuple(row) for row in native+exact_rows}
                assert len(values) == len(native)+len(exact_rows), 'duplicate completed row'
                assert tuple(original) in values, 'positive fixture lost at the width transition'
                trailer = (rf'DONE nodes=(\d+) profiles={len(native)} exact={len(handoffs)}'
                           if kind == 'square' else
                           rf'rank {rank}, q {task["q"]}, B {task["B"]}, multiplicity {task["m"]}, nodes (\d+), candidates {len(native)}, exact {len(handoffs)}')
                match = re.fullmatch(trailer,run.stderr.strip())
                assert match, 'completion trailer or sanitizer diagnostic failed'
                assert all(obstruction(list(row),task) is not None for row in values)
                if define:
                    assert handoffs and not native
                else:
                    assert native and not handoffs
                observations.append(dict(mode=label, native_rows=len(native),
                    handoffs=len(handoffs), exact_rows=len(exact_rows),
                    exact_generators_exhausted=True, nodes=int(match.group(1)),
                    completion=run.stderr.strip()))
                rowsets.append(values)
            assert rowsets[0] == rowsets[1]
            # Exercise the actual native integer-square-root routine at
            # large 64/128-bit boundaries, independently of the small fixture.
            vectors = {0,1,2,3,4,(1<<127)-1}
            for exponent in (24,48,58,61,63):
                for delta in (-1,0,1):
                    value = ((1<<exponent)+delta)**2
                    vectors.update((value-1,value,value+1))
            vectors = sorted(vectors)
            wrapper = build/(kind+'-roots.cpp')
            typename = 'Wide' if kind=='square' else 'W'
            function = 'square_root' if kind=='square' else 'root'
            wrapper.write_text(
                '#define main rank_enumerator_main\n#include "'+str(source)+'"\n#undef main\n'
                'int main(){std::string t; while(std::cin>>t){'+typename+' n=0;'
                'for(char c:t)n=n*10+unsigned(c-\'0\');'
                'std::cout<<'+function+'(n)<<"\\n";}}\n')
            executable = build/(kind+'-roots')
            subprocess.run([*base_options,str(wrapper),'-o',str(executable)],
                           check=True,capture_output=True,text=True,timeout=30)
            result = subprocess.run([str(executable)],
                input=''.join(str(value)+'\n' for value in vectors),
                check=True,capture_output=True,text=True,timeout=5)
            assert not result.stderr
            assert list(map(int,result.stdout.splitlines())) == [isqrt(value) for value in vectors]
            records.append(dict(backend=kind, rank=rank, task=task,
                original_row=original, prefix_lcm=lcm(*task['prefix']),
                modes=observations, completed_rowsets_equal=True,
                scalar_root_checks=len(vectors),
                scalar_input_sha256=hashlib.sha256(json.dumps(vectors).encode()).hexdigest()))
    result = dict(status='passed',
        scope='Two bounded positive native/exact fixtures; no partition generation or rank census.',
        sanitizer='undefined behavior, nonrecovering',
        reference_archive_sha256=sha(archive_path),
        source_sha256={str(path.relative_to(root)):sha(path) for path in (
            Path(__file__),census/'engine_src/square44.cpp',census/'engine_src/nonsquare.cpp',
            census/'engine_src/square.py',census/'engine_src/nonsquare.py',
            census/'engine_src/exact_square_legacy.py',census/'engine_src/exact_square44.py',
            census/'engine/nonsquare_exact.py',census/'engine/denominators.py')},
        fixtures=records)
    text = json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.write:
        (here/'certificate.json').write_text(text)
    else:
        assert (here/'certificate.json').read_text() == text, 'native boundary receipt differs'
    print('Both native/exact boundary fixtures passed, with identical completed row sets;')
    print('native square roots agree with unbounded integers on all large-width boundary inputs.')


if __name__ == '__main__':
    main()
