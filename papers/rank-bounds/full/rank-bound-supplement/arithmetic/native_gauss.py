"""Exact native continuation, with independently maintained Python terminal checks."""
from pathlib import Path
import os,subprocess
import ordinary_tail

class NativeUnresolved(RuntimeError):
    pass


def applicable(c):
    return (os.environ.get('ODD_NATIVE')=='1'
            and ((c['support']==(3,5,7,11) and c['marked']=={11:(1,22,20)})
                 or (c['support']==(3,5,7) and c['marked']=={}))
            and c['isolated']==(7,) and c['qcaps'][0] in (12,24,36)
            and (len(c['support'])==3 or c['neutral']==14)
            and c['z'] in (32,48,64))



def possible(c,charges):
    # A wide union weakens its count envelopes. Bounded batches retain
    # every exact state while still sharing denominator prefixes. Keeping
    # the original caps and neutral reservations only enlarges each batch.
    states=tuple(sorted(set(c['states'])))
    for offset in range(0,len(states),8):
        context=dict(c)
        context['states']=states[offset:offset+8]
        if _run(context,charges):
            return True
    return False


def _run(c,charges):
    executable=Path(__file__).resolve().parents[1]/'gauss_native'
    if not executable.is_file():
        raise NativeUnresolved('native Gauss executable is missing')
    states=sorted(set(c['states']))
    header=[len(c['support']),c.get('wq',0),c['omit3'],c['w5'],c['nines'],c['z'],c['big_z'],c['neutral'],c['qcaps'][0],c['slots'],c['root'],c['scale_root'],
            c['initial_mass'],c['initial_charge'],len(states)]
    rows=[' '.join(map(str,header))]
    rows.extend(' '.join(map(str,state)) for state in states)
    rows.append(str(int(c['caps'] is not None)))
    if c['caps'] is not None:
        for side in c['caps']:
            records=[(a,x) for a,totals in side.items() for x in totals]
            rows.append(str(len(records)))
            rows.extend(f'{a} {x}' for a,x in records)
    rows.append(str(len(charges)))
    for charge,imbalances in charges:
        rows.append(' '.join(map(str,[charge,len(imbalances),*imbalances])))
    process=subprocess.Popen([str(executable),os.environ.get('ODD_NATIVE_NODES','0')],
                             stdin=subprocess.PIPE,stdout=subprocess.PIPE,
                             text=True,bufsize=1)
    try:
        process.stdin.write('\n'.join(rows)+'\n');process.stdin.flush()
        for line in process.stdout:
            fields=line.split()
            if not fields:
                raise NativeUnresolved('empty native result')
            if fields[0]=='T':
                if len(fields)!=9:
                    raise NativeUnresolved('malformed terminal handoff')
                M,L,k,mass,scale,first,o3,o5=map(int,fields[1:])
                answer=ordinary_tail.possible(M,c['support'],c['marked'],L,c['nines'],
                    k,c['neutral'],o3,o5,mass,scale,c['witnesses'],(first,),c['qexact'])
                process.stdin.write(str(int(answer))+'\n');process.stdin.flush()
            elif fields[0]=='DONE':
                if len(fields)!=4 or fields[1] not in ('0','1'):
                    raise NativeUnresolved('malformed native completion')
                if process.wait()!=0:
                    raise NativeUnresolved('native completion exited unsuccessfully')
                if os.environ.get('ODD_NATIVE_TRACE')=='1':
                    print('NATIVE',c['nines'],c['shape'],len(states),
                          'nodes',fields[2],'terminals',fields[3],flush=True)
                return fields[1]=='1'
            else:
                raise NativeUnresolved(line.strip())
        raise NativeUnresolved('native process ended without completion')
    finally:
        if process.poll() is None:
            process.kill()
        process.wait()
        process.stdin.close();process.stdout.close()
