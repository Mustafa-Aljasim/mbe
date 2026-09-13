"""Compare fresh acceptance evidence with the pre-edit post-5C baseline."""
from pathlib import Path
import json
import runpy

ROOT=Path(__file__).resolve().parents[1]

if __name__=='__main__':
    for phase in ('5c','5b'):
        result=runpy.run_path(str(ROOT/f'examples/run_phase{phase}_acceptance.py'))['evidence']()
        normalized=json.loads(json.dumps(result,default=str,allow_nan=False))
        baseline=json.loads((ROOT/f'docs/post5c_baseline/phase_{phase}_numerical_results.json').read_text())
        (ROOT/f'docs/post5c_phase_{phase}_results.json').write_text(json.dumps(normalized,indent=2,allow_nan=False)+'\n',encoding='utf-8')
        assert normalized==baseline,f'Phase {phase} numerical output changed'
        print(f'Phase {phase}: EXACT BASELINE MATCH',flush=True)
    print('Phase 5B recursively verifies Phase 5A and prior acceptance outputs, including Phase 4C.',flush=True)
