from pathlib import Path
import os, shutil, subprocess
root=Path(__file__).resolve().parent
expected=[root/'figures/src'/f'{i:02d}.tex' for i in range(13)]
missing=[str(p.relative_to(root)) for p in expected if not p.exists()]
if missing: raise SystemExit('Figures not implemented: '+', '.join(missing))
if not shutil.which('pdflatex'):raise SystemExit('pdflatex not installed')
out=root/'build/figures';out.mkdir(parents=True,exist_ok=True)
env=os.environ.copy();env['TEXINPUTS']=str(root/'figures/shared')+'//:'+env.get('TEXINPUTS','')
for src in expected:
    subprocess.run(['pdflatex','-interaction=nonstopmode','-halt-on-error','-output-directory='+str(out),str(src)],cwd=root,env=env,check=True)
for _ in range(2):subprocess.run(['pdflatex','-interaction=nonstopmode','-halt-on-error','-output-directory='+str(root/'build'),'scaffold/outline.tex'],cwd=root,env=env,check=True)
