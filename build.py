#!/usr/bin/env python3
"""Clean-state build of the plates and the manuscript.

    python3 build.py            # data -> checks -> plates -> text-collision report -> manuscript
    python3 build.py --only 03  # one plate (data and checks still run; no manuscript)

Requires python3 with SymPy, pdflatex and bibtex with standalone, TikZ, PGFPlots and titlesec, and poppler (pdftoppm
for the 600 dpi plates, pdftotext for the text-collision report). GAP, when present, runs an independent cross-check.
Any failing step stops the build.
"""
import os, re, shutil, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / 'figures' / 'src'
BUILD = ROOT / 'build'
FIGS = BUILD / 'figures'

def run(cmd, cwd=ROOT, env=None, quiet=True):
    r = subprocess.run(cmd, cwd=cwd, env=env, capture_output=quiet, text=True)
    if r.returncode != 0:
        if quiet:
            sys.stderr.write(r.stdout[-4000:] + r.stderr[-2000:])
        raise SystemExit(f'failed: {" ".join(str(c) for c in cmd)}')
    return r

def deliver(pdfs):
    """Copy each plate PDF to figures/pdf and render it at 600 dpi to figures/png."""
    out, png = ROOT / 'figures' / 'pdf', ROOT / 'figures' / 'png'
    out.mkdir(exist_ok=True); png.mkdir(exist_ok=True)
    for pdf in pdfs:
        shutil.copy2(pdf, out / pdf.name)
        if shutil.which('pdftoppm'):
            run(['pdftoppm', '-r', '600', '-png', '-singlefile', str(pdf), str(png / pdf.stem)])
    r = subprocess.run([sys.executable, 'checks/collisions.py'] + [str(p) for p in pdfs],
                       cwd=ROOT, capture_output=True, text=True)
    print(r.stdout.strip())
    return r.stdout

def check_numbering(ms):
    """Figure k sits in Section k and is numbered k (k = 1..12); the guide and the closing plate are unnumbered."""
    src = (ROOT / 'docs' / 'manuscript.tex').read_text()
    placed = [re.findall(r'\\label\{fig:(\d+)\}', s) for s in src.split('\\subsection{')[1:]]
    numbers = dict(re.findall(r'\\newlabel\{fig:(\d+)\}\{\{(\d+)\}', (ms / 'manuscript.aux').read_text()))
    want = [str(k) for k in range(1, 13)]
    figures = re.findall(r'\\begin\{figure\}.*?\\end\{figure\}', src, re.S)
    if placed != [[k] for k in want] or numbers != {k: k for k in want} or sum(f.count('\\caption{') for f in figures) != 12:
        raise SystemExit(f'figure numbering mismatch: sections hold {placed}, numbers {numbers}')
    print('    figures 1-12 numbered as their sections; the guide and the closing plate unnumbered')

def main():
    only = None
    if '--only' in sys.argv:
        only = sys.argv[sys.argv.index('--only') + 1]
    print('[1/5] regenerate data'); run([sys.executable, 'compute/make_data.py'])
    print('[2/5] checks')
    for script in ('verify.py', 'verify_spin.py', 'verify_methane.py', 'verify_symbolic.py'):
        run([sys.executable, 'checks/' + script], quiet=False)
    if shutil.which('gap'):
        run(['gap', '-q', '-b', '--quitonbreak', 'checks/verify.g'], quiet=False)
    else:
        print('      GAP cross-check: gap not found, skipped')
    if not shutil.which('pdflatex'):
        raise SystemExit('pdflatex not installed')
    if FIGS.exists():
        shutil.rmtree(FIGS)
    FIGS.mkdir(parents=True)
    env = os.environ.copy()
    env['TEXINPUTS'] = os.pathsep.join([str(ROOT / 'figures' / 'shared') + '//',
                                        str(ROOT / 'figures' / 'data') + '//',
                                        env.get('TEXINPUTS', '')])
    sources = sorted(SRC.glob('fig*.tex'))
    missing = [f'fig{i:02d}' for i in range(14) if f'fig{i:02d}' not in [s.stem.split('-')[0] for s in sources]]
    if missing and not only:
        raise SystemExit('plates missing: ' + ', '.join(missing))
    print('[3/5] compile plates')
    for src in sources:
        if only and not src.stem.startswith('fig' + only):
            continue
        run(['pdflatex', '-interaction=nonstopmode', '-halt-on-error',
             '-output-directory=' + str(FIGS), str(src)], env=env)
        log = (FIGS / (src.stem + '.log')).read_text(errors='replace')
        over = re.findall(r'Overfull \\hbox', log)
        print(f'    {src.stem}.pdf', f'(overfull boxes: {len(over)})' if over else '')
    print('[4/5] plates to figures/pdf and, at 600 dpi, figures/png; text-collision report')
    report = deliver(sorted(FIGS.glob('fig*.pdf')))
    if only:
        return
    (BUILD / 'collisions.txt').write_text(report)
    print('[5/5] manuscript')
    ms = BUILD / 'ms'; ms.mkdir(exist_ok=True)
    tex = ['pdflatex', '-interaction=nonstopmode', '-halt-on-error', '-output-directory=' + str(ms), 'docs/manuscript.tex']
    run(tex)
    benv = os.environ.copy(); benv['BIBINPUTS'] = str(ROOT / 'docs') + '//' + os.pathsep
    run(['bibtex', 'manuscript'], cwd=ms, env=benv)
    run(tex); run(tex)
    log = (ms / 'manuscript.log').read_text(errors='replace')
    undefined = len(re.findall(r"(?:Reference|Citation) `[^']*' on page \d+ undefined", log))
    bad = len(re.findall(r'(?:Overfull|Underfull) \\[hv]box', log))
    pages = re.search(r'Output written on .*?\((\d+) pages', log)
    shutil.copy2(ms / 'manuscript.pdf', BUILD / 'manuscript.pdf')
    print(f'    build/manuscript.pdf: {pages.group(1) if pages else "?"} pages, {undefined} undefined references, {bad} over- or underfull boxes')
    if undefined:
        raise SystemExit('undefined references in the manuscript')
    check_numbering(ms)
    run([sys.executable, 'checks/palette.py'], quiet=False)   # the sixteen tones, on every plate
    print('done:', BUILD / 'manuscript.pdf')

if __name__ == '__main__':
    main()
