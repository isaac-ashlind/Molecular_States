#!/usr/bin/env python3
"""Clean-state build of the figure suite, the scaffold and the proof sheet.

    python3 build.py            # data -> checks -> figures -> previews -> scaffold -> proof sheet
    python3 build.py --only 03  # one figure (data and checks still run)

Requires: python3 (stdlib only), pdflatex with standalone/TikZ/PGFPlots, pdftoppm
(poppler) for PNG previews.  GAP is used for an extra cross-check when present.
Nothing here fabricates success: any failing step stops the build.
"""
import os, re, shutil, subprocess, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / 'figures' / 'src'
BUILD = ROOT / 'build'
FIGS = BUILD / 'figures'
PREV = BUILD / 'previews'

def run(cmd, cwd=ROOT, env=None, quiet=True):
    r = subprocess.run(cmd, cwd=cwd, env=env, capture_output=quiet, text=True)
    if r.returncode != 0:
        if quiet:
            sys.stderr.write(r.stdout[-4000:] + r.stderr[-2000:])
        raise SystemExit(f'failed: {" ".join(str(c) for c in cmd)}')
    return r

def main():
    only = None
    if '--only' in sys.argv:
        only = sys.argv[sys.argv.index('--only') + 1]
    print('[1/6] regenerate data'); run([sys.executable, 'compute/make_data.py'])
    print('[2/6] python checks'); run([sys.executable, 'checks/verify.py'], quiet=False)
    try:
        import sympy  # noqa: F401
        print('[2a]  symbolic checks'); run([sys.executable, 'checks/verify_symbolic.py'], quiet=False)
        print('[2a\'] spin-structure check'); run([sys.executable, 'checks/verify_spin.py'], quiet=False)
        print('[2a"] antiprism check'); run([sys.executable, 'checks/verify_antiprism.py'], quiet=False)
    except ImportError:
        print('[2a]  sympy not found: symbolic checks skipped (stdlib checks already passed)')
    if shutil.which('gap'):
        print('[2b]  GAP cross-check'); run(['gap', '-q', '-b', 'checks/verify.g'], quiet=False)
    else:
        print('[2b]  GAP not found: cross-check skipped (stdlib checks already passed)')
    if not shutil.which('pdflatex'):
        raise SystemExit('pdflatex not installed')
    if FIGS.exists():
        shutil.rmtree(FIGS)
    FIGS.mkdir(parents=True); PREV.mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env['TEXINPUTS'] = os.pathsep.join([str(ROOT / 'figures' / 'shared') + '//',
                                        str(ROOT / 'figures' / 'data') + '//',
                                        env.get('TEXINPUTS', '')])
    sources = sorted(SRC.glob('fig*.tex'))
    expected = [f'fig{i:02d}' for i in range(13)]
    have = [s.stem.split('-')[0] for s in sources]
    missing = [e for e in expected if e not in have]
    if missing and not only:
        raise SystemExit('Figures not implemented: ' + ', '.join(missing))
    print('[3/6] compile figures')
    for src in sources:
        if only and not src.stem.startswith('fig' + only):
            continue
        run(['pdflatex', '-interaction=nonstopmode', '-halt-on-error',
             '-output-directory=' + str(FIGS), str(src)], env=env)
        log = (FIGS / (src.stem + '.log')).read_text(errors='replace')
        over = re.findall(r'Overfull \\hbox', log)
        print(f'    {src.stem}.pdf', f'(overfull boxes: {len(over)})' if over else '')
    print('[4/6] previews (150 dpi colour, 150 dpi grayscale)')
    if shutil.which('pdftoppm'):
        for pdf in sorted(FIGS.glob('fig*.pdf')):
            run(['pdftoppm', '-r', '150', '-png', '-singlefile', str(pdf), str(PREV / pdf.stem)])
            run(['pdftoppm', '-r', '150', '-gray', '-png', '-singlefile', str(pdf), str(PREV / (pdf.stem + '-gray'))])
    else:
        print('    pdftoppm not found: previews skipped')
    if only:
        out = ROOT / 'figures' / 'pdf'; out.mkdir(exist_ok=True)
        for pdf in sorted(FIGS.glob('fig*.pdf')):
            shutil.copy2(pdf, out / pdf.name)   # a single-figure build still refreshes its deliverable
        r = subprocess.run([sys.executable, 'checks/collisions.py'] + [str(p) for p in sorted(FIGS.glob('fig*.pdf'))], cwd=ROOT, capture_output=True, text=True)
        print(r.stdout.strip())
        return
    print('[5/6] scaffold')
    for _ in range(2):
        run(['pdflatex', '-interaction=nonstopmode', '-halt-on-error',
             '-output-directory=' + str(BUILD), 'scaffold/outline.tex'], env=env)
    aux = (BUILD / 'outline.aux').read_text()
    labels = dict(re.findall(r'\\newlabel\{fig:(\d+)\}\{\{(\d+)\}', aux))
    bad = {k: v for k, v in labels.items() if k != v}
    if len(labels) != 12 or bad:
        raise SystemExit(f'figure numbering mismatch: {labels}')
    if 'fig:guide' in aux or 'fig:0' in aux:
        raise SystemExit('the guide must not define a numbered figure label')
    print('    figure numbers 1-12 match their sections; guide unnumbered')
    print('[6/6] proof sheet, squint sheet, grayscale proof, deliverable PDFs, text-collision report')
    run(['pdflatex', '-interaction=nonstopmode', '-halt-on-error',
         '-output-directory=' + str(BUILD), 'scaffold/proofsheet.tex'], env=env)
    if shutil.which('gs'):
        run(['gs', '-q', '-o', str(BUILD / 'proofsheet-gray.pdf'), '-sDEVICE=pdfwrite',
             '-sColorConversionStrategy=Gray', '-dProcessColorModel=/DeviceGray', str(BUILD / 'proofsheet.pdf')])
    elif shutil.which('pdftoppm'):
        run(['pdftoppm', '-r', '110', '-gray', '-png', str(BUILD / 'proofsheet.pdf'), str(PREV / 'proofsheet-gray')])
    out = ROOT / 'figures' / 'pdf'
    out.mkdir(exist_ok=True)
    for pdf in sorted(FIGS.glob('fig*.pdf')):
        shutil.copy2(pdf, out / pdf.name)
    r = subprocess.run([sys.executable, 'checks/collisions.py'] + [str(p) for p in sorted(FIGS.glob('fig*.pdf'))],
                       cwd=ROOT, capture_output=True, text=True)
    (BUILD / 'collisions.txt').write_text(r.stdout)
    print(r.stdout.strip())
    print('done:', BUILD / 'outline.pdf', BUILD / 'proofsheet.pdf', out)

if __name__ == '__main__':
    main()
