"""Colour-vision-deficiency previews of the rendered figures (QA tooling, not a build dependency).

Usage:  python3 checks/cvd_preview.py build/previews/fig03-physical-hilbert-space.png [...]
Writes <name>-protan.png, -deutan.png, -tritan.png next to each input, using the Brettel, Vienot
and Mollon (1997) dichromat simulation as implemented in the daltonlens package (sRGB input,
severity 1.0).  This is an established published model, not an ad hoc RGB filter; like every
such model it approximates average dichromat perception and does not certify any individual's.
Requires: daltonlens, numpy, Pillow (pip install daltonlens).
"""
import sys
from PIL import Image
import numpy as np
from daltonlens import simulate

sim = simulate.Simulator_Brettel1997()
for path in sys.argv[1:]:
    im = np.asarray(Image.open(path).convert('RGB'))
    for name, deficiency in (('protan', simulate.Deficiency.PROTAN),
                             ('deutan', simulate.Deficiency.DEUTAN),
                             ('tritan', simulate.Deficiency.TRITAN)):
        out = sim.simulate_cvd(im, deficiency, severity=1.0)
        Image.fromarray(out).save(path.replace('.png', f'-{name}.png'))
    print('cvd previews written for', path)
