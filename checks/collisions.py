"""Text-collision check for the plates: parse `pdftotext -bbox` output and report overlapping word boxes
that belong to different lines (same-line neighbours are expected to touch).  Usage: collisions.py file.pdf ...
Exit status 1 when any overlap larger than the tolerance is found."""
import subprocess, sys, re, itertools

TOL = 1.0   # points of allowed overlap

def words(pdf):
    xml = subprocess.run(['pdftotext', '-bbox', pdf, '-'], capture_output=True, text=True).stdout
    out = []
    for m in re.finditer(r'<word xMin="([\d.]+)" yMin="([\d.]+)" xMax="([\d.]+)" yMax="([\d.]+)">([^<]*)</word>', xml):
        x0, y0, x1, y1 = map(float, m.groups()[:4])
        out.append((x0, y0, x1, y1, m.group(5)))
    return out

def overlap(a, b):
    dx = min(a[2], b[2]) - max(a[0], b[0])
    dy = min(a[3], b[3]) - max(a[1], b[1])
    return dx, dy

def check(pdf):
    ws = words(pdf)
    bad = []
    for a, b in itertools.combinations(ws, 2):
        dx, dy = overlap(a, b)
        if dx > TOL and dy > TOL:
            # same baseline (within 2 pt) and horizontally adjacent: a kerning touch, not a collision
            if abs(a[3] - b[3]) < 2.0 and dx < 3.0:
                continue
            bad.append((a, b, dx, dy))
    return ws, bad

if __name__ == '__main__':
    status = 0
    for pdf in sys.argv[1:]:
        ws, bad = check(pdf)
        name = pdf.split('/')[-1]
        if bad:
            status = 1
            print(f'{name}: {len(bad)} overlapping word pairs')
            for a, b, dx, dy in bad[:12]:
                print(f'   "{a[4]}" ({a[0]:.0f},{a[1]:.0f}) with "{b[4]}" ({b[0]:.0f},{b[1]:.0f})  overlap {dx:.1f} x {dy:.1f} pt')
        else:
            print(f'{name}: no text collisions ({len(ws)} words)')
    sys.exit(status)
