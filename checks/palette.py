#!/usr/bin/env python3
"""The sixteen-tone check: every color the plates and any extra PDF paint is one of the sixteen tones of
figures/shared/figurestyle.sty. The manuscript is not checked, since its text and math are plain black.

    python3 checks/palette.py                 # figures/pdf/fig*.pdf
    python3 checks/palette.py extra.pdf ...   # and these (the poster art, say)
    python3 checks/palette.py --sty S --plates D   # another palette or plate folder

The palette is every \\definecolor{name}{HTML}{RRGGBB} (or {rgb}, {gray}) in figurestyle.sty, and there must be
sixteen. Each PDF is read with the standard library only: its objects (object streams included) are parsed, its
content streams and Form XObjects are decompressed (FlateDecode) and interpreted, and the fill and stroke color in
force at every painting operator is recorded (g G rg RG k K cs CS sc SC scn SCN, color spaces tracked through the
resources, q and Q honored, black by default, text painted in its render mode).

The check fails when
  (1) a painted color is not within 1/255 per channel of one of the sixteen tones (plus 1e-5 for the decimals the PDF
      writes), or the PDF paints with a shading, a pattern or a raster image, which bring their own colors;
  (2) an ExtGState sets a blend mode (any /BM but /Normal) or an opacity (/CA or /ca below 1, a soft mask).
It prints a table, the tones as rows and the PDFs as columns, each cell the number of painting operations.
"""
import re, sys, zlib
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
STY = ROOT / 'figures' / 'shared' / 'figurestyle.sty'
PLATES = ROOT / 'figures' / 'pdf'
TOL = 1 / 255 + 1e-5
COUNT = 16   # the number of tones figurestyle.sty must define

# ------------------------------------------------------------------------------------------------ the palette
def read_palette(sty):
    tones = {}
    for name, model, spec in re.findall(r'^[^%\n]*?\\definecolor\{([^}]+)\}\{(HTML|rgb|RGB|gray)\}\{([^}]+)\}',
                                        sty.read_text(), re.M):
        spec = spec.strip()
        if model == 'HTML':
            rgb = tuple(int(spec[k:k + 2], 16) / 255 for k in (0, 2, 4))
        elif model == 'rgb':
            rgb = tuple(float(v) for v in spec.split(','))
        elif model == 'RGB':
            rgb = tuple(int(v) / 255 for v in spec.split(','))
        else:
            rgb = (float(spec),) * 3
        tones[name] = rgb
    return tones

def hexof(rgb):
    return '#' + ''.join(f'{min(255, max(0, round(c * 255))):02X}' for c in rgb)

# ------------------------------------------------------------------------------------------------ PDF objects
WS = set(b' \t\r\n\f\x00')
DELIM = set(b'()<>[]{}/%')

class Name(str):
    pass

class Ref(tuple):
    pass

NUM_RE = re.compile(rb'[+-]?(?:\d+\.\d*|\.\d+|\d+)')
REF_RE = re.compile(rb'(\d+)\s+(\d+)\s+R(?![^\s()<>\[\]{}/%])')
KW_RE = re.compile(rb'[^\s()<>\[\]{}/%]+')
OBJ_RE = re.compile(rb'(?<![0-9])(\d+)\s+(\d+)\s+obj\b')

def skip_ws(d, i):
    n = len(d)
    while i < n:
        c = d[i]
        if c in WS:
            i += 1
        elif c == 37:   # a comment
            while i < n and d[i] not in (10, 13):
                i += 1
        else:
            break
    return i

def parse_string(d, i):
    """A literal string starting at d[i] == '('; returns (bytes, end)."""
    out, depth, i = bytearray(), 1, i + 1
    while True:
        c = d[i]
        if c == 92:   # backslash
            e = d[i + 1]
            if e in b'01234567':
                j = i + 1
                while j < i + 4 and d[j] in b'01234567':
                    j += 1
                out.append(int(d[i + 1:j], 8) & 255); i = j; continue
            out += {110: b'\n', 114: b'\r', 116: b'\t', 98: b'\b', 102: b'\f'}.get(e, bytes([e]) if e not in (10, 13) else b'')
            i += 2; continue
        if c == 40:
            depth += 1
        elif c == 41:
            depth -= 1
            if depth == 0:
                return bytes(out), i + 1
        out.append(c); i += 1

def parse(d, i):
    """One PDF object at d[i:]; returns (object, end)."""
    i = skip_ws(d, i)
    c = d[i]
    if c == 60 and d[i + 1] == 60:   # <<
        i += 2; out = {}
        while True:
            i = skip_ws(d, i)
            if d[i] == 62 and d[i + 1] == 62:
                return out, i + 2
            k, i = parse(d, i)
            v, i = parse(d, i)
            out[str(k)] = v
    if c == 91:   # [
        i += 1; out = []
        while True:
            i = skip_ws(d, i)
            if d[i] == 93:
                return out, i + 1
            v, i = parse(d, i); out.append(v)
    if c == 47:   # /
        j = i + 1
        while j < len(d) and d[j] not in WS and d[j] not in DELIM:
            j += 1
        raw = re.sub(rb'#([0-9A-Fa-f]{2})', lambda m: bytes([int(m.group(1), 16)]), d[i + 1:j])
        return Name(raw.decode('latin-1')), j
    if c == 40:
        return parse_string(d, i)
    if c == 60:   # a hex string
        j = d.index(b'>', i)
        h = re.sub(rb'\s', b'', d[i + 1:j])
        return bytes.fromhex((h + b'0' * (len(h) % 2)).decode()), j + 1
    m = REF_RE.match(d, i)
    if m:
        return Ref((int(m.group(1)), int(m.group(2)))), m.end()
    m = NUM_RE.match(d, i)
    if m:
        t = m.group()
        return (float(t) if b'.' in t else int(t)), m.end()
    m = KW_RE.match(d, i)
    if m:
        t = m.group()
        return {b'true': True, b'false': False, b'null': None}.get(t, Name(t.decode('latin-1'))), m.end()
    raise ValueError(f'cannot parse at {i}')

class PDF:
    def __init__(self, path):
        self.path = Path(path)
        self.data = self.path.read_bytes()
        self.objs = {}      # number -> object
        self.streams = {}   # number -> raw stream bytes
        self._scan()
        self._expand_object_streams()

    def _length(self, v):
        if isinstance(v, Ref):
            m = re.search(rb'(?<![0-9])%d\s+%d\s+obj\s*(\d+)\s*endobj' % (v[0], v[1]), self.data)
            return int(m.group(1)) if m else None
        return v

    def _scan(self):
        d, pos = self.data, 0
        while True:
            m = OBJ_RE.search(d, pos)
            if not m:
                break
            num = int(m.group(1))
            try:
                val, p = parse(d, m.end())
            except (ValueError, IndexError):
                pos = m.end(); continue
            self.objs[num] = val
            p = skip_ws(d, p)
            if d.startswith(b'stream', p) and isinstance(val, dict):
                s = p + 6
                s += 2 if d[s:s + 2] == b'\r\n' else 1
                n = self._length(val.get('Length'))
                if not isinstance(n, int) or not d.startswith(b'endstream', skip_ws(d, s + n)):
                    n = d.index(b'endstream', s) - s
                self.streams[num] = d[s:s + n]
                pos = s + n
            else:
                pos = p

    def _expand_object_streams(self):
        for num, val in list(self.objs.items()):
            if isinstance(val, dict) and val.get('Type') == 'ObjStm':
                data = self.decode(num)
                if data is None:
                    continue
                n, first = val['N'], val['First']
                head = [int(t) for t in data[:first].split()]
                for k in range(n):
                    onum, off = head[2 * k], head[2 * k + 1]
                    if onum not in self.objs:
                        try:
                            self.objs[onum] = parse(data, first + off)[0]
                        except (ValueError, IndexError):
                            pass

    def get(self, v):
        """Resolve a reference (repeatedly)."""
        seen = 0
        while isinstance(v, Ref) and seen < 32:
            v = self.objs.get(v[0]); seen += 1
        return v

    def num(self, v):
        return v[0] if isinstance(v, Ref) else None

    def decode(self, num):
        """The decoded data of stream object num, or None when a filter is not FlateDecode."""
        raw, d = self.streams.get(num), self.objs.get(num)
        if raw is None:
            return None
        filters = self.get(d.get('Filter'))
        filters = [] if filters is None else (filters if isinstance(filters, list) else [filters])
        parms = self.get(d.get('DecodeParms'))
        parms = parms if isinstance(parms, list) else [parms] * max(1, len(filters))
        for f, p in zip(filters, parms):
            f = self.get(f)
            if f in ('FlateDecode', 'Fl'):
                try:
                    raw = zlib.decompress(raw)
                except zlib.error:
                    raw = zlib.decompressobj().decompress(raw)
                p = self.get(p) or {}
                if p.get('Predictor', 1) >= 10:
                    raw = unpredict(raw, p.get('Columns', 1) * p.get('Colors', 1) * p.get('BitsPerComponent', 8) // 8)
            else:
                return None
        return raw

def unpredict(data, cols):
    """PNG predictors (xref and object streams of some writers)."""
    out, prev, row = bytearray(), bytearray(cols), cols + 1
    for r in range(0, len(data) // row):
        ft, line = data[r * row], bytearray(data[r * row + 1:(r + 1) * row])
        for k in range(cols):
            a = line[k - 1] if k else 0
            b, c = prev[k], (prev[k - 1] if k else 0)
            if ft == 1: line[k] = (line[k] + a) & 255
            elif ft == 2: line[k] = (line[k] + b) & 255
            elif ft == 3: line[k] = (line[k] + (a + b) // 2) & 255
            elif ft == 4:
                pa, pb, pc = abs(b - c), abs(a - c), abs(a + b - 2 * c)
                line[k] = (line[k] + (a if pa <= pb and pa <= pc else b if pb <= pc else c)) & 255
        out += line; prev = line
    return bytes(out)

# ------------------------------------------------------------------------------------------------ content streams
TOKEN_RE = re.compile(rb'''(?P<ws>[\s\x00]+)|(?P<com>%[^\r\n]*)|(?P<dopen><<)|(?P<dclose>>>)|(?P<hex><[0-9A-Fa-f\s]*>)
  |(?P<aopen>\[)|(?P<aclose>\])|(?P<str>\()|(?P<name>/[^\s()<>\[\]{}/%]*)|(?P<num>[+-]?(?:\d+\.\d*|\.\d+|\d+)(?![^\s()<>\[\]{}/%]))
  |(?P<op>[^\s()<>\[\]{}/%]+)''', re.X)

def tokens(d):
    """Yield (operands, operator) pairs of a content stream; inline images come back as ('BI', dict)."""
    i, n, stack = 0, len(d), []
    while i < n:
        m = TOKEN_RE.match(d, i)
        if not m:
            i += 1; continue
        kind = m.lastgroup; i = m.end()
        if kind in ('ws', 'com'):
            continue
        if kind == 'num':
            t = m.group(); stack.append(float(t))
        elif kind == 'name':
            stack.append(Name(m.group()[1:].decode('latin-1')))
        elif kind == 'str':
            s, i = parse_string(d, m.start()); stack.append(s)
        elif kind == 'hex':
            stack.append(b'')
        elif kind in ('aopen', 'dopen'):
            stack.append(kind)
        elif kind == 'aclose':
            k = len(stack) - 1
            while k >= 0 and stack[k] != 'aopen':
                k -= 1
            arr = stack[k + 1:]; del stack[max(k, 0):]; stack.append(arr)
        elif kind == 'dclose':
            k = len(stack) - 1
            while k >= 0 and stack[k] != 'dopen':
                k -= 1
            items = stack[k + 1:]; del stack[max(k, 0):]
            stack.append({str(items[j]): items[j + 1] for j in range(0, len(items) - 1, 2)})
        else:
            op = m.group().decode('latin-1')
            if op == 'BI':   # inline image: key-value pairs up to ID, then binary data up to EI
                j = d.index(b'ID', i)
                parms = {}
                try:
                    items = []
                    for ops, o in tokens(d[i:j]):
                        items += ops + ([Name(o)] if o not in ('',) else [])
                    parms = {str(items[k]): items[k + 1] for k in range(0, len(items) - 1, 2)}
                except Exception:
                    pass
                e = re.compile(rb'[\s\x00]EI(?=[\s\x00]|$)').search(d, j + 3)
                i = e.end() if e else n
                yield [parms], 'BI'
                stack = []; continue
            yield stack, op
            stack = []

PAINT_FILL = {'f', 'F', 'f*'}
PAINT_STROKE = {'S', 's'}
PAINT_BOTH = {'B', 'B*', 'b', 'b*'}
TEXT_SHOW = {'Tj', 'TJ', "'", '"'}

class Walker:
    """Interpret the content of one PDF and tally the colors painted."""
    def __init__(self, pdf):
        self.pdf = pdf
        self.uses = {}        # rgb (rounded) -> [fill, stroke, text]
        self.problems = []    # (1): shadings, patterns, images, color spaces that cannot be read
        self.blends = []      # (2): graphics states in use that blend or thin the paint
        self.seen_forms = set()

    def space(self, name, res):
        """A color space as (family, number of components, extra)."""
        pdf = self.pdf
        if name in ('DeviceGray', 'G'): return ('gray', 1, None)
        if name in ('DeviceRGB', 'RGB'): return ('rgb', 3, None)
        if name in ('DeviceCMYK', 'CMYK'): return ('cmyk', 4, None)
        if name == 'Pattern': return ('pattern', 0, None)
        cs = pdf.get((pdf.get(res.get('ColorSpace')) or {}).get(name)) if isinstance(name, str) else pdf.get(name)
        return self.space_of(cs, res)

    def space_of(self, cs, res):
        pdf = self.pdf
        if isinstance(cs, Name) or isinstance(cs, str):
            return self.space(Name(cs), res)
        if not isinstance(cs, list) or not cs:
            return ('unknown', 0, cs)
        fam = pdf.get(cs[0])
        if fam == 'ICCBased':
            n = (pdf.get(cs[1]) or {}).get('N', 3)
            return {1: ('gray', 1, None), 3: ('rgb', 3, None), 4: ('cmyk', 4, None)}.get(n, ('unknown', n, cs))
        if fam == 'CalRGB': return ('rgb', 3, None)
        if fam == 'CalGray': return ('gray', 1, None)
        if fam == 'Pattern': return ('pattern', 0, None)
        if fam in ('Indexed', 'I'):
            base = self.space_of(pdf.get(cs[1]), res)
            look = pdf.get(cs[3])
            if isinstance(cs[3], Ref) and cs[3][0] in pdf.streams:
                look = pdf.decode(cs[3][0])
            return ('indexed', 1, (base, look))
        if fam in ('Separation', 'DeviceN'):
            alt = self.space_of(pdf.get(cs[2]), res)
            fn = pdf.get(cs[3])
            if isinstance(fn, dict) and fn.get('FunctionType') == 2:
                return ('tint', 1, (alt, fn))
            return ('unknown', 1, fam)
        return ('unknown', 0, fam)

    def to_rgb(self, sp, vals):
        fam, n, extra = sp
        try:
            if fam == 'gray': return (vals[0],) * 3
            if fam == 'rgb': return tuple(vals[:3])
            if fam == 'cmyk':
                c, m, y, k = vals[:4]
                return ((1 - c) * (1 - k), (1 - m) * (1 - k), (1 - y) * (1 - k))
            if fam == 'indexed':
                base, look = extra
                k = int(vals[0]) * base[1]
                return self.to_rgb(base, [b / 255 for b in look[k:k + base[1]]])
            if fam == 'tint':
                alt, fn = extra
                c0, c1, e = fn.get('C0', [0]), fn.get('C1', [1]), fn.get('N', 1)
                t = vals[0] ** e
                return self.to_rgb(alt, [a + t * (b - a) for a, b in zip(c0, c1)])
        except (IndexError, TypeError):
            pass
        return None

    def initial(self, sp):
        fam = sp[0]
        if fam == 'cmyk': return (0.0, 0.0, 0.0)
        if fam == 'indexed': return self.to_rgb(sp, [0])
        if fam == 'tint': return self.to_rgb(sp, [1])
        if fam in ('gray', 'rgb'): return (0.0, 0.0, 0.0)
        return None

    def tally(self, rgb, slot, where):
        if rgb == 'pattern':
            return
        if rgb is None:
            self.problems.append(f'{where}: paints in a color that cannot be read (a pattern or an unknown space)')
            return
        key = tuple(round(c, 5) for c in rgb)
        self.uses.setdefault(key, [0, 0, 0])[slot] += 1

    def run_page(self, contents, res, where):
        data = b''
        for c in (contents if isinstance(contents, list) else [contents]):
            if isinstance(c, Ref):
                part = self.pdf.decode(c[0])
                if part is None:
                    self.problems.append(f'{where}: a content stream that is not FlateDecode (object {c[0]})')
                    continue
                data += part + b'\n'
        self.run(data, res, where, self.fresh())

    @staticmethod
    def fresh():
        return {'fcs': ('gray', 1, None), 'scs': ('gray', 1, None), 'fill': (0.0, 0.0, 0.0), 'stroke': (0.0, 0.0, 0.0), 'Tr': 0}

    def run(self, data, res, where, gs):
        pdf, stack = self.pdf, []
        res = pdf.get(res) or {}
        for ops, op in tokens(data):
            if op == 'q':
                stack.append(dict(gs))
            elif op == 'Q':
                if stack: gs = stack.pop()
            elif op in ('g', 'G', 'rg', 'RG', 'k', 'K'):
                sp = {'g': ('gray', 1, None), 'rg': ('rgb', 3, None), 'k': ('cmyk', 4, None)}[op.lower()]
                rgb = self.to_rgb(sp, ops)
                if op.islower(): gs['fcs'], gs['fill'] = sp, rgb
                else: gs['scs'], gs['stroke'] = sp, rgb
            elif op in ('cs', 'CS'):
                sp = self.space(ops[-1], res) if ops else ('unknown', 0, None)
                if sp[0] == 'unknown':
                    self.problems.append(f'{where}: color space {ops[-1] if ops else "?"} cannot be read ({sp[2]})')
                if op == 'cs': gs['fcs'], gs['fill'] = sp, self.initial(sp)
                else: gs['scs'], gs['stroke'] = sp, self.initial(sp)
            elif op in ('sc', 'scn', 'SC', 'SCN'):
                fill = op.islower()
                sp = gs['fcs'] if fill else gs['scs']
                nums = [v for v in ops if isinstance(v, float)]
                if sp[0] == 'pattern' or any(isinstance(v, Name) for v in ops):
                    self.problems.append(f'{where}: paints with a pattern ({op} {" ".join(map(str, ops))})')
                    rgb = 'pattern'   # reported here, not again at each paint
                else:
                    rgb = self.to_rgb(sp, nums)
                if fill: gs['fill'] = rgb
                else: gs['stroke'] = rgb
            elif op == 'Tr':
                gs['Tr'] = int(ops[-1]) if ops else 0
            elif op in PAINT_FILL:
                self.tally(gs['fill'], 0, where)
            elif op in PAINT_STROKE:
                self.tally(gs['stroke'], 1, where)
            elif op in PAINT_BOTH:
                self.tally(gs['fill'], 0, where); self.tally(gs['stroke'], 1, where)
            elif op in TEXT_SHOW:
                if gs['Tr'] in (0, 2, 4, 6): self.tally(gs['fill'], 2, where)
                if gs['Tr'] in (1, 2, 5, 6): self.tally(gs['stroke'], 2, where)
            elif op == 'sh':
                self.problems.append(f'{where}: a shading ({ops[-1] if ops else "?"}), which paints its own colors')
            elif op == 'gs':
                eg = pdf.get((pdf.get(res.get('ExtGState')) or {}).get(ops[-1])) if ops else None
                bad = extgstate_problem(eg, pdf)
                if bad:
                    self.blends.append(f'{where}: graphics state /{ops[-1]} sets {bad}')
            elif op == 'BI':
                p = ops[0] if ops else {}
                if p.get('IM') is True or p.get('ImageMask') is True:
                    self.tally(gs['fill'], 0, where)
                else:
                    self.problems.append(f'{where}: an inline raster image')
            elif op == 'Do' and ops:
                ref = (pdf.get(res.get('XObject')) or {}).get(ops[-1])
                xo = pdf.get(ref)
                if not isinstance(xo, dict):
                    continue
                if xo.get('Subtype') == 'Form':
                    data2 = pdf.decode(pdf.num(ref)) if pdf.num(ref) is not None else None
                    if data2 is None:
                        self.problems.append(f'{where}: a form that is not FlateDecode (object {pdf.num(ref)})')
                        continue
                    self.run(data2, xo.get('Resources') or res, where, dict(gs))
                elif xo.get('Subtype') == 'Image':
                    if xo.get('ImageMask') is True:
                        self.tally(gs['fill'], 0, where)
                    else:
                        self.problems.append(f'{where}: a raster image (object {pdf.num(ref)})')

def extgstate_problem(eg, pdf):
    if not isinstance(eg, dict):
        return None
    bad = []
    bm = pdf.get(eg.get('BM'))
    if bm is not None:
        names = bm if isinstance(bm, list) else [bm]
        if any(pdf.get(b) not in ('Normal', 'Compatible') for b in names):
            bad.append('blend mode /' + '/'.join(str(pdf.get(b)) for b in names))
    for k in ('CA', 'ca'):
        v = pdf.get(eg.get(k))
        if isinstance(v, (int, float)) and v < 1:
            bad.append(f'opacity /{k} {v}')
    sm = pdf.get(eg.get('SMask'))
    if sm is not None and sm != 'None':
        bad.append('a soft mask')
    return ', '.join(bad) or None

def all_extgstates(pdf):
    """Every ExtGState in the file: dicts typed /ExtGState and every entry of an /ExtGState resource dict."""
    found = {}
    def visit(v, depth=0):
        if depth > 40: return
        if isinstance(v, dict):
            if v.get('Type') == 'ExtGState':
                found[id(v)] = v
            eg = pdf.get(v.get('ExtGState'))
            if isinstance(eg, dict):
                for e in eg.values():
                    e = pdf.get(e)
                    if isinstance(e, dict): found[id(e)] = e
            for w in v.values(): visit(w, depth + 1)
        elif isinstance(v, list):
            for w in v: visit(w, depth + 1)
    for v in pdf.objs.values():
        visit(v)
    return list(found.values())

def pages(pdf):
    """The page dicts in order with their inherited resources."""
    out = []
    cat = next((v for v in pdf.objs.values() if isinstance(v, dict) and v.get('Type') == 'Catalog'), None)
    def walk(node, res, depth=0):
        node = pdf.get(node)
        if not isinstance(node, dict) or depth > 64: return
        res = node.get('Resources', res)
        if node.get('Type') == 'Pages' or 'Kids' in node:
            for k in pdf.get(node.get('Kids')) or []:
                walk(k, res, depth + 1)
        else:
            out.append((node, res))
    if cat is not None:
        walk(cat.get('Pages'), None)
    if not out:   # no catalog found: every page object, resources up the parent chain
        for v in pdf.objs.values():
            if isinstance(v, dict) and v.get('Type') == 'Page':
                res, p = v.get('Resources'), v
                while res is None and isinstance(p, dict) and 'Parent' in p:
                    p = pdf.get(p['Parent']); res = p.get('Resources') if isinstance(p, dict) else None
                out.append((v, res))
    return out

def inspect(path):
    pdf = PDF(path)
    w = Walker(pdf)
    for k, (page, res) in enumerate(pages(pdf), 1):
        w.run_page(page.get('Contents'), res, f'page {k}')
    blends = list(w.blends)
    for eg in all_extgstates(pdf):
        bad = extgstate_problem(eg, pdf)
        if bad: blends.append(f'an ExtGState sets {bad}')
    return w.uses, w.problems, blends

# ------------------------------------------------------------------------------------------------ the report
def nearest(rgb, tones):
    name = min(tones, key=lambda n: max(abs(a - b) for a, b in zip(rgb, tones[n])))
    return name, max(abs(a - b) for a, b in zip(rgb, tones[name]))

def short(path):
    s = Path(path).stem
    m = re.match(r'fig(\d\d)', s)
    return m.group(1) if m else s[:8]

def main(argv):
    global STY, PLATES
    import argparse
    ap = argparse.ArgumentParser(description='The sixteen-tone check.')
    ap.add_argument('pdfs', nargs='*', help='extra PDFs to check (the poster art, say)')
    ap.add_argument('--sty', type=Path, default=STY)
    ap.add_argument('--plates', type=Path, default=PLATES)
    a = ap.parse_args(argv)
    STY, PLATES, argv = a.sty, a.plates, a.pdfs
    tones = read_palette(STY)
    failures = []
    if len(tones) != COUNT:
        failures.append(f'figurestyle.sty defines {len(tones)} tones, not {COUNT}: {", ".join(tones)}')
    if len(set(hexof(v) for v in tones.values())) != len(tones):
        failures.append('two tones of figurestyle.sty are the same color')
    paths = sorted(PLATES.glob('fig*.pdf')) + [Path(a) for a in argv]
    if not paths:
        failures.append(f'no plates in {PLATES}')
    cols, cells, offrows = [], {}, {}
    for p in paths:
        col = short(p); cols.append(col)
        uses, problems, blends = inspect(p)
        for rgb, (f, s, t) in uses.items():
            name, dist = nearest(rgb, tones)
            if dist <= TOL:
                cells[(name, col)] = cells.get((name, col), 0) + f + s + t
            else:
                h = hexof(rgb)
                row = offrows.setdefault(h, {'near': (name, dist), 'cols': {}})['cols']
                row[col] = row.get(col, 0) + f + s + t
                failures.append(f'(1) {p.name}: {h} painted {f + s + t} times (fill {f}, stroke {s}, text {t}), '
                                f'nearest {name} {hexof(tones[name])} at {dist * 255:.1f}/255')
        for pr in dict.fromkeys(problems):
            failures.append(f'(1) {p.name}: {pr}')
        for b in dict.fromkeys(blends):
            failures.append(f'(2) {p.name}: {b}')
    # the table
    w = max([len(n) for n in tones] + [10]) + 17
    print(f'palette: the {len(tones)} tones of figurestyle.sty, painting operations per PDF')
    cw = [max(6, len(c) + 2) for c in cols]
    print(' ' * w + ''.join(f'{c:>{k}}' for c, k in zip(cols, cw)))
    for name, rgb in tones.items():
        print(f'{name:<{w - 8}}{hexof(rgb):>8}' + ''.join(f'{cells.get((name, c), 0) or "·":>{k}}' for c, k in zip(cols, cw)))
    for h, row in sorted(offrows.items()):
        n, dist = row['near']
        lab = f'!! ~{n} +{dist * 255:.0f}'
        print(f'{lab[:w - 8]:<{w - 8}}{h:>8}' + ''.join(f'{row["cols"].get(c, 0) or "·":>{k}}' for c, k in zip(cols, cw)))
    if failures:
        print(f'\npalette: {len(failures)} failures')
        for f in failures:
            print('  ' + f)
        return 1
    print(f'\npalette: ok, {len(paths)} PDFs paint only the {COUNT} tones')
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
