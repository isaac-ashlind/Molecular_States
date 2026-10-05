"""diffsheet.py BEFORE_DIR AFTER_DIR OUTPREFIX plate1 plate2 ... : before/after crops of every changed region at 300 dpi."""
import sys, os, subprocess
from PIL import Image, ImageChops, ImageFilter, ImageDraw
bd, ad, out = sys.argv[1:4]; plates = sys.argv[4:]
tmp = out + '_tmp'; os.makedirs(tmp, exist_ok=True)
crops = []
for n in plates:
    for tag, d in (('b', bd), ('a', ad)):
        subprocess.run(['pdftoppm', '-r', '300', '-png', '-singlefile', f'{d}/{n}.pdf', f'{tmp}/{tag}-{n}'], check=True)
    A = Image.open(f'{tmp}/b-{n}.png').convert('RGB'); B = Image.open(f'{tmp}/a-{n}.png').convert('RGB')
    if A.size != B.size: print(n, 'size changed', A.size, B.size); B = B.crop((0, 0) + A.size)
    d = ImageChops.difference(A.convert('L'), B.convert('L')).point(lambda v: 255 if v > 40 else 0).filter(ImageFilter.MaxFilter(41))
    small = d.resize((d.size[0] // 10, d.size[1] // 10)); W, H = small.size; px = small.load(); seen = set(); k = 0
    for y in range(H):
        for x in range(W):
            if px[x, y] and (x, y) not in seen:
                st = [(x, y)]; seen.add((x, y)); xs = []; ys = []
                while st:
                    cx, cy = st.pop(); xs.append(cx); ys.append(cy)
                    for nx, ny in ((cx+1, cy), (cx-1, cy), (cx, cy+1), (cx, cy-1)):
                        if 0 <= nx < W and 0 <= ny < H and px[nx, ny] and (nx, ny) not in seen: seen.add((nx, ny)); st.append((nx, ny))
                box = (max(0, min(xs)*10-40), max(0, min(ys)*10-40), min(A.size[0], max(xs)*10+50), min(A.size[1], max(ys)*10+50))
                crops.append((n, A.crop(box), B.crop(box))); k += 1
    print(n, k, 'regions')
if not crops: print('no changes'); sys.exit()
sheets = []; cur = []; h = 0
for c in crops:
    ch = c[1].size[1] + 26
    if cur and h + ch > 1800: sheets.append(cur); cur = []; h = 0
    cur.append(c); h += ch
sheets.append(cur)
for i, sh in enumerate(sheets):
    W = max(c[1].size[0]*2 + 24 for c in sh); H = sum(c[1].size[1] + 26 for c in sh)
    img = Image.new('RGB', (W, H), 'white'); dr = ImageDraw.Draw(img); y = 0
    for n, a, b in sh:
        dr.text((3, y+3), n + '  before | after', fill=(0, 0, 0)); img.paste(a, (0, y+20)); img.paste(b, (a.size[0]+24, y+20))
        dr.rectangle([a.size[0]+8, y+20, a.size[0]+14, y+20+a.size[1]], fill=(200, 0, 0)); y += a.size[1] + 26
    img.save(f'{out}{i}.png'); print(f'{out}{i}.png', img.size)
