"""Inline src/*.css and src/*.js into self-contained pages (no network needed):
index.html -> metamaterial-workbench.html (FE workbench)."""
import re, pathlib
root = pathlib.Path(__file__).resolve().parent.parent
for src, dst in [('index.html', 'metamaterial-workbench.html')]:
    html = (root / src).read_text(encoding='utf-8')
    html = re.sub(r'<link rel="stylesheet" href="(src/[^"]+)">', lambda m: '<style>\n' + (root / m.group(1)).read_text(encoding='utf-8') + '\n</style>', html)
    html = re.sub(r'<script src="(src/[^"]+)"></script>', lambda m: '<script>\n' + (root / m.group(1)).read_text(encoding='utf-8').replace('</script', '<\\/script') + '\n</script>', html)
    (root / dst).write_text(html, encoding='utf-8')
    print('wrote', root / dst, len(html), 'bytes')
