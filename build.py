"""Build web/index.html: inlines web/beam.js into web/template.html so the page is a single file."""
import os
ROOT = os.path.dirname(os.path.abspath(__file__))
tpl = open(os.path.join(ROOT, "web", "template.html"), encoding="utf-8").read()
js = open(os.path.join(ROOT, "web", "beam.js"), encoding="utf-8").read()
body = tpl.replace("/*BEAMJS*/", js).replace("<!--SITEBAR-->", "").replace('/*FALLBACKS*/["../examples/"]', '["../examples/"]')
page = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n<meta name="viewport" content="width=device-width,initial-scale=1,viewport-fit=cover">\n'
        '<style>body{margin:0}</style>\n</head>\n<body>\n' + body + "\n</body>\n</html>\n")
open(os.path.join(ROOT, "web", "index.html"), "w", encoding="utf-8").write(page)
print("built web/index.html")

# also make the page work when opened straight from disk (file://)
import sys as _sys
_sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from offline_snapshot import inline as _inline
_R = os.path.dirname(os.path.abspath(__file__))
_inline(os.path.join(_R, "web/index.html"), ['presets.json'], "../examples/", os.path.join(_R, "examples"))
