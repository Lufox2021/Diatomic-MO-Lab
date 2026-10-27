"""Assemble the self-contained index.html from the source files.

Run from the repository root:  python3 src/build.py
"""
import json, os
here = os.path.dirname(os.path.abspath(__file__))
read = lambda name: open(os.path.join(here, name), encoding='utf-8').read()
sto = json.loads(read('sto3g.json'))
hf = read('hf.js').replace('/*STO3G*/', 'const STO3G = ' + json.dumps(sto, separators=(',', ':')) + ';')
open(os.path.join(here, 'hf.built.js'), 'w', encoding='utf-8').write(hf)   # used by the Node tests
hf_page = hf.replace("if (typeof module !== 'undefined') module.exports = HF;", '')
body = read('template.html').replace('/*ENGINE*/', read('engine.js')).replace('/*HF*/', hf_page).replace('/*APP*/', read('app.js'))
head = ('<!doctype html>\n<html lang="en">\n<head>\n<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
        '<style>:root{color-scheme:light;padding-top:env(safe-area-inset-top,0px);padding-bottom:env(safe-area-inset-bottom,0px)}'
        'body{margin:0;font-size:14px}img{max-width:100%}[hidden]{display:none!important}</style>\n')
# move the <title>, meta and font links from the template into <head>
split = body.index('<style>')
body = body.replace('<meta name="claude:viewport-fit-sides" content="cover">\n', '')
split = body.index('<style>')
page = head + body[:split] + body[split:body.index('</style>') + 8] + '\n</head>\n<body>\n' + body[body.index('</style>') + 8:] + '\n</body>\n</html>\n'
open(os.path.join(here, '..', 'index.html'), 'w', encoding='utf-8').write(page)
print('wrote index.html', len(page), 'bytes')
