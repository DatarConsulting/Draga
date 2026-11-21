"""Optional reproducible workaround for the present unindexed TeX runtime.

Normal TeX Live users should run latexmk directly. This script creates only
local temporary runtime files and invokes installed tools; it installs nothing.
"""
from pathlib import Path
import os
import subprocess
import tempfile

root=Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix='hankel-tex-') as directory:
    temp=Path(directory)
    env=dict(os.environ)
    env['TEXMF']='{/usr/share/texlive/texmf-dist,/usr/share/texmf,/etc/texmf}'
    subprocess.run(['pdftex','-ini','-etex','-jobname=pdflatex','-progname=pdflatex',
                    '-interaction=nonstopmode','pdflatex.ini'],cwd=temp,env=env,check=True)
    maps=[Path('/usr/share/texmf/fonts/map/dvips/lm/lm.map')]
    maps+=sorted(Path('/usr/share/texlive/texmf-dist/fonts/map/dvips/amsfonts').glob('*.map'))
    (temp/'pdftex.map').write_text('\n'.join(p.read_text() for p in maps))
    env['TEXFORMATS']=str(temp)+':'
    env['TEXFONTMAPS']=str(temp)+':'
    subprocess.run(['latexmk','-pdf','-interaction=nonstopmode','-halt-on-error',
                    'Hankel_Unified_Manuscript.tex'],cwd=root,env=env,check=True)
