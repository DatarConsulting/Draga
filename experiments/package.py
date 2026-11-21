"""Build three PDFs and a single-entry-point arXiv source archive."""
from pathlib import Path
import os, shutil, subprocess, tempfile, zipfile, hashlib
from datetime import datetime, timezone

SOURCE_FILES=("main.tex","body.tex","preamble.tex","references.tex","proofs.tex",
              "measurements.tex","benchmark-details.tex","neurips_2026.sty",
              "figures/joint.pdf","figures/planner.pdf")

def compile_tex(source, target):
    source=Path(source);target=Path(target)
    env=dict(os.environ, SOURCE_DATE_EPOCH=str(int(datetime(2026,11,21,tzinfo=timezone.utc).timestamp())),
             FORCE_SOURCE_DATE="1",TZ="UTC")
    with tempfile.TemporaryDirectory() as tmp:
        tmp=Path(tmp)
        # One log per pass is kept in memory; print the tail on failure.
        for _ in range(3):
            res=subprocess.run(["pdflatex","-file-line-error","-interaction=nonstopmode","-halt-on-error",
                                "-output-directory="+str(tmp),source.name],cwd=source.parent,env=env,
                               stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,errors="replace")
            if res.returncode: raise RuntimeError(res.stdout[-6500:])
        log=(tmp/(source.stem+".log")).read_text(errors="replace")
        for marker in ("There were undefined references","There were undefined citations","multiply defined",
                       "Missing character:","Overfull \\hbox","Overfull \\vbox"):
            if marker in log: raise RuntimeError(f"{source.name}: {marker}")
        shutil.copyfile(tmp/(source.stem+".pdf"),target)
    return target

def write_zip(path, files):
    with zipfile.ZipFile(path,"w") as z:
        for arc,src in sorted(files.items()):
            item=zipfile.ZipInfo(arc,(2026,11,21,0,0,0));item.compress_type=zipfile.ZIP_DEFLATED
            item.external_attr=0o100644<<16
            z.writestr(item,Path(src).read_bytes())

def build(root,validate=True):
    root=Path(root).resolve(); p=root/"paper"; dest=root/"artifacts";dest.mkdir(exist_ok=True)
    for src,pdf in [('main.tex','DRAGA_Article.pdf'),('supplement-wrapper.tex','DRAGA_Mathematical_Supplement.pdf'),
                    ('review-wrapper.tex','DRAGA_Anonymous_Review.pdf')]:
        compile_tex(p/src,dest/pdf)
    arxiv=dest/"DRAGA_arXiv_Source.zip"
    include={name:p/name for name in SOURCE_FILES}
    include.update({"anc/qwen.json":root/"results/qwen.json",
                    "anc/planner.py":root/"experiments/planner.py",
                    "anc/figures.py":root/"experiments/figures.py",
                    "anc/REPRODUCING.txt":p/"ancillary-guide.txt"})
    write_zip(arxiv,include)
    if validate:
        # Rebuild using only files from the outgoing ZIP.
        with tempfile.TemporaryDirectory() as d:
            with zipfile.ZipFile(arxiv) as z: z.extractall(d)
            check=compile_tex(Path(d)/"main.tex",Path(d)/"rebuilt.pdf")
            if check.read_bytes()!=(dest/"DRAGA_Article.pdf").read_bytes():
                # pdfTeX's trailer ID may depend on the temporary location. Page rendering is invariant.
                def rendered(path,out):
                    subprocess.run(["pdftoppm","-r","40","-gray",str(path),str(out)],check=True)
                    return [hashlib.sha256(x.read_bytes()).hexdigest() for x in sorted(out.parent.glob(out.name+"-*.pgm"))]
                if rendered(check,Path(d)/"a")!=rendered(dest/"DRAGA_Article.pdf",Path(d)/"b"):
                    raise RuntimeError("Rebuilding the source ZIP changed the rendered pages")
    return dest

if __name__=="__main__":
    import sys
    for f in sorted(build(sys.argv[1]).iterdir()):print(f.name,f.stat().st_size)
