"""Verify source identity, reference integrity, catalogue coverage and compiled numbering.

The three fixed chapter digests identify the book before it was split. Replacing
each proof input with its precise byte sequence must reproduce those digests.
--aux checks the numbered destinations of a compiled book. --pdf compares its
pages with the published book at 72 dpi and compares the extracted text as well.
"""
import argparse
from collections import Counter
import json
from pathlib import Path, PurePosixPath
import re
from archive_tools import Catalogue, sha256

ROOT = Path(__file__).resolve().parents[1]
BASE = {
    "operator-approximation": "33ae4a36e3d55d798edc3f2a6f2bdc4ab07960e88c40701874c64533e93e1622",
    "probability": "901505ff1dbea84964158362426802deab20c198f1e6480a294974d5c6b458ab",
    "executable-memory": "b92f01b6781d3bd4287dd66ec1a991de83b9c0f8ecfc6af6820430dc01d59992",
}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def corpus():
    inv = json.loads((ROOT / "inventory.json").read_text())
    require(BASE == inv["chapter_sha256"], "Historical chapter fingerprint changed")
    records = {p["path"]: p for p in inv["proofs"]}
    statements = {s["number"]: s for s in inv["statements"]}
    seen = []
    texts = {}
    indexed = []
    lean = (ROOT / "lean/Draga.lean").read_text()
    formal_names = set(re.findall(r"^theorem (\w+)", lean, re.M))
    real = (ROOT / "lean-real/FibonacciMemory.lean").read_text()
    formal_names.update("FibonacciMemory." + n for n in re.findall(r"^theorem (\w+)", real, re.M))
    require(len(records) == len(inv["proofs"]) == 197, "Explicit proof count")
    require(len(statements) == len(inv["statements"]) == 206, "Statement count")
    for slug, digest in BASE.items():
        wrapper = (ROOT.parent / "chapters" / (slug + ".tex")).read_text().strip()
        require(wrapper == r"\input{proofs/" + slug + "/chapter.tex}", f"Chapter entry {slug}")
        chapter = (ROOT / slug / "chapter.tex").read_text()
        for s in statements.values():
            if s["chapter"] != slug:
                continue
            line = chapter.splitlines()[s["line"] - 1]
            require(line.startswith(r"\begin{" + s["kind"] + "}"), "Bad anchor " + s["number"])
            if "formal_component" in s:
                part = s["formal_component"]
                require(bool(part["scope"]) and set(part["lean"]) <= formal_names, "Invalid Lean link " + s["number"])
            if s["proofs"]:
                for name in s["proofs"]:
                    require(records[name]["statement"] == s["number"], f"Wrong owner {name}")
                    indexed.append(name)
            else:
                evidence = s["evidence"]
                require(bool(evidence["statements"] or evidence.get("lines")), "Missing proof location")
                require(all(n in statements for n in evidence["statements"]), "Unresolved proof pointer")
                if "lines" in evidence:
                    first, last = evidence["lines"]
                    require(1 <= first <= last <= len(chapter.splitlines()), "Invalid derivation interval")
        def inline(match):
            name = match[1]
            block = (ROOT / name).read_bytes()
            require(block.endswith(b"\n"), "Fragment terminator " + name)
            body = block[:-1]
            rec = records[name]
            require(sha256(body) == rec["sha256"] and len(body) == rec["bytes"], "Changed proof " + name)
            require(body.startswith(br"\begin{proof}") and body.endswith(br"\end{proof}"), "Proof bounds " + name)
            refs = sorted(set(re.findall(r"\\(?:eqref|ref)\{([^}]+)\}", body.decode())))
            require(refs == rec["references"], "Fragment references " + name)
            seen.append(name)
            return body.decode()
        full = re.sub(r"\\input\{proofs/([^}]+)\}", inline, chapter)
        require(sha256(full.encode()) == digest, "Reconstructed chapter " + slug)
        texts[slug] = full
    require(sorted(indexed) == sorted(seen) == sorted(records), "Orphan or duplicate proof")
    text = "\n".join(texts.values())
    labels = re.findall(r"\\label\{([^}]+)\}", text)
    require(len(labels) == len(set(labels)) == 590, "Chapter labels")
    uses = set(re.findall(r"\\(?:ref|eqref|pageref|autoref)\{([^}]+)\}", text))
    require(uses <= set(labels), f"Unresolved labels {uses-set(labels)}")
    cites = {k.strip() for m in re.findall(r"\\cite\w*(?:\[[^\]]*\])*\{([^}]+)\}", text) for k in m.split(",")}
    bibliography = (ROOT.parent / "references.tex").read_text()
    keys = re.findall(r"\\bibitem(?:\[[^\]]*\])?\{([^}]+)\}", bibliography)
    require(cites <= set(keys) and len(keys) == len(set(keys)), "Bibliographic key integrity")
    result = dict(statements=len(statements), kinds=dict(Counter(s["kind"] for s in statements.values())),
                  explicit_proofs=len(seen), associated_statements=len(set(indexed_s["statement"] for indexed_s in records.values())),
                  chapter_labels=len(labels), referenced_labels=len(uses), cited_keys=len(cites), bibliography_entries=len(keys))
    return statements, result


def software():
    with Catalogue() as cat:
        seen, source_ids = set(), set()
        n_members = 0
        n_source = 0
        for origin in cat.origins.values():
            for nested in origin["embedded"].values():
                require(nested in cat.origins, "Unresolved nested release")
            for name, rec in origin["members"].items():
                n_members += 1
                path = PurePosixPath(name)
                require(not path.is_absolute() and ".." not in path.parts and "\\" not in name, "Unsafe archive path")
                if rec["sha256"] not in seen:
                    cat.data(rec)
                    seen.add(rec["sha256"])
                if rec["role"] == "source":
                    source_ids.add(rec["sha256"])
                    n_source += 1
        require(seen == set(cat.objects), "Orphan archive object")
        rows = json.loads((ROOT / "archive/catalogue.json").read_text())
        require(len(rows) == len(source_ids) and {r["sha256"] for r in rows} == source_ids, "Readable source catalogue")
        extension_counts = Counter()
        for row in rows:
            path = ROOT / row["path"]
            data = path.read_bytes()
            require(sha256(data) == row["sha256"], "Readable copy changed " + row["path"])
            require(cat.origins[row["origin"]]["members"][row["member"]]["sha256"] == row["sha256"], "Invalid source origin")
            if path.suffix == ".py":
                compile(data, str(path), "exec")
            extension_counts[path.suffix] += 1
        for item in json.loads((ROOT / "checks/entrypoints.json").read_text()):
            require(cat.origins[item["origin"]]["members"][item["entry"]]["sha256"] == item["sha256"], "Invalid replay entry")
        return dict(releases=len(cat.origins), members=n_members, distinct_objects=len(seen),
                    source_occurrences=n_source, unique_sources=len(rows), source_types=dict(extension_counts))


def compiled(statements, aux, pdf):
    result = {}
    if aux:
        text = aux.read_text()
        numbered = dict(re.findall(r"\\newlabel\{([^}]+)\}\{\{([^}]+)\}", text))
        checked = 0
        for number, s in statements.items():
            if s["label"]:
                require(numbered.get(s["label"]) == number, f"Compiled number {number} / {s['label']}")
                checked += 1
        result["numbered_source_labels"] = checked
    if pdf:
        import pymupdf
        with pymupdf.open(pdf) as new, pymupdf.open(ROOT.parent / "monograph.pdf") as reference:
            require(len(new) == len(reference) == 289, "Page count")
            require(new.get_toc() == reference.get_toc(), "PDF outline mismatch")
            for i, (p, q) in enumerate(zip(new, reference), 1):
                require(p.rect == q.rect and p.get_text() == q.get_text(), f"Page {i} text or geometry")
                require(p.get_pixmap().samples == q.get_pixmap().samples, f"Page {i} at 72 dpi")
            destinations = new.resolve_names()
            for n in statements:
                require("theorem."+n in destinations, f"Missing PDF destination {n}")
        result.update(pages=289, numbered_pdf_destinations=len(statements), raster_dpi=72,
                      equal_text=True, equal_rasters=True, equal_outline=True)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aux", type=Path)
    parser.add_argument("--pdf", type=Path)
    args = parser.parse_args()
    records, structure = corpus()
    print(json.dumps(dict(corpus=structure, archive=software(),
                          typesetting=compiled(records,args.aux,args.pdf)), indent=2))
