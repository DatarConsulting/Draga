"""Deterministic compiler for supplied symbolic memory contracts.

No API, shell execution, automatic semantic labels, or benchmark answers.
Guarantee: exact future-output preservation relative to the supplied total
Moore machine, with identical unknown-event handling. Not extraction accuracy.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import re
from collections import defaultdict, deque
from pathlib import Path, PurePosixPath

VERSION = "0.1.0"
MAX_BYTES = 8 * 1024 * 1024
SAFE_ID = re.compile(r"[a-z][a-z0-9_-]{0,63}\Z")


def encoded(obj):
    return (json.dumps(obj, ensure_ascii=False, sort_keys=True, indent=2) + "\n").encode("utf-8")


def text(value, name, limit=20000):
    if not isinstance(value, str) or not value.strip() or len(value) > limit:
        raise ValueError(f"invalid {name}")
    if any(ord(c) < 32 and c not in "\n\t" for c in value):
        raise ValueError(f"control character in {name}")
    return value


def ident(value):
    if not isinstance(value, str) or not SAFE_ID.fullmatch(value):
        raise ValueError(f"invalid identifier {value!r}")
    return value


def checked_evidence(items, inputs):
    if not isinstance(items, list) or not items:
        raise ValueError("every record requires source spans")
    root = Path(inputs).resolve()
    if not root.is_dir():
        raise ValueError("missing corpus input directory")
    result = []
    for item in items:
        rel = PurePosixPath(item["trajectory"])
        if rel.is_absolute() or ".." in rel.parts or not rel.parts:
            raise ValueError("unsafe evidence path")
        cursor = root
        for part in rel.parts:
            cursor = cursor / part
            if cursor.is_symlink():
                raise ValueError("symlink evidence forbidden")
        path = root / rel
        if not path.resolve().is_relative_to(root) or not path.is_file():
            raise ValueError("evidence missing or outside inputs")
        # A citation is a byte anchor, not a semantic verification.
        raw = path.read_bytes()
        raw.decode("utf-8")
        lines = raw.splitlines(keepends=True)
        first, last = item["first_line"], item["last_line"]
        if (type(first) is not int or type(last) is not int
                or not 1 <= first <= last <= len(lines)):
            raise ValueError("invalid evidence line interval")
        result.append({"trajectory": str(rel), "first_line": first, "last_line": last,
                       "file_sha256": hashlib.sha256(raw).hexdigest(),
                       "span_sha256": hashlib.sha256(b"".join(lines[first-1:last])).hexdigest()})
    return sorted(result, key=lambda e: (e["trajectory"], e["first_line"], e["last_line"]))


def validate_machine(skill):
    ident(skill["id"])
    for name in ("title", "scope", "trigger", "unknown_action"):
        text(skill[name], name)
    events = skill["events"]
    if not isinstance(events, dict) or not events or len(events) > 64:
        raise ValueError("1..64 explicitly guarded event symbols required")
    for key, value in events.items():
        ident(key); text(value, "event guard")
    states = skill["states"]
    if not isinstance(states, dict) or not 1 <= len(states) <= 2048:
        raise ValueError("1..2048 states required")
    for key, state in states.items():
        ident(key); text(state["output"], "state output")
        if set(state["transitions"]) != set(events):
            raise ValueError("machine must be total on its declared alphabet")
        if any(dest not in states for dest in state["transitions"].values()):
            raise ValueError("unknown transition destination")
    if skill["initial"] not in states:
        raise ValueError("invalid initial state")


def minimize(skill):
    """Reachability plus Moore partition refinement, with checked witness."""
    validate_machine(skill)
    states, alphabet = skill["states"], sorted(skill["events"])
    queue = deque([skill["initial"]]); reachable = set(queue)
    while queue:
        s = queue.popleft()
        for a in alphabet:
            t = states[s]["transitions"][a]
            if t not in reachable:
                reachable.add(t); queue.append(t)
    order = sorted(reachable)
    def label(signatures):
        ids = {}; out = {}
        for state in order:
            sig = signatures[state]
            if sig not in ids:
                ids[sig] = len(ids)
            out[state] = ids[sig]
        return out
    block = label({s: states[s]["output"] for s in order})
    rounds = 0
    while True:
        nxt = label({s: (states[s]["output"], tuple(block[states[s]["transitions"][a]]
                                                    for a in alphabet)) for s in order})
        rounds += 1
        if nxt == block:
            break
        block = nxt
    quotient = {}
    for s in order:
        name = f"q{block[s]}"
        quotient[name] = {"output": states[s]["output"],
                          "transitions": {a: f"q{block[states[s]['transitions'][a]]}" for a in alphabet}}
    mapping = {s: f"q{block[s]}" for s in order}
    witness = {"state_map": mapping, "input_states": len(states),
               "reachable_states": len(reachable), "output_states": len(quotient),
               "refinement_rounds": rounds}
    machine = {key: skill[key] for key in ("id", "title", "scope", "trigger", "events", "unknown_action")}
    machine.update(initial=mapping[skill["initial"]], states=quotient)
    # Do not trust the minimizer: verify the commuting diagram independently.
    for s, q in mapping.items():
        if states[s]["output"] != quotient[q]["output"]:
            raise ValueError("quotient witness changes an output")
        for a in alphabet:
            if mapping[states[s]["transitions"][a]] != quotient[q]["transitions"][a]:
                raise ValueError("quotient witness does not commute with a transition")
    if machine["initial"] != mapping[skill["initial"]]:
        raise ValueError("quotient witness changes the initial state")
    return machine, witness


def replay(machine, events):
    s = machine["initial"]; outputs = [machine["states"][s]["output"]]
    for a in events:
        if a not in machine["events"]:
            # Unknown inputs suspend the procedure; never silently continue.
            return outputs + [machine["unknown_action"]]
        s = machine["states"][s]["transitions"][a]
        outputs.append(machine["states"][s]["output"])
    return outputs


def build(spec, inputs, budget=MAX_BYTES):
    if type(budget) is not int or not 1 <= budget <= MAX_BYTES:
        raise ValueError("invalid persistent byte budget")
    if spec.get("schema_version") != "dhmc-0.1":
        raise ValueError("unsupported candidate schema")
    facts, skills, residuals = (spec.get(k, []) for k in ("facts", "skills", "residuals"))
    if not all(isinstance(x, list) for x in (facts, skills, residuals)) or not any((facts, skills, residuals)):
        raise ValueError("at least one grounded memory record required")
    seen = set(); checked = {}; audit = {}; sources = {}
    for item in facts + skills + residuals:
        rid = ident(item["id"])
        if rid in seen:
            raise ValueError("record IDs must be globally unique")
        seen.add(rid)
        checked[rid] = checked_evidence(item["evidence"], inputs)
        for ev in checked[rid]:
            sources[ev["trajectory"]] = ev["file_sha256"]
        text(item["scope"], "scope")
    files = {}; fact_groups = defaultdict(list)
    for item in facts:
        for key in ("key", "value"):
            text(item[key], key)
        fact_groups[item["scope"], item["key"]].append(item)
    fact_text = ["# Scoped facts and constraints", "", "Conflicts are unresolved unless source semantics settle them.", ""]
    for (scope, key), group in sorted(fact_groups.items()):
        values = {i["value"] for i in group}
        fact_text += [f"## {scope}: {key}", ""]
        if len(values) > 1:
            fact_text += ["**UNRESOLVED CONFLICT — establish the applicable context before acting.**", ""]
        for item in sorted(group, key=lambda i: i["id"]):
            fact_text += [f"- [{item['id']}] {item['value']}"]
        fact_text += [""]
    files["facts.md"] = "\n".join(fact_text).encode()
    index = ["# Datar memory", "", "Use only when the current scope and trigger match. Current user instructions take precedence.",
             "Read facts.md and exceptions.md before a procedure. Treat unresolved conflicts as questions to resolve.",
             "Event order matters. Select a transition only after establishing its guard; otherwise stop that procedure.",
             "Every new task starts each procedure at its initial state; no state is shared between tasks.", "",
             "- [Scoped facts](facts.md)", "- [Exceptions and uncertain lessons](exceptions.md)", ""]
    for skill in sorted(skills, key=lambda i: i["id"]):
        # Two distinct files is a minimum evidence gate, not proof of independence.
        if len({e["trajectory"] for e in checked[skill["id"]]}) < 2:
            raise ValueError("a reusable skill needs at least two distinct source trajectories")
        machine, witness = minimize(skill)
        sid = skill["id"]; audit[sid] = witness
        files[f"skills/{sid}/machine.json"] = encoded(machine)
        lines = [f"# {skill['title']}", "", f"Use when: {skill['trigger']}", f"Scope: {skill['scope']}",
                 f"Initial state: {machine['initial']}", f"If an event cannot be classified: {skill['unknown_action']}",
                 "", "Check exceptions.md first. Follow this table only while its scope holds.", "", "## Event guards", ""]
        lines += [f"- `{a}`: {machine['events'][a]}" for a in sorted(machine["events"])]
        for state, value in machine["states"].items():
            lines += ["", f"## State {state}", value["output"], ""]
            lines += [f"- On `{a}` → `{dest}`." for a, dest in value["transitions"].items()]
        lines += ["", f"Source anchors: {', '.join(e['trajectory']+':'+str(e['first_line'])+'-'+str(e['last_line']) for e in checked[sid])}",
                  "These citations do not prove extraction correctness. See provenance.json for hashes.", ""]
        files[f"skills/{sid}/SKILL.md"] = "\n".join(lines).encode()
        index += [f"- [{skill['title']}](skills/{sid}/SKILL.md): {skill['trigger']}"]
    exceptions = ["# Exceptions and uncertain lessons", "", "A failed trajectory does not identify a causal failure by itself.", ""]
    for item in sorted(residuals, key=lambda i: i["id"]):
        for key in ("trigger", "instruction"):
            text(item[key], key)
        exceptions += [f"## {item['id']} — {item['scope']}", f"When: {item['trigger']}", item["instruction"], ""]
    files["exceptions.md"] = "\n".join(exceptions).encode()
    files["MEMORY.md"] = ("\n".join(index)+"\n").encode()
    files["provenance.json"] = encoded({"sources": sources, "records": checked,
                                         "meaning": "Reference integrity only; no causal or semantic certificate."})
    files["compiler-certificate.json"] = encoded({"version": VERSION, "machines": audit,
        "scope": "Exact equivalence of supplied symbolic transducers on all declared event words; not of the LLM.",
        "semantic_extraction_verified": False, "format_limit_bytes": MAX_BYTES})
    # All retained data counts, including JSON, index, evidence anchors and certificate.
    total = sum(map(len, files.values()))
    if total > budget:
        raise ValueError(f"{total} retained bytes exceed {budget}; revise the candidate; no silent deletion")
    return files, {"bytes": total, "budget_bytes": budget, "file_count": len(files), "machines": audit}


def publish(spec_path, inputs, out, budget):
    out = Path(out)
    if out.is_symlink() or Path(spec_path).is_symlink():
        raise ValueError("symlink output or spec forbidden")
    files, report = build(json.loads(Path(spec_path).read_text()), inputs, budget)
    out.mkdir(parents=True, exist_ok=True)
    allowed = {Path(spec_path).resolve()} if Path(spec_path).parent.resolve() == out.resolve() else set()
    if any(p.resolve() not in allowed for p in out.iterdir()):
        raise ValueError("output must be empty or contain only the input candidate")
    for name, data in files.items():
        path = out / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)
    if allowed:
        Path(spec_path).unlink()
    actual = {p.relative_to(out).as_posix(): p.read_bytes() for p in out.rglob("*") if p.is_file()}
    if actual != files:
        raise ValueError("published byte tree differs from compiled export")
    return report


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--spec", type=Path, required=True)
    p.add_argument("--inputs", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--budget", type=int, default=MAX_BYTES)
    a = p.parse_args()
    print(json.dumps(publish(a.spec, a.inputs, a.out, a.budget), indent=2))


if __name__ == "__main__":
    main()
