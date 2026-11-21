"""Exhaustively compare packed C and exact integer updates in 72 archived machines.

For each policy we construct the binary transition graph reachable from its
prescribed initial packet. Equality on every edge gives equality on every
finite word by induction. Then all prefixes through depth 16 and the six
original length-100000 traces are evaluated against the binary64 reference.
This avoids repeatedly calling the same integer kernel for an already proved edge.
"""
import ctypes as C
from fractions import Fraction as F
from hashlib import sha256
import importlib.util
import json
import math
from pathlib import Path
import subprocess
from tempfile import TemporaryDirectory
import numpy as np
from archive_tools import Catalogue

HERE = Path(__file__).resolve().parent
array = np.ctypeslib.ndpointer


def require(test, description):
    if not test:
        raise ArithmeticError(description)


def libraries(code, tmp):
    native, diagnostic = tmp / "kernel.so", tmp / "trace.so"
    for src, dst in [(code / "kernel.c", native), (HERE / "trace.c", diagnostic)]:
        subprocess.run(["cc","-std=gnu11","-O2","-ffp-contract=off","-fPIC","-shared",
                        str(src),"-lm","-o",str(dst)],check=True)
    step = C.CDLL(str(native)).step
    i32 = array(dtype=np.int32,flags="C_CONTIGUOUS")
    f64 = array(dtype=np.float64,flags="C_CONTIGUOUS")
    step.argtypes = [C.POINTER(C.c_uint8),C.c_int,C.c_int,i32,i32,i32,C.c_int,C.c_int]
    step.restype = C.c_int
    trace = C.CDLL(str(diagnostic)).trace_maximum
    trace.argtypes = [C.c_size_t,array(dtype=np.uint8,flags="C_CONTIGUOUS"),C.c_size_t,
                      f64,f64,f64,i32,f64]
    trace.restype = C.c_double
    return step, trace


def graph(model, policy, original, step):
    alloc = model["best"][policy]
    frame = model["frames"][alloc["basis"]]
    bits = np.array(alloc["bits"],np.int32)
    shifts = np.array(alloc["shifts"],np.int32)
    coeff = np.array([[[int(F(x)*2**24) for x in row] for row in a] for a in frame["A"]],np.int32)
    initial = original.pack(alloc["initial_integers"],bits)
    states = [initial]
    ids = {initial:0}
    transition, answers = [], []
    vector = original.read_matrix(frame["beta"]).ravel()
    for raw in states:
        old = original.unpack(raw,bits)
        answers.append(float((np.asarray(old)*np.exp2(shifts.astype(float)))@vector))
        edges = []
        for letter in (0,1):
            buf = (C.c_uint8*len(raw)).from_buffer_copy(raw)
            require(step(buf,len(raw),len(bits),coeff,bits,shifts,letter,24)==0,"Native overflow")
            new = bytes(buf)
            exact = original.py_step(old,coeff,bits.tolist(),shifts.tolist(),letter,24)
            require(new==original.pack(exact,bits) and original.unpack(new,bits)==exact,
                    "Integer transition mismatch or loss of signed range")
            if new not in ids:
                ids[new] = len(states)
                states.append(new)
            edges.append(ids[new])
        transition.append(edges)
    return np.array(transition,np.int32), np.array(answers,np.float64)


def traces():
    rng = np.random.default_rng(991)
    # Preserve the historical generator's int8 sampling before converting the buffer.
    random = rng.integers(0,2,100000,dtype=np.int8)
    common = rng.integers(0,2,99998,dtype=np.int8)
    seqs = [np.zeros(100000),np.ones(100000),np.arange(100000)%2,random,
            np.r_[0,1,common],np.r_[1,0,common]]
    return [np.asarray(w,dtype=np.uint8) for w in seqs]


def replay(model, policy, original, step, trace, seqs, base):
    transitions, outputs = graph(model,policy,original,step)
    archive = json.loads((base/"results"/f"replay_{model['id']}_{policy}.json").read_text())
    require(len(transitions)==archive["reachable_machine_states"],"Reachable graph size")
    packets = list((base/"models").glob(f"{model['id']}_{policy}.bin"))
    require(len(packets)==1 and sha256(packets[0].read_bytes()).hexdigest()==archive["packet_sha256"],
            "Archived packet digest")
    a = np.array([original.read_matrix(v) for v in model["source_A"]])
    start = original.read_matrix(model["source_alpha"]).ravel()
    readout = original.read_matrix(model["source_beta"]).ravel()
    require(a.shape==(2,len(start),len(start)) and readout.shape==start.shape,"Reference dimensions")
    states, indices = start[None,:], np.zeros(1,np.int32)
    maximum = abs(float(start@readout)-outputs[0])
    visits = 0
    for _ in range(16):
        states = np.concatenate([states@operator for operator in a],axis=0)
        indices = np.concatenate([transitions[indices,i] for i in (0,1)])
        maximum = max(maximum,float(np.max(np.abs(states@readout-outputs[indices]))))
        visits += len(indices)
    for word in seqs:
        value = trace(len(word),word,len(start),start,a,readout,transitions,outputs)
        require(math.isfinite(value),"Non-finite reference trajectory")
        maximum = max(maximum,value)
        visits += len(word)
    require(visits==archive["exact_C_python_steps"],"Trace length or prefix count")
    require(math.isclose(maximum,archive["max_observed_abs_error"],rel_tol=1e-10,abs_tol=1e-12),
            f"Historical maximum error differs for {model['id']} / {policy}")
    require(maximum<=float(F(model["best"][policy]["total_bound"]))+1e-11,"Exceeded reference bound")
    return len(transitions), visits


def main():
    with TemporaryDirectory(prefix="draga-closure-") as temp, Catalogue() as cat:
        tmp = Path(temp)
        cat.extract("P4_LaTeX_Preuve_Relecture_Complet",tmp/"release")
        base = tmp/"release/P4_Manuscrit/supporting_code/breakthrough"
        code = base/"code"
        spec = importlib.util.spec_from_file_location("draga_original_replay",code/"replay.py")
        original = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(original)
        step, trace = libraries(code,tmp)
        machines = states = visits = 0
        seqs = traces()
        for path in sorted((base/"models").glob("*.json")):
            model = json.loads(path.read_text())
            for policy in ("uniform","mixed"):
                if model.get("best",{}).get(policy):
                    n,v = replay(model,policy,original,step,trace,seqs,base)
                    states += n
                    visits += v
                    machines += 1
        require(machines==72,"Policy count")
        print(json.dumps(dict(policies=machines,reachable_states_summed=states,
                              exact_binary_edges=2*states,reference_steps=visits),indent=2))


if __name__ == "__main__":
    main()
