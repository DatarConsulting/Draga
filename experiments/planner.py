"""Teacher-forced policy scoring on the complete binary gate-suite. CPU ONNX."""
from pathlib import Path
from itertools import product
import argparse, hashlib, json, math, platform, time, importlib.metadata
import numpy as np
from scipy.special import logsumexp, xlogy

MODEL = "onnx-community/Qwen2.5-0.5B-Instruct"
REVISION = "22942cb7d7ba4cc81bb4673549ca4d4614469b5e"
DIGESTS = {"model_q4.onnx": "09235a3b1c135cd04ef570e5053b5c079e028078a2cc5f76ba6a251e91bf3296",
           "tokenizer.json": "a8506e7111b80c6d8635951a02eab0f4e1a8e4e5772da83846579e97b16f61bf"}

def digest(path):
    with open(path, "rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()

def instance(m, descending):
    order = list(range(2**m))
    if descending:
        order.reverse()
    target = "".join(format(i, f"0{m}b") for i in order)
    direction = "descending" if descending else "ascending"
    prompt = (f"The robot must unlock all {2**m} gates in {direction} numerical order, "
              f"starting at gate {order[0]} and ending at gate {order[-1]}. "
              f"To unlock gate j, encode j as exactly {m} binary digits, most significant first "
              "(pad with leading zeros). For digit 0 send A; for digit 1 send B. "
              "A wrong action fails the mission. Output only the actions separated by spaces, "
              "then end the response immediately after the last gate. Do not print gate numbers.")
    return prompt, tuple(map(int, target))

def verify(actions, m, descending):
    if len(actions) != m * 2**m: return False
    rooms = [int("".join(str(a) for a in actions[j:j+m]), 2)
             for j in range(0, len(actions), m)]
    return rooms == sorted(range(2**m), reverse=descending)

def lower_kl(logp, logcap):
    if logcap >= 0 or logp <= logcap: return 0.0
    p, c = math.exp(logp), math.exp(logcap)
    return float(p*(logp-logcap) + xlogy(1-p, 1-p) - (1-p)*math.log1p(-c))

def bounds(m):
    d, K, L = 2, 2**m, m*2**m
    # The scalar bound includes termination and letter frequencies.
    one = L*(math.log(L)-math.log(L+1)-math.log(d)) - math.log(L+1)
    return {str(n): (one if n==1 else min(0., K*math.log(n/K))) for n in (1,2,4)}

class Reader:
    def __init__(self, model, tokenizer):
        import onnxruntime as ort
        from tokenizers import Tokenizer
        self.tokenizer=Tokenizer.from_file(str(tokenizer))
        options=ort.SessionOptions()
        options.intra_op_num_threads=3
        options.inter_op_num_threads=1
        self.session=ort.InferenceSession(str(model),sess_options=options,providers=["CPUExecutionProvider"])
        self.tokens=[self.tokenizer.encode(t,add_special_tokens=False).ids for t in (" A"," B")]
        assert all(len(t)==1 for t in self.tokens)
        self.tokens=[t[0] for t in self.tokens]+[self.tokenizer.token_to_id("<|im_end|>")]
        assert len(set(self.tokens))==3 and None not in self.tokens

    def logprobs(self, prompt, actions):
        text=("<|im_start|>system\nYou are a helpful assistant.<|im_end|>\n"
              "<|im_start|>user\n"+prompt+"<|im_end|>\n<|im_start|>assistant\n")
        prefix=self.tokenizer.encode(text,add_special_tokens=False).ids
        seq=prefix+[self.tokens[a] for a in actions]
        ids=np.asarray([seq],dtype=np.int64)
        feed={"input_ids":ids,"position_ids":np.arange(len(seq),dtype=np.int64)[None,:],
              "attention_mask":np.ones_like(ids)}
        for inp in self.session.get_inputs():
            if inp.name.startswith("past_key_values."):
                feed[inp.name]=np.empty((1,2,0,64),np.float32)
        logits=self.session.run(["logits"],feed)[0][0,len(prefix)-1:].astype(np.float64)
        selected=logits[:,self.tokens]
        return selected-logsumexp(selected,axis=1)[:,None], selected-logsumexp(logits,axis=1)[:,None]

def capped_hankel(reader,prompt,T=3):
    transitions={}
    for path in product(range(2),repeat=T-1):
        restricted,_=reader.logprobs(prompt,path)
        for j in range(len(path)+1):
            transitions.setdefault(path[:j],np.exp(restricted[j]))
    words=[p for j in range(T+1) for p in product(range(2),repeat=j)]
    mass={}
    for w in words:
        v=math.prod(transitions[w[:j]][a] for j,a in enumerate(w))
        mass[w]=v*(transitions[w][2] if len(w)<T else 1.0)
    assert abs(sum(mass.values())-1)<1.e-5
    h=np.array([[mass.get(u+v,0.) for v in words] for u in words])
    sv=np.linalg.svd(h,compute_uv=False)
    skew=np.linalg.norm(h-h.T,ord=2)/2
    return {"horizon":T,"words":{"".join(map(str,w)):float(p) for w,p in mass.items()},
            "sigma1":float(sv[0]),"sigma2":float(sv[1]),"skew_floor":float(skew),
            "ratio":float(skew/sv[1])}

def run(model,tokenizer,out):
    assert digest(model)==DIGESTS["model_q4.onnx"] and digest(tokenizer)==DIGESTS["tokenizer.json"]
    reader=Reader(model,tokenizer)
    rows=[]; t=time.monotonic()
    for m in (1,2,3,4):
        for desc in (False,True):
            prompt,actions=instance(m,desc)
            assert verify(actions,m,desc)
            cond,raw=reader.logprobs(prompt,actions)
            target=list(actions)+[2]
            logp=float(cond[np.arange(len(target)),target].sum())
            rawp=float(raw[np.arange(len(target)),target].sum())
            row={"m":m,"order":"descending" if desc else "ascending","prompt":prompt,
                 "target":"".join("AB"[a] for a in actions),"logp":logp,"raw_logp":rawp,
                 "greedy_success":bool(np.array_equal(np.argmax(cond,axis=1),target)),
                 "logcaps":bounds(m),"kl_floor":{k:lower_kl(logp,v) for k,v in bounds(m).items()},
                 "step_logp":cond[np.arange(len(target)),target].tolist()}
            if m==1: row["capped"]=capped_hankel(reader,prompt)
            rows.append(row)
            print(m,row["order"],"ln P =",round(logp,4),"KL",row["kl_floor"],flush=True)
    payload={"model":MODEL,"revision":REVISION,"sha256":DIGESTS,"protocol_tokens":reader.tokens,
             "threads":3,"platform":platform.platform(),"seconds":time.monotonic()-t,
             "versions":{p:importlib.metadata.version(p) for p in ["onnxruntime","tokenizers","numpy","scipy"]},
             "instances":rows}
    Path(out).write_text(json.dumps(payload,indent=2)+"\n")

if __name__=="__main__":
    p=argparse.ArgumentParser(); p.add_argument("model");p.add_argument("tokenizer");p.add_argument("output")
    a=p.parse_args(); run(a.model,a.tokenizer,a.output)
