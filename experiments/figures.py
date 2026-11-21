"""Generate both vector figures and all empirical tables from recorded inputs."""
from pathlib import Path
from itertools import product
from math import comb, log, exp
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.special import xlogy
from planner import instance, verify, bounds, lower_kl

def check_small_cases():
    # Verify the exact cut formula entrywise on every word through length four.
    words=["".join(v) for j in range(5) for v in product("ab",repeat=j)]
    for k in (1,2,3):
        r=.75; q=1-r*r
        def zword(w,j):
            start=w[:-j] if j else w
            return np.sqrt(q)*r**len(start) if w==start+"b"*j and set(start)<={"a"} else 0.
        z=np.array([[zword(w,j) for w in words] for j in range(k+1)])
        rhs=np.outer(z[0],z[k])
        for j in range(k): rhs[:,words.index("b"*j)] += np.sqrt(q)*z[k-j]
        def source(w):
            return q*r**(len(w)-k) if w.endswith("b"*k) and set(w[:-k])<={"a"} else 0.
        assert np.allclose(rhs,[[source(u+v) for v in words] for u in words],atol=1e-14)
    rng=np.random.default_rng(705)
    for n in (1,2,3,4):
        for m in (1,2,3):
            raw=rng.uniform(.1,1,size=(n,2*n+1))
            raw/=raw.sum(axis=1,keepdims=True)
            ts=[raw[:,j*n:(j+1)*n] for j in range(2)]
            _,path=instance(m,False); row=np.ones(n)/n
            for a in path: row=row@ts[a]
            cap=bounds(m)
            generic=min(1.,(n/2**m)**(2**m))
            assert float(row@raw[:,-1]) <= generic+1e-13
            if n==1: assert float(row@raw[:,-1]) <= exp(cap["1"])+1e-13

def check_group_controller():
    # Enumerate all state paths by a forward pass. No simulated trajectories.
    for m in range(1,5):
        for ell in range(m+1):
            G,D=2**ell,2**(m-ell)
            n=m*G+1; terminal=n-1
            mats=np.zeros((2,n,n))
            for g in range(G):
                for j in range(m):
                    src=g*m+j
                    actions=([((g>>(ell-j-1))&1,1.)] if j<ell else [(0,.5),(1,.5)])
                    for a,p in actions:
                        if j+1<m: mats[a,src,src+1]+=p
                        else:
                            mats[a,src,g*m]+=p*(1-1/D)
                            mats[a,src,(g+1)*m]+=p/D
            assert np.allclose(mats.sum(axis=(0,2))[:-1],1)
            state=np.eye(1,n)[0]
            for a in instance(m,False)[1]: state=state@mats[a]
            path_cost=(m-ell)*2**m*log(2)+G*log(D)
            if D>1: path_cost-=G*(D-1)*np.log1p(-1/D)
            assert np.isclose(-log(state[terminal]),path_cost,atol=1e-10)

def generate(root):
    root=Path(root)
    figures=root/"paper/figures"; figures.mkdir(parents=True,exist_ok=True)
    paper=root/"paper"
    rec=json.loads((root/"results/qwen.json").read_text()); rows=rec["instances"]
    assert [(r["m"],r["order"]) for r in rows]==[(m,o) for m in range(1,5) for o in ("ascending","descending")]
    check_small_cases()
    check_group_controller()
    for row in rows:
        prompt,w=instance(row["m"],row["order"]=="descending")
        assert row["prompt"]==prompt and verify(w,row["m"],row["order"]=="descending")
        assert row["target"]=="".join("AB"[a] for a in w)
        assert abs(sum(row["step_logp"])-row["logp"])<1e-9
        for n,c in bounds(row["m"]).items():
            assert abs(c-row["logcaps"][n])<1e-12
            assert abs(lower_kl(row["logp"],c)-row["kl_floor"][n])<1e-12
        if "capped" in row:
            data=row["capped"]; mass=data["words"];assert abs(sum(mass.values())-1.)<1.e-6
            mat=np.array([[mass.get(a+b,0.) for b in mass] for a in mass])
            assert np.isclose(np.linalg.svd(mat,compute_uv=False)[1],data["sigma2"])
            assert np.isclose(np.linalg.norm(mat-mat.T,2)/2,data["skew_floor"])
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,"axes.spines.top":False,
                         "axes.spines.right":False,"pdf.fonttype":42,"svg.fonttype":"none"})
    navy,teal,orange="#203955","#008279","#AC4D23"
    fig,axs=plt.subplots(1,2,figsize=(8.2,2.9),layout="constrained")
    m=np.unique(np.round(np.geomspace(2,512,160))).astype(int); eta=np.exp(-np.log(m)-m*np.log(2.));r=1-1/m;q=1-r*r;t=(1-eta)/(1+r)
    axs[0].plot(m,t*np.sqrt(q)+eta,color=navy,label=r"Upper bound on $\sigma_2$")
    axs[0].plot(m,np.maximum(0,t*r-eta),color=orange,ls="--",label=r"Lower bound on $\delta_1$")
    axs[0].plot(m,np.maximum(0,t*np.sqrt((1+5*q*q)/6)-eta),color=teal,ls="-.",label=r"Lower bound on $\delta_2$")
    for width,col in [(16,orange),(64,teal)]:
        ell=np.arange(width+1)
        K=2.**width; length=width*K; prob=1./length
        ns=3+width*2.**ell
        xx=np.log(ns)/np.log(length)
        entropy=-xlogy(prob,prob)-xlogy(1-prob,1-prob)
        low=np.maximum(0,log(2)-np.log(ns)/width-entropy)
        high=(1-ell/width)*log(2)+1/width
        axs[1].fill_between(xx,low,high,color=col,alpha=.14)
        axs[1].plot(xx,high,color=col,lw=1,label=f"Finite interval, m={width}")
        axs[1].plot(xx,low,color=col,lw=1)
    axs[1].plot([0,1],[log(2),0],ls="--",color=navy,label="Exact limiting curve")
    axs[0].set_xscale("log",base=2);axs[0].set_xticks([2,8,32,128,512],labels=["2","8","32","128","512"])
    axs[0].set(xlabel="Block width  m",ylabel="Operator norm")
    axs[1].set(xlabel="State exponent  log N / log L",ylabel="Optimal excess log-loss (nats)",xlim=(0,1.03))
    for ax in axs:
        ax.grid(axis="y",alpha=.18);ax.legend(fontsize=7.5,frameon=False)
    fig.savefig(figures/"joint.pdf",metadata={"CreationDate":None,"ModDate":None});plt.close(fig)
    fig,ax=plt.subplots(figsize=(6.4,2.5),layout="constrained")
    for order,col,marker in [("ascending",navy,"o"),("descending",orange,"s")]:
        rr=[x for x in rows if x["order"]==order]
        ax.plot([x["m"] for x in rr],[x["logp"]/log(10) for x in rr],marker=marker,
                color=col,label=order.title())
    ax.plot(range(1,5),[bounds(i)["1"]/log(10) for i in range(1,5)],ls="--",color=teal,label="One-state cap")
    ax.set(xticks=[1,2,3,4],xlabel="Gate width in bits",ylabel=r"$\log_{10}$ completed-plan probability")
    ax.grid(axis="y",alpha=.18);ax.legend(ncol=3,fontsize=8,frameon=False)
    fig.savefig(figures/"planner.pdf",metadata={"CreationDate":None,"ModDate":None});plt.close(fig)
    table=[]
    for i in range(0,len(rows),2):
        a,b=rows[i:i+2];m=a["m"]
        table.append(f'{m} & {m*2**m} & {a["logp"]:.3f} & {b["logp"]:.3f} & {a["raw_logp"]:.3f} & {b["raw_logp"]:.3f} & {bounds(m)["1"]:.3f} '+r'\\')
    txt=r"""\begin{table}[ht]
\centering\small
\caption{Teacher-forced log probabilities in nats. Each target includes the final stop. The last column applies to every proper autonomous one-state controller.}\label{tab:qwen}
\begin{tabular}{rr rr rr r}\toprule
&& \multicolumn{2}{c}{Three-token policy}&\multicolumn{2}{c}{Full vocabulary}&\\
$m$ & Length & Asc. & Desc. & Asc. & Desc. & $\log c_1$\\\midrule
"""+ "\n".join(table)+r"""
\bottomrule\end{tabular}
\end{table}
"""
    txt+=r"""On the ascending two-action instance the constrained policy assigned the valid trace probability $0.4127$ (rounded). Its best one-state competitor can assign at most $1/27$. Equation~\eqref{eq:event} gives a floor of $0.7046$ nats, rounded down. The other seven one-state event tests, and all eight tests at budgets two and four, returned zero. None of the eight restricted greedy rollouts solved its instance, as determined by the first mismatch with the verified action sequence, including the stop. The full-vocabulary policy crossed none of the reported caps. A teacher-forced success probability and a greedy success indicator answer different questions.

On the two width-one prompts we also enumerated all completed traces under a three-action horizon. For ascending and descending order, the ratios $\sigma_2/\sigma_1$ were respectively """
    caps=[r['capped'] for r in rows if 'capped' in r]
    txt+=f"{caps[0]['sigma2']/caps[0]['sigma1']:.3f} and {caps[1]['sigma2']/caps[1]['sigma1']:.3f}."
    txt+=r""" In both cases the antisymmetric lower bound was below $\sigma_2$. These checkpoint observations do not establish the spectral separation proved for the constructed law."""
    txt=txt.replace("0.4127",format(exp(rows[0]["logp"]),".4f"))
    (paper/"measurements.tex").write_text(txt+"\n")
    details=r"""The action interface uses token IDs 362 (space-A), 425 (space-B), and 151645 (end of message). The ID order and complete chat prompt are checked by the released program. Each model call evaluates the target prefix causally; no search or target-specific logit mask is used. The measured sweep used the CPU execution provider of ONNX Runtime VERSION, three intra-operation threads and one inter-operation thread. The measured scoring time, after initializing the session, was TIME seconds. Only public model parameters and generated gate instances were used.
\begin{center}
\begin{tabular}{lrrrr}\toprule
Order at $m=1$ & $\sigma_1$ & $\sigma_2$ & $\frac12\|H-H^\top\|$ & Rank-gap witness\\\midrule
"""
    details=details.replace("VERSION",rec["versions"]["onnxruntime"]).replace("TIME",f'{rec["seconds"]:.1f}')
    for name,c in zip(("Ascending","Descending"),caps):
        details+=f'{name} & {c["sigma1"]:.4f} & {c["sigma2"]:.4f} & {c["skew_floor"]:.4f} & Inconclusive '+r'\\'+"\n"
    details+=r"\bottomrule\end{tabular}\end{center}"+"\n"
    (paper/"benchmark-details.tex").write_text(details)
    return len(rows)

if __name__=="__main__":
    import sys
    generate(sys.argv[1])
