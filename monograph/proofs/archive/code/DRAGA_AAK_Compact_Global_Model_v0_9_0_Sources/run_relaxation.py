"""Symmetry-exact, order-two moment relaxation of the compact PMI model.

Numerical discovery only. Solver status and numerical duals are not certificates.
The lower_v08 constraint relies on the previously verified v0.8 certificate.
"""
import argparse,json,time
from pathlib import Path
import numpy as np
import scipy.sparse as sp
import cvxpy as cp
from compact_model import build,monomials,character,polynomial,multiply
ROOT=Path(__file__).resolve().parent

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--iterations',type=int,default=30000)
    parser.add_argument('--eps',type=float,default=2e-6)
    parser.add_argument('--stein-equalities',action='store_true')
    parser.add_argument('--trace-cuts',action='store_true')
    parser.add_argument('--warm-start-from',default='')
    parser.add_argument('--tag',default='')
    parser.add_argument('--time-limit',type=float,default=0);args=parser.parse_args()
    suffix=('_stein' if args.stein_equalities else '')+('_trace' if args.trace_cuts else '')
    if args.tag:suffix+='_'+args.tag
    tic=time.time();d=build(args.stein_equalities,args.trace_cuts);xs=d['variables'];chars=d['characters'];n=len(xs)
    bases={k:monomials(n,k) for k in (0,1,2,4)}
    moms=[m for m in bases[4] if character(m,chars)==0];index={m:i for i,m in enumerate(moms)}
    y=cp.Variable(len(moms));constraints=[y[index[()]]==1];blocks=[];eqmeta=[]
    def localize(name,entries,states,degree):
        pol=[[polynomial(p,xs) for p in row] for row in entries]
        basis=bases[(4-degree)//2]
        for ch in range(4):
            pairs=[(i,m) for i,c in enumerate(states) for m in basis if c^character(m,chars)==ch]
            if not pairs:continue
            size=len(pairs);rows=[];cols=[];vals=[]
            for i,(a,ma) in enumerate(pairs):
                for j,(b,mb) in enumerate(pairs):
                    ab=ma+mb
                    for e,c in pol[a][b].items():
                        key=tuple(sorted(ab+e));assert character(key,chars)==0
                        rows.append(i*size+j);cols.append(index[key]);vals.append(float(c))
            M=sp.coo_matrix((vals,(rows,cols)),shape=(size*size,len(moms))).tocsr()
            con=cp.reshape(M@y,(size,size),order='C')>>0;constraints.append(con)
            blocks.append({'name':name,'char':ch,'size':size,'pairs':pairs,'constraint':con})
    localize('moment',[[1]],[0],0)
    for name,M,states in d['matrices']:
        degree=max((max(map(len,polynomial(p,xs)),default=0) for p in M),default=0)
        localize(name,M.tolist(),states,degree)
    er=[];ec=[];ev=[]
    for i,p in enumerate(d['equalities']):
        terms=polynomial(p,xs);degree=max(map(len,terms));ch=character(next(iter(terms)),chars)
        assert all(character(e,chars)==ch for e in terms)
        for m in bases[4-degree]:
            if character(m,chars)!=ch:continue
            row=len(eqmeta);eqmeta.append((i,m))
            for e,c in terms.items():er.append(row);ec.append(index[multiply(m,e)]);ev.append(float(c))
    eqmat=sp.coo_matrix((ev,(er,ec)),shape=(len(eqmeta),len(moms))).tocsr()
    eqcon=eqmat@y==0;constraints.append(eqcon)
    info={'version':'0.9.0','status':'STARTED','certified':False,'stein_equalities':args.stein_equalities,
          'trace_cuts':args.trace_cuts,'variables':n,'moments':len(moms),
          'equalities':len(d['equalities']),'localized_equalities':len(eqmeta),
          'blocks':[{k:b[k] for k in ('name','char','size')} for b in blocks]}
    (ROOT/('relaxation_discovery'+suffix+'.json')).write_text(json.dumps(info,indent=2))
    print('MODEL',json.dumps({k:v for k,v in info.items() if k!='blocks'}),flush=True)
    ui=xs.index(d['u']);problem=cp.Problem(cp.Minimize(y[index[(ui,)]]),constraints)
    try:
        if args.warm_start_from:
            seed=np.load(ROOT/('relaxation_numeric'+args.warm_start_from+'.npz'))
            prev=json.loads((ROOT/('relaxation_discovery'+args.warm_start_from+'.json')).read_text())
            pdata,chain,inverse=problem.get_problem_data(cp.SCS)
            dual=[float(prev['normalization_dual'])]+seed['equality_dual'].tolist()
            for i,b in enumerate(blocks):
                D=seed['dual_'+str(i)].copy();D[np.triu_indices(len(D),k=1)]*=np.sqrt(2)
                dual.extend(D[np.triu_indices(len(D))].tolist())
            rawx=seed['moments'];rawy=np.array(dual);raws=pdata['b']-pdata['A']@rawx
            assert pdata['A'].shape==(len(rawy),len(rawx))
            residual=float(np.linalg.norm(pdata['A'].T@rawy+pdata['c'],ord=np.inf))
            assert residual<.05,(residual,'Warm-start ordering check failed')
            problem._solver_cache['SCS']={'x':rawx,'y':rawy,'s':raws}
            info['warm_start_stationarity_residual']=residual
        settings={'time_limit_secs':args.time_limit} if args.time_limit else {}
        val=problem.solve(solver='SCS',eps=args.eps,max_iters=args.iterations,verbose=False,
                          acceleration_lookback=10,normalize=True,**settings)
        info.update(status=problem.status,objective=val,
                    normalization_dual=float(constraints[0].dual_value),
                    iterations=problem.solver_stats.num_iters,wall_seconds=time.time()-tic,
                    certified=False)
        np.savez_compressed(ROOT/('relaxation_numeric'+suffix+'.npz'),moments=y.value,
                            equality_dual=eqcon.dual_value,
                            raw_x=problem.solver_stats.extra_stats['x'],
                            raw_y=problem.solver_stats.extra_stats['y'],
                            raw_s=problem.solver_stats.extra_stats['s'],
                            **{'dual_'+str(i):b['constraint'].dual_value for i,b in enumerate(blocks)})
    except Exception as e:info.update(status='ERROR',error=str(e),wall_seconds=time.time()-tic)
    (ROOT/('relaxation_discovery'+suffix+'.json')).write_text(json.dumps(info,indent=2))
    print('RESULT',json.dumps({k:v for k,v in info.items() if k!='blocks'}),flush=True)

if __name__=='__main__':main()
