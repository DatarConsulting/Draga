from pathlib import Path
import struct,json
from uniform import weights
ROOT=Path(__file__).resolve().parents[1]
rows=json.loads((ROOT/'results/VERIFIED.json').read_text())['configs']
for r in rows:
    p=r['p'];m,J,R,nxt,cost=weights(p)
    vals=[p,m,J,R,r['states']]+nxt+cost
    (ROOT/'results'/f'model_p{p}.bin').write_bytes(struct.pack('<'+'I'*len(vals),*vals))
