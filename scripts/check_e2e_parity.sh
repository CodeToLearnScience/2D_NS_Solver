#!/usr/bin/env bash
#
# End-to-end parity check: runs the SAME case through the legacy binary and
# the new ns_solver, then reports byte-identity and per-variable deviation.
#
#   scripts/check_e2e_parity.sh [iterations]     (default 200)
#
# Requires: build/release/2D_NS_Solver, build/release/ns_solver,
#           configs/RampEuler.toml. Run from the repo root.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
ITERS="${1:-200}"
RUN="$ROOT/build"
LEG_OUT="$ROOT/output"; NEW_OUT=$(mktemp -d)

# --- legacy -----------------------------------------------------------------
rm -rf "$LEG_OUT" && mkdir -p "$LEG_OUT"
python3 - "$ROOT/input/InputRampEuler.cfg" "$ROOT/input/generated_e2e.cfg" RegRampEuler \
          RampGrid24080.dat GridTopRampEuler24080.dat "240 80" "$ITERS" $((ITERS/2)) <<'PY'
import sys
src,dst,tc,mesh,topo,dims,mi,of = sys.argv[1:9]
K={"#Test case name":("s",tc),"#mesh file":("s",mesh),
   "#Mesh topology file containing boundary informations":("s",topo),
   "#grid size (Nx*Ny)":("d",dims),"#Maximum number of iterations":("s",mi),
   "#Output frequency. Output will be written":("s",of)}
lines=open(src).read().splitlines(); out=[]; i=0
while i<len(lines):
    ln=lines[i]; out.append(ln)
    hit=next((v for k,v in K.items() if ln.startswith(k)),None)
    if hit and i+1<len(lines):
        i+=1
        out.append(" ".join(hit[1].split()) if hit[0]=="d" else str(hit[1]))
    i+=1
open(dst,"w").write("\n".join(out)+"\n")
PY
(cd "$RUN" && ./release/2D_NS_Solver ../input/generated_e2e.cfg >/dev/null 2>&1)

# --- new stack --------------------------------------------------------------
./build/release/ns_solver configs/RampEuler.toml input/RampGrid24080.dat "$NEW_OUT" \
    "$ITERS" > /dev/null

# --- compare -----------------------------------------------------------------
python3 - "$ROOT" "$NEW_OUT" "$ITERS" <<'PY'
import sys, os
root,new_out,iters=sys.argv[1],sys.argv[2],sys.argv[3]
pairs=[(f"{root}/output/Residue_RegRampEuler240_80.dat",
        f"{new_out}/Residue_Ramp_KFDS_maxU_24080_2nd240_80.dat"),
       (f"{root}/output/RegRampEuler_IsoCont240_80Iter{iters}.dat",
        f"{new_out}/Ramp_KFDS_maxU_24080_2nd_IsoCont240_80Iter{iters}.dat")]
names=["xc","yc","rho","u","v","p","T","Mach"]
def load(p,skip):
    rows=[]
    for ln in open(p):
        t=ln.split()
        try: rows.append([float(x) for x in t])
        except ValueError:
            if not skip: continue
    return rows
res_l,res_n=pairs[0]
rl=[[float(x) for x in ln.split()] for ln in open(res_l) if ln.split() and ln.split()[0].isdigit()]
rn=[[float(x) for x in ln.split()] for ln in open(res_n) if ln.split() and ln.split()[0].isdigit()]
print(f"residue rows: legacy={len(rl)} new={len(rn)}")
if rl and rn and rl[-1][3]>0:
    print(f"final residual: legacy={rl[-1][3]:.4e} new={rn[-1][3]:.4e} ratio={rn[-1][3]/rl[-1][3]:.6f}")
sl,sn=pairs[1]
A=[[float(x) for x in ln.split()] for ln in open(sl) if not ln.startswith(('TITLE','VARIABLES','Zone'))]
B=[[float(x) for x in ln.split()] for ln in open(sn) if not ln.startswith(('TITLE','VARIABLES','Zone'))]
print(f"solution rows: {len(A)} vs {len(B)}")
for c in range(2,8):
    scale=max(max(abs(a[c]),abs(b[c])) for a,b in zip(A,B))
    worst=max((abs(a[c]-b[c]),a[0],a[1]) for a,b in zip(A,B))
    print(f"  {names[c]:5s} max-abs-dev={worst[0]:.3e} (column scale {scale:.3e}, "
          f"rel {worst[0]/max(scale,1e-300):.3e}) at node x={worst[1]:.4f} y={worst[2]:.4f}")
PY
rm -rf "$NEW_OUT"
echo "done."
