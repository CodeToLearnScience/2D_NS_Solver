#!/usr/bin/env bash
#
# Phase-2 gate: convert every legacy config that is compatible with the
# current positional reader into the canonical TOML schema, then verify each
# result loads and validates.
#
#   scripts/migrate_configs.sh
#
# Known-incompatible legacy files (they predate the current reader or carry
# unsupported flags) are skipped with an explanation -- see
# REFACTORING_PLAN.md, section "Legacy input-file inconsistencies".
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PRESET="${PRESET:-release}"
TOOL="$ROOT/build/$PRESET/cfg2toml"

if [[ ! -x "$TOOL" ]]; then
    echo "cfg2toml not built; run: cmake --build --preset $PRESET" >&2
    exit 1
fi

CONVERTIBLE=(
    InputBlasius.cfg            Blasius.toml
    InputFfsEuler.cfg           FfsEuler.toml
    InputHypersonicFlowEuler.cfg HypersonicFlowEuler.toml
    InputObliqueShockRefEuler.dat ObliqueShockRefEuler.toml
    InputRampEuler.cfg          RampEuler.toml
    InputSWBLI_NS.dat           SWBLI_NS.toml
)

KNOWN_INCOMPATIBLE=(
    "InputHypersonicFlowNS.dat     space_accuracy=2 is not supported by the legacy solver itself"
    "InputSWBLI_NS_Dimen.dat       field layout does not match the current positional reader"
    "InputSWBLI_NS_ND1.dat         old pre-disp/outp/tot_time format"
    "InputWedgeReflection.dat      old pre-disp/outp/tot_time format"
    "InputRiemannProbWedgeRefEuler.dat uses the separate Riemann-problem parser"
    "InputRiemannSlipFlowEuler.dat uses the separate Riemann-problem parser"
    "InputFlowParND.dat            flow-parameter reference only, not a solver input"
    "InputFlowParDimensional.dat   flow-parameter reference only, not a solver input"
)

pass=0
mkdir -p "$ROOT/configs"
for ((i = 0; i < ${#CONVERTIBLE[@]}; i += 2)); do
    src="$ROOT/input/${CONVERTIBLE[i]}"
    dst="$ROOT/configs/${CONVERTIBLE[i + 1]}"
    if "$TOOL" "$src" "$dst" 2>/dev/null; then
        echo "OK    ${CONVERTIBLE[i]} -> configs/${CONVERTIBLE[i + 1]}"
        pass=$((pass + 1))
    else
        echo "FAIL  ${CONVERTIBLE[i]}"
        "$TOOL" "$src" /dev/null 2>&1 | sed 's/^/      /'
        exit 1
    fi
done

echo
echo "Converted $pass/${#CONVERTIBLE[@]} convertible configs into configs/"
echo
echo "Known-incompatible legacy inputs (not converted):"
for line in "${KNOWN_INCOMPATIBLE[@]}"; do
    printf '  - %s\n' "$line"
done
