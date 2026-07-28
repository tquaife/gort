#!/bin/sh
#
# Golden-file regression harness for gortt.
#
#   ./run.sh                compare every case against golden/, exit 1 on any diff
#   ./run.sh -b             (re)generate golden/ from the current binary
#   ./run.sh NAME...        restrict to the named cases
#   ./run.sh -b NAME...
#
# Every case captures stdout, stderr AND exit status, so the user-facing error
# paths are covered as well as the model output.  See ./cmds for the registry
# and ./numdiff.awk for reading a numerical failure.
#
# gortt is a deterministic stdin->stdout filter, so a byte comparison of its
# output is a complete check of the numerical path.  The one thing that can
# defeat it is print precision: "%f" gives 4 significant figures, which will
# not show a small regression.  See the GORTT_PREC block below.

cd "$(dirname "$0")" || exit 2

GORTT=../gortt
BLESS=0

[ "$1" = "-b" ] && { BLESS=1; shift; }

[ -x "$GORTT" ] || { echo "run.sh: $GORTT not built -- run make" >&2; exit 2; }

mkdir -p actual golden || exit 2

# ---------------------------------------------------------------------------
# Print precision.  17 significant digits round-trips an IEEE double exactly,
# so the goldens hold the model's actual values rather than a 4-figure
# rounding of them.  Measured: with the default "%f", a 1e-7 perturbation of
# favd was visible in only 1 of 30 BRF cases, and then only in the last digit.
# The gap-probability (-W) cases are unaffected -- they already print %0.40f.
# ---------------------------------------------------------------------------
GORTT_PREC=17
export GORTT_PREC

# ---------------------------------------------------------------------------
# Shared gap-probability file for the -P cases (@PROBS_DEFAULT@ in cmds).
# Default canopy geometry, so only spectral options may vary in those cases.
# ---------------------------------------------------------------------------
PROBS=actual/probs_default.dat
"$GORTT" -W < /dev/null > "$PROBS" 2> /dev/null || {
	echo "run.sh: could not generate $PROBS" >&2
	exit 2
}

npass=0
nfail=0
failed=

while read -r name input args; do

	case "$name" in '' | \#*) continue ;; esac

	# Restrict to the cases named on the command line, if any.
	if [ $# -gt 0 ]; then
		wanted=0
		for w in "$@"; do [ "$w" = "$name" ] && wanted=1; done
		[ "$wanted" = 1 ] || continue
	fi

	if [ "$input" = "-" ]; then
		stdin=/dev/null
	else
		stdin=cases/$input
		[ -f "$stdin" ] || { echo "run.sh: no such input: $stdin" >&2; exit 2; }
	fi

	args=$(printf '%s' "$args" | sed "s|@PROBS_DEFAULT@|$PROBS|g")

	# Unquoted $args on purpose: it must split into separate arguments.
	"$GORTT" $args < "$stdin" > "actual/$name.out" 2> "actual/$name.err"
	echo $? > "actual/$name.rc"

	if [ "$BLESS" = 1 ]; then
		cp "actual/$name.out" "golden/$name.out"
		cp "actual/$name.err" "golden/$name.err"
		cp "actual/$name.rc"  "golden/$name.rc"
		echo "bless $name"
		npass=$((npass + 1))
		continue
	fi

	diffs=
	for ext in out err rc; do
		if [ ! -f "golden/$name.$ext" ]; then
			diffs="$diffs $ext(no-golden)"
		elif ! cmp -s "golden/$name.$ext" "actual/$name.$ext"; then
			diffs="$diffs $ext"
		fi
	done

	if [ -z "$diffs" ]; then
		npass=$((npass + 1))
	else
		nfail=$((nfail + 1))
		failed="$failed $name"
		echo "FAIL $name --$diffs"
		if [ -f "golden/$name.out" ]; then
			awk -f numdiff.awk "golden/$name.out" "actual/$name.out" \
				| sed 's/^/       /'
		fi
	fi

done < cmds

if [ "$BLESS" = 1 ]; then
	echo "blessed $npass case(s)"
	exit 0
fi

if [ "$nfail" = 0 ]; then
	echo "PASS $npass case(s)"
	exit 0
fi

echo "FAILED $nfail of $((npass + nfail)) case(s):$failed"
echo "  inspect with: diff tests/golden/NAME.out tests/actual/NAME.out"
exit 1
