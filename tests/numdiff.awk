# numdiff.awk -- summarise the numerical difference between two gortt outputs.
#
#   awk -f numdiff.awk golden/NAME.out actual/NAME.out
#
# Invoked by run.sh only after cmp has already reported a difference, to say
# how big that difference is: a change in the last printed digit is a rounding
# artefact, a change in the third significant figure is a broken model.
#
# Fields that do not look like numbers (nan, inf, the {} and [] delimiters
# around the component spectra) are compared as strings and reported
# separately -- a BRF that has turned into "-nan" must never be dismissed as
# a small numerical drift.

function isnum( x ) {
	return x ~ /^[-+]?([0-9]+\.?[0-9]*|\.[0-9]+)([eE][-+]?[0-9]+)?$/
}

function abs( x ) { return x < 0 ? -x : x }

FNR == NR { g[FNR] = $0; ng = FNR; next }

{
	na = FNR
	if ( FNR > ng ) { extra++; next }
	nf = split( g[FNR], G )

	if ( nf != NF ) { shape = shape sprintf( "line %d: %d fields -> %d\n", FNR, nf, NF ) }

	for ( i = 1; i <= NF && i <= nf; i++ ) {

		if ( $i == G[i] ) continue

		if ( !isnum($i) || !isnum(G[i]) ) {
			nnonnum++
			if ( nnonnum <= 3 )
				nonnum = nonnum sprintf( "line %d field %d: %s -> %s\n", FNR, i, G[i], $i )
			continue
		}

		d = abs( $i - G[i] )
		if ( d > maxabs ) { maxabs = d; abl = FNR; abf = i; abg = G[i]; aba = $i }

		if ( G[i] != 0 ) {
			r = d / abs( G[i] )
			if ( r > maxrel ) { maxrel = r; rll = FNR; rlf = i; rlg = G[i]; rla = $i }
		}
	}
}

END {
	if ( na < ng ) printf "output is %d line(s) short of golden\n", ng - na
	if ( extra )   printf "output has %d line(s) more than golden\n", extra
	if ( shape )   printf "%s", shape

	if ( nnonnum ) {
		printf "%d non-numeric difference(s) -- NOT a rounding artefact:\n", nnonnum
		printf "%s", nonnum
		if ( nnonnum > 3 ) printf "  ... and %d more\n", nnonnum - 3
	}

	if ( maxabs > 0 )
		printf "max abs diff %.3g at line %d field %d: %s -> %s\n", maxabs, abl, abf, abg, aba
	if ( maxrel > 0 )
		printf "max rel diff %.3g at line %d field %d: %s -> %s\n", maxrel, rll, rlf, rlg, rla

	if ( !nnonnum && !maxabs && !shape && !extra && na >= ng )
		printf "differs only in whitespace or field layout\n"
}
