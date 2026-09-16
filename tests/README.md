# Regression tests

`gortt` is a deterministic stdin→stdout filter, so a byte comparison of its
output is a complete check of the numerical path. These tests capture that
output for a set of representative runs and compare against it.

```sh
make test        # or: cd tests && ./run.sh
make bless       # re-capture the goldens -- see below
```

`run.sh NAME...` restricts the run to named cases. A failure prints the
largest absolute and relative deviation, so you can tell a last-digit rounding
artefact from a broken model.

## Layout

| path | |
|---|---|
| `cmds` | the case registry: name, input file, arguments |
| `cases/*.in` | model input (stdin) |
| `golden/*.{out,err,rc}` | expected stdout, stderr and exit status — **tracked** |
| `actual/` | last run's output — gitignored |
| `run.sh` | the harness |
| `bless.sh` | re-capture goldens |
| `numdiff.awk` | quantifies a numerical failure |

Exit status and stderr are captured as well as stdout, so the user-facing
error paths are covered too.

## Precision

`run.sh` sets `GORTT_PREC=17`, which makes `gortt` print reflectances with
enough significant digits to round-trip a double exactly. This matters: with
the default `%f` (4 significant figures), perturbing `favd` by 1 part in 10⁷
showed up in only **1 of 30** BRF cases, and then only in the last digit. At
`GORTT_PREC=17` the same perturbation was caught by **25 of 36**.

The gap-probability cases (`gortt -W`) print `%.17g` and need no knob; they
are the most sensitive probe in the suite and cover the whole expensive side
of the model. That format is deliberate — see the note on `useprobs-extreme`
below for what the previous `%0.40f` cost.

## Blessing

Re-capturing goldens is how a regression gets silently accepted, so only do it
deliberately: when adding a case, or on a commit that is *meant* to change the
numbers. In the latter case run `./run.sh` first and put the reported
deviation in the commit message.

## Known-bad behaviour pinned here

These cases record bugs, not correct behaviour. They exist so that fixing the
bug shows up as a deliberate, reviewable change.

- **`badanth`** — `gortt -u` documents `-Anth`, but the parser only matches
  `-canth` (at 3 characters, i.e. `-ca`), so the documented spelling is
  rejected as an unknown option. `prospect-anth` exercises the spelling that
  does work.

## Low sun

`lowsun` is a dense canopy (`-lambda 1.0`) with the sun at 0 and 5 degrees zenith, printing the
proportions. It exists because a change to `Po` in `gortt_kc_fFbeta` made Kc negative in
backscatter at low sun, and no existing case had such a row. Its sza 0 rows also hold β at its
overhead limit. See `KC-INVESTIGATION.md`, section 9.

## The -P round trip

`useprobs-extreme` and `useprobs-grazing` exist to hold `-P` to reproducing
the direct path exactly. They used to pin the opposite: `-P` returned `-nan`
above 89°, because the probability file was written with `%0.40f`, a *fixed*
point format, and the gap probabilities reach 4e-65 at grazing incidence.
Everything below 1e-40 was written as zero and `gortt_kuusk` took its log.

If you change what `-W` writes, these are the cases that will notice.

## Soil spectra

`soilspec` reads `cases/soil.spec`, which is deliberately on an irregular grid and overhangs
both ends of gortt's 400–2500 nm range — both of which `gortt_read_soil_lut` documents as
allowed. `soilwl.in` asks for two grid knots, one wavelength that falls between them, and both
endpoints.

`soilflat` and `soilconst` are the same run expressed two ways: a constant 0.25 file, and
`-alb_soil 0.25`. **Their goldens must stay byte-identical** — that is the check that the file
path and the constant path agree.

n.b. the soil file reader takes two numeric columns and has no comment syntax, so the `.spec`
fixtures carry no header.

## The black-soil limit

`blacksoil` (`-alb_soil 0`) and `blacksoil-limit` (`-alb_soil 1e-12`) bracket the one place the
energy balance divides by the soil reflectance. Their outputs agree to about 3e-12 relative,
which is what makes the black-soil case a limit rather than a special case. `blacksoil` used to
print `-nan` for the absorbed fractions.

## Lidar

`-lidar` does not read stdin, so every lidar case registers `-` as its input. It is also much
slower than a BRF run, because a waveform needs the gap probabilities at every canopy layer
rather than only at the ground, and it is the only mode whose cost rises with `-nlayers`
(roughly 0.4 s at 8 layers, 0.6 s at 15 and 2.0 s at 60, against 0.1 s for a BRF run at any
resolution). The cases that are not about vertical resolution therefore use a coarse
`-nlayers 8`.

`lidar-waveform` and `lidar-pulse` pin an invariant that is easy to break and invisible in
a byte comparison: **the pulse convolution conserves energy exactly.** Summing the second
column and multiplying by the `binsize` on the header line gives the total returned energy,
and that number must not change with `-lidar_pulse` or `-lidar_binsize`. It must also equal
`J0*rho_v*(1-P(0)) + Rg`, which is eq. 22 of Ni-Meister et al. (2001). Checked by hand:

```sh
GORTT_PREC=17 ../gortt -lidar -lidar_waveform -lidar_pulse 2.5 -lidar_binsize 0.1 \
  | awk 'NR==1{db=$6;next}{s+=$2*db} END{printf "%.15g\n", s}'
```

`P(z)` used to be not quite monotonic in height on a coarse grid — the default 15 layers gave
a small negative `dP/dz` around z = 7-9 m — and the profile wobbled from layer to layer by
about a quarter of the peak `dP/dz`. Both were quadrature artefacts, and both are gone: the
quadratures behind `epgap` no longer step with `nlayers`. `lidar-nlayers` exists to keep an
eye on it. The measure worth re-running after any change to `gortt_pn_kopen.c` is the rms
second difference of `dP/dz` over its peak, which `figures/quadrature_convergence.py` prints:

| `-nlayers` | 20 | 40 | 80 | 160 |
|---|---:|---:|---:|---:|
| before | 24.1% | 22.3% | 14.3% | 8.9% |
| now | 4.1% | 1.0% | 0.3% | 0.1% |

The `badlidaropt*` cases cover a hang: a `-lidar_*` option given without `-lidar` used to be
inert, and an inert `gortt` falls through to the BRF path and blocks on stdin, producing no
output and no message. The golden suite can only see that as an exit status, which is why
there is one case per option rather than one for the family — each sets its own flag, and a
missed one would be silent. `lidar-optorder` pins the other half: the check is made after the
argv loop, so `-lidar_wl 1550 ... -lidar` must give byte-identical output to
`-lidar -lidar_wl 1550 ...`, i.e. `golden/lidar-optorder.out` must equal `golden/lidar-wl.out`.

## -nlayers and -nquad

`-nlayers` is the reporting grid, not an accuracy setting: it says how finely the lidar
profile is written out and nothing else. No quadrature in the model steps with it any more,
so **the BRF must not move at all when it changes**. `nlayers-noop-a`, `-b` and `-c` are the
same run at 15, 31 and 120 layers and their goldens must stay byte-identical to each other;
`probs-nlayers` says the same thing about the gap probability tables, which are the sharpest
probe in the suite, and must stay byte-identical to `probs`. These replaced the older
`nlayers-noop`/`nlayers-brf` pair, which only pinned the much weaker "`-nlayers 15` is the
built-in value".

`nlayers-vs-N` proves `-N` still works alongside it — `-N` is matched case-insensitively at
two characters, so if `-nlayers` or `-nquad` is ever moved down the if-chain it gets swallowed
and that case changes.

`-nquad` multiplies the internal quadrature orders. `nquad-converged` and `nquad-probs` pin
the high-order answer; the claim they support is that the defaults are converged, which is a
*tolerance* statement rather than a byte one and so cannot be expressed as a golden:

```sh
awk -f numdiff.awk golden/default.out golden/nquad-converged.out
```

is 1.4e-5 relative today, dominated by `NQ_CONTACT`. If that grows, the default orders in
`include/gortt.h` have become too low for this canopy.

## Input-validation cases

`badmissingval`, `badmissingval2`, `badnowl` and `badwlcount` cover crashes
that were fixed once the suite exposed them: a value-taking option given as
the last word on the command line used to read `argv[argc] == NULL`, and a
header line listing more wavebands than it declared used to write past the end
of the wavelength array before the count check downstream noticed.
