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

`-lidar` does not read stdin, so every lidar case registers `-` as its input. It is also
about eleven times slower than a BRF run at the same `-nlayers`, because a waveform needs
the gap probabilities at every canopy layer rather than only at the ground; the cases that
are not about vertical resolution therefore use a coarse `-nlayers 8`.

`lidar-waveform` and `lidar-pulse` pin an invariant that is easy to break and invisible in
a byte comparison: **the pulse convolution conserves energy exactly.** Summing the second
column and multiplying by the `binsize` on the header line gives the total returned energy,
and that number must not change with `-lidar_pulse` or `-lidar_binsize`. It must also equal
`J0*rho_v*(1-P(0)) + Rg`, which is eq. 22 of Ni-Meister et al. (2001). Checked by hand:

```sh
GORTT_PREC=17 ../gortt -lidar -lidar_waveform -lidar_pulse 2.5 -lidar_binsize 0.1 \
  | awk 'NR==1{db=$6;next}{s+=$2*db} END{printf "%.15g\n", s}'
```

`P(z)` is not quite monotonic in height on a coarse grid — the default 15 layers give a
non-monotonic patch around z = 7-9 m, i.e. a small negative `dP/dz`. That is a
discretisation artefact of the `ds = dz` path-length binning behind `epgap`, not a modelling
error: it shrinks with `-nlayers` and is gone by 60 layers. `lidar-nlayers` exists to keep an
eye on it.

## -nlayers

`nlayers-noop` is the case that pins "the default is unchanged": `-nlayers 15` names the
built-in value explicitly, so `golden/nlayers-noop.out` must be **byte-identical** to
`golden/default.out`. `nlayers-brf` proves the option reaches the BRF and not only the lidar,
and `nlayers-vs-N` proves `-N` still works alongside it — `-N` is matched case-insensitively
at two characters, so if `-nlayers` is ever moved down the if-chain it gets swallowed and
that case changes.

## Input-validation cases

`badmissingval`, `badmissingval2`, `badnowl` and `badwlcount` cover crashes
that were fixed once the suite exposed them: a value-taking option given as
the last word on the command line used to read `argv[argc] == NULL`, and a
header line listing more wavebands than it declared used to write past the end
of the wavelength array before the count check downstream noticed.
