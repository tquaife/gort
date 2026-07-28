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

The gap-probability cases (`gortt -W`) already print `%0.40f` and need no
knob; they are the most sensitive probe in the suite and cover the whole
expensive side of the model.

## Blessing

Re-capturing goldens is how a regression gets silently accepted, so only do it
deliberately: when adding a case, or on a commit that is *meant* to change the
numbers. In the latter case run `./run.sh` first and put the reported
deviation in the commit message.

## Known-bad behaviour pinned here

These cases record bugs, not correct behaviour. They exist so that fixing the
bug shows up as a deliberate, reviewable change.

- **`useprobs-extreme`** — `-P` returns `-nan` at 89° and 89.5° zenith, where
  the direct path returns finite values. `gortt.c` writes only 90 of the
  `nth = 91` theta rows to the probability file, so index 90 is read back as
  zero. Compare with the `extreme` case, which is the same input run directly.
- **`badanth`** — `gortt -u` documents `-Anth`, but the parser only matches
  `-canth` (at 3 characters, i.e. `-ca`), so the documented spelling is
  rejected as an unknown option. `prospect-anth` exercises the spelling that
  does work.

## Input-validation cases

`badmissingval`, `badmissingval2`, `badnowl` and `badwlcount` cover crashes
that were fixed once the suite exposed them: a value-taking option given as
the last word on the command line used to read `argv[argc] == NULL`, and a
header line listing more wavebands than it declared used to write past the end
of the wavelength array before the count check downstream noticed.
