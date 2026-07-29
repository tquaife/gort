#!/bin/sh
#
# Re-capture the golden files from the current binary.
#
# Only ever do this deliberately: either when adding a case, or on a commit
# that is MEANT to change the numbers.  In the latter case, record the size of
# the change (run.sh prints it before you bless) in the commit message.
#
#   ./bless.sh              re-capture every case
#   ./bless.sh NAME...      re-capture only the named cases

exec "$(dirname "$0")/run.sh" -b "$@"
