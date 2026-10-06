# common.sh -- helpers sourced by scripts/*.sh, cells/*/build.sh, cells/*/decode.sh
# and oracle/*.sh. Bash 3.2+ (Linux, macOS, MSYS2 / Git Bash).
#
# Environment (all optional):
#   BIN      directory with the built tools          (default: <repo>/build)
#   DATA     downloaded files: the hosted cells, the oracle, the text
#                                                     (default: <repo>/data)
#   WORK     intermediates, rebuilt files, scratch     (default: <repo>/work)
#   ORC1     the oracle                               (default: $DATA/<the file named
#                                                     in oracle/ORACLE>)
#   THREADS  threads for extraction and factoring     (default: all hardware threads)
# Large files never go into the repository tree unless DATA / WORK point there
# (both are in .gitignore).

set -u
export LC_ALL=C

[ -n "${ROOT:-}" ] || { echo "common.sh: the calling script must set ROOT" >&2; exit 2; }
BIN=${BIN:-$ROOT/build}
DATA=${DATA:-$ROOT/data}
WORK=${WORK:-$ROOT/work}
if [ -z "${THREADS:-}" ]; then
    THREADS=$( (nproc || getconf _NPROCESSORS_ONLN || sysctl -n hw.ncpu) 2>/dev/null | head -n 1)
    [ -n "$THREADS" ] || THREADS=1
fi

die() { echo "error: $*" >&2; exit 2; }

# tool NAME : path of a built tool (NAME.exe on Windows)
tool() {
    if [ -f "$BIN/$1.exe" ]; then printf '%s' "$BIN/$1.exe"
    elif [ -x "$BIN/$1" ]; then printf '%s' "$BIN/$1"
    else die "tool $1 not found in $BIN (run make)"; fi
}

# sha256_of FILE : lowercase hex SHA-256
if command -v sha256sum >/dev/null 2>&1; then
    sha256_of() { sha256sum < "$1" | awk '{print $1}'; }
elif command -v shasum >/dev/null 2>&1; then
    sha256_of() { shasum -a 256 < "$1" | awk '{print $1}'; }
elif command -v openssl >/dev/null 2>&1; then
    sha256_of() { openssl dgst -sha256 < "$1" | awk '{print $NF}'; }
else
    die "no sha256sum, shasum or openssl found"
fi
bytes_of() { wc -c < "$1" | tr -d ' \t\r\n'; }

# mget MANIFEST KEY : value of "KEY = value" in a CELL / ORACLE manifest
mget() {
    [ -f "$1" ] || die "manifest $1 not found"
    awk -v k="$2" '{ sub(/\r$/, "") } /^[ \t]*#/ { next }
        { i = index($0, "="); if (i == 0) next
          key = substr($0, 1, i - 1); gsub(/[ \t]+$/, "", key); gsub(/^[ \t]+/, "", key)
          if (key == k) { v = substr($0, i + 1); gsub(/^[ \t]+|[ \t]+$/, "", v); print v; exit } }' "$1"
}

# the oracle: the file named in oracle/ORACLE, in $DATA, unless ORC1 is set
ORC1=${ORC1:-$DATA/$(mget "$ROOT/oracle/ORACLE" file)}

# run CMD... : echo the command, then run it (exit status passed through)
run() { printf '+ %s\n' "$*" >&2; "$@"; }

# check_file FILE BYTES SHA WHAT : size and SHA-256 against the manifest
check_file() {
    local f=$1 b=$2 h=$3 what=$4 gb gh
    [ -f "$f" ] || { echo "MISSING   $what: $f" >&2; return 1; }
    gb=$(bytes_of "$f")
    if [ "$gb" != "$b" ]; then
        echo "MISMATCH  $what: $f is $gb bytes, expected $b" >&2; return 1
    fi
    printf 'hashing   %s (%s bytes) ...\n' "$f" "$gb"
    gh=$(sha256_of "$f")
    if [ "$gh" = "$h" ]; then
        printf 'MATCH     %s  sha256 %s\n' "$what" "$gh"; return 0
    fi
    echo "MISMATCH  $what: sha256 $gh, expected $h" >&2; return 1
}

need_oracle() {
    [ -f "$ORC1" ] || die "oracle not found: $ORC1 (make fetch-oracle, or set ORC1=...)"
}

# cell_file MANIFEST : the cell's file, from $FILE, else $DATA/<file>, else $WORK/<file>
cell_file() {
    local name
    name=$(mget "$1" file)
    if [ -n "${FILE:-}" ]; then printf '%s' "$FILE"
    elif [ -f "$DATA/$name" ]; then printf '%s' "$DATA/$name"
    else printf '%s' "$WORK/$name"; fi
}
