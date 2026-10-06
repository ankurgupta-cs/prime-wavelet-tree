#!/usr/bin/env bash
# cell.sh ACTION TAG -- the driver behind `make ACTION-TAG`.
#
#   ACTION  fetch | verify | certify | rebuild | decode | paths
#   TAG     000 001 100 101 010 011 110 111 111x, or oracle / text
#
#   fetch    download the cell's file into $DATA and check its size and SHA-256
#            (URL from the manifest, or $DATA_BASE_URL/<file> when that is set;
#            the hosted files are in the v1/ folder of the data directory)
#   verify   size + SHA-256 + a structural decode (cells/<dir>/decode.sh verify)
#   certify  decode and compare every Carmichael number with the oracle
#   rebuild  rebuild the file from the oracle and require the published bytes
#            (cells/<dir>/build.sh)
#   decode   write the decoded values to $WORK (see cells/<dir>/decode.sh)
#   paths    trees only: export the stored paths to $WORK
ROOT=$(cd "$(dirname "$0")/.." && pwd)
. "$ROOT/scripts/common.sh"
[ $# -eq 2 ] || die "usage: cell.sh fetch|verify|certify|rebuild|decode|paths TAG"
action=$1 tag=$2

case $tag in
    oracle) dir=$ROOT/oracle; man=$dir/ORACLE ;;
    text)   dir=$ROOT/oracle; man=$dir/TEXT ;;
    *)
        set -- "$ROOT"/cells/"$tag"-*
        [ $# -eq 1 ] && [ -d "$1" ] || die "unknown cell $tag (cells: 000 001 100 101 010 011 110 111 111x)"
        dir=$1; man=$dir/CELL ;;
esac

fetch() {
    local file url bytes sha out
    file=$(mget "$man" file); url=$(mget "$man" url)
    bytes=$(mget "$man" bytes); sha=$(mget "$man" sha256)
    if [ "$(mget "$man" hosted)" = no ] && [ "$tag" != text ]; then
        die "$tag is not hosted; rebuild it with: make rebuild-$tag"
    fi
    if [ -n "${DATA_BASE_URL:-}" ] && [ "$tag" != text ]; then url=${DATA_BASE_URL%/}/$file; fi
    case $url in
        http://*|https://*) ;;
        *) die "$tag: no download URL in $man (url = '$url'); set DATA_BASE_URL=... to fetch from a mirror" ;;
    esac
    mkdir -p "$DATA"
    out=$DATA/$file
    if [ -f "$out" ] && check_file "$out" "$bytes" "$sha" "$tag ($file)"; then return 0; fi
    command -v curl >/dev/null 2>&1 || die "curl not found"
    run curl -fL --retry 3 -C - -o "$out.partial" "$url" || die "download failed: $url"
    check_file "$out.partial" "$bytes" "$sha" "$tag ($file)" || exit 1
    mv -f "$out.partial" "$out"
}

if [ "$tag" = text ]; then
    case $action in
        fetch) fetch ;;
        verify|certify) check_file "$DATA/$(mget "$man" file)" "$(mget "$man" bytes)" "$(mget "$man" sha256)" "text table" ;;
        rebuild|decode) exec bash "$dir/decode.sh" decode ;;   # the text, regenerated from the oracle
        *) die "unknown action $action for the text" ;;
    esac
    exit $?
fi

case $action in
    fetch) fetch ;;
    verify|certify|decode|paths)
        exec bash "$dir/decode.sh" "$action" ;;
    rebuild)
        exec bash "$dir/build.sh" ;;
    *) die "unknown action $action" ;;
esac
