# VBG1: sorted lists as VByte gaps

The byte-aligned list cells use this container:

- **cell 000:** the n list (every Carmichael number, ascending);
- **cell 100:** the d_min list (d_min(n) for every n, ascending).

Writer: `list_build vbyte <in.raw> <out.vbg>` (`src/list_build.cpp`).
Reader: `list_decode` (`src/list_decode.cpp`). The input `in.raw` is a sorted
list of 16-byte little-endian values (`n.raw` or `d.raw`, see
`results/INTERMEDIATES`).

## Layout

All integers little-endian.

    offset 0   8 B  magic "CNVBG1\0\0"
    offset 8   8 B  u64 N, the record count
    offset 16  8 B  u64 P, the payload byte count
    offset 24  P B  VByte(v_0), VByte(v_1 - v_0), ..., VByte(v_{N-1} - v_{N-2})

- VByte (`src/vbyte.hpp`): base-128, low 7 bits first, the high bit set on
  every byte but the last of a value.
- The file is exactly 24 + P bytes. There is no checksum: the file's
  SHA-256 (in `cells/*/CELL`) identifies it.

## Decoder obligations

A decoder must produce exactly N values, every gap after the first >= 1
(strictly ascending values), consuming exactly P bytes. `list_decode` also
requires canonical VBytes (no zero-padded encodings), no 128-bit overflow,
at most one record per payload byte, and every value odd.

The CNs of cell 000 are the values themselves. For cell 100, each value is a
divisor d and n = d * (d^{-1} mod lambda(d)), which needs d factored
(`list_decode --kind d --emit-n` or `--certify`).
