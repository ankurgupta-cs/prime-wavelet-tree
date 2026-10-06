#!/usr/bin/env python3
"""Fetch every page of a running local server and check that each one renders.

    py -m http.server 8765 --bind 127.0.0.1 --directory docs     (in another shell)
    py site/serve_check.py http://127.0.0.1:8765/ docs

For every HTML file under the docs directory: HTTP 200, an HTML content type,
a body equal to the file on disk, and a clean parse (tags closed in order, one
h1, th scope, unique ids). Every local link target and every asset is fetched
too. Python 3.8+, standard library only.
"""

import os
import sys
from pathlib import Path
from urllib.parse import urljoin, urlsplit
from urllib.request import urlopen

sys.path.insert(0, str(Path(__file__).resolve().parent))
import build_site as B  # noqa: E402


def fetch(url):
    with urlopen(url, timeout=10) as r:
        return r.status, r.headers.get("Content-Type", ""), r.read()


def main():
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    base, docs = sys.argv[1].rstrip("/") + "/", Path(sys.argv[2])
    problems = []
    pages = sorted(p.relative_to(docs).as_posix() for p in docs.rglob("*.html"))
    targets = set()
    for rel in pages:
        status, ctype, body = fetch(base + rel)
        if status != 200 or "text/html" not in ctype:
            problems.append("%s: HTTP %s, %s" % (rel, status, ctype))
            continue
        if body != (docs / rel).read_bytes():
            problems.append("%s: served bytes differ from the file" % rel)
        text = body.decode("utf-8")
        s = B.scan_html(text)
        problems += ["%s: %s" % (rel, m) for m in s.problems]
        for _, _, href, _ in s.links:
            if not urlsplit(href).scheme:
                targets.add(urljoin(base + rel, href).split("#")[0])
        print("ok  %-22s %6d bytes  %-28s %3d links  title %r" % (rel, len(body), ctype, len(s.links), s.title))
    targets.add(base)                # the site root must serve index.html
    for url in sorted(targets):
        status, ctype, body = fetch(url)
        if status != 200 or not body:
            problems.append("%s: HTTP %s" % (url, status))
    print("fetched %d pages and %d distinct local link targets (assets and the site root included)" % (len(pages), len(targets)))
    for p in problems:
        print("PROBLEM: " + p)
    print("serve check: %d problems" % len(problems))
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
