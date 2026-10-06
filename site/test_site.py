#!/usr/bin/env python3
"""Run build_site.py under three configs and check what each one links to.

    py site/test_site.py [--keep DIR]

  a  github_owner, data_folder_url and data_base_url null   (no public links yet)
  b  site/site.json as it is                                (the release config)
  c  github_owner as in site.json, no data URLs

Each build must pass its own link check with zero problems. On top of that:
  a  no link to GitHub, to GitHub Pages or into the data directory; every such
     link is a marker
  b  repository, website, the data directory, every hosted file, every index file
     of the data folder, every specification and every cell folder are linked;
     no marker is left; the data folder files are written, and their SHA256SUMS
     equals the repository's byte for byte
  c  GitHub links present, nothing in the data directory, every download is a marker
Then negative tests: a broken local link, an unknown host, a bad owner, a plain-http
data URL, a website link in a non-canonical form, a README block without its data
link, a data_base_url that disagrees with results/cells.tsv and a stale page must
each make build_site.py fail; a null arXiv id must build cleanly; the other value of
repo_public must change the data folder wording (a GitHub link once public, plain
'will be released at' before) and leave docs/ as it is. The repository's
own docs/, README.md and site/README_header.md are never written. Python 3.8+,
standard library only.
"""

import argparse
import json
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(HERE))
import build_site as B  # noqa: E402

BUILD = HERE / "build_site.py"


def run(args):
    p = subprocess.run([sys.executable, str(BUILD)] + args, capture_output=True, text=True)
    return p.returncode, p.stdout + p.stderr


def outs(d):
    return ["--config", str(d / "site.json"), "--out", str(d / "docs"), "--readme-out", str(d / "README_header.md"),
            "--readme", str(d / "README.md"), "--data-folder", str(d / "data-folder")]


def hrefs(docs):
    out = {}
    for path in sorted(docs.rglob("*.html")):
        s = B.scan_html(path.read_text(encoding="utf-8"))
        out[path.relative_to(docs).as_posix()] = [h for _, _, h, _ in s.links]
    return out


def data_folder_ok(dfolder, repo, public):
    """The repository in the data folder files: a link ('is at') once public, plain text
    ('will be released at') and no GitHub link before."""
    pages = [dfolder / "index.html", dfolder / "v1" / "index.html"]
    texts = [p.read_text(encoding="utf-8") for p in pages] + [(dfolder / "v1" / "README.txt").read_text(encoding="utf-8")]
    links = [h for p in pages for _, _, h, _ in B.scan_html(p.read_text(encoding="utf-8")).links]
    gh = [h for h in links if "github" in h]
    if public:
        return repo in links and not any("will be released" in t for t in texts)
    return not gh and all("will be released at" in t for t in texts)


def markers(docs):
    return {p.relative_to(docs).as_posix(): p.read_text(encoding="utf-8").count('class="pending"')
            for p in sorted(docs.rglob("*.html"))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--keep", type=Path, help="keep the test builds in this directory")
    args = ap.parse_args()
    work = args.keep or Path(tempfile.mkdtemp(prefix="pwt-site-test-"))
    work.mkdir(parents=True, exist_ok=True)
    base = json.loads((HERE / "site.json").read_text(encoding="utf-8"))
    owner, folder, data = base["github_owner"], base["data_folder_url"], base["data_base_url"].rstrip("/")
    cfg0 = B.load_config(HERE / "site.json")
    model = B.load_data(ROOT, HERE / "cells.json", cfg0)
    hosted = [c["file"] for c in model["cells"] if c["hosted"]] + [f["file"] for f in model["files"]]
    index_files = base["data_index_files"]
    specs = sorted(set(base["repo_paths"]["specs"].values()))
    configs = {
        "a": {"github_owner": None, "data_folder_url": None, "data_base_url": None},
        "b": {},
        "c": {"data_folder_url": None, "data_base_url": None},
    }
    failures = []
    summary = []

    def expect(cond, name, what):
        if not cond:
            failures.append("%s: %s" % (name, what))

    for name, over in configs.items():
        cfg = dict(base)
        cfg.update(over)
        d = work / name
        if d.exists():
            shutil.rmtree(d)
        d.mkdir(parents=True)
        (d / "site.json").write_text(json.dumps(cfg, indent=2), encoding="utf-8")
        shutil.copy(ROOT / "README.md", d / "README.md")
        rc, log = run(outs(d))
        (d / "build.log").write_text(log, encoding="utf-8")
        expect(rc == 0, name, "build_site.py exited %d" % rc)
        expect("check: 0 problems" in log, name, "link check reported problems")
        docs = d / "docs"
        links = hrefs(docs)
        allh = [h for hs in links.values() for h in hs]
        readme = (d / "README_header.md").read_text(encoding="utf-8")
        rlinks = re.findall(r"\]\(([^)]*)\)", readme)
        ext = sorted({re.sub(r"^https://([^/]+).*$", r"\1", h) for h in allh + rlinks if h.startswith("http")})
        mk = markers(docs)
        nmk = sum(mk.values())
        rmk = readme.count("*(link at release)*")
        spliced = (d / "README.md").read_text(encoding="utf-8")
        expect(readme.rstrip("\n") in spliced, name, "README.md does not hold the generated links block")
        dfolder = d / "data-folder"
        if name == "a":
            bad = [h for h in allh + rlinks if "github" in h or h.startswith(folder)]
            expect(not bad, name, "links to GitHub or the data directory: %s" % bad[:5])
            expect(set(ext) <= {"arxiv.org", "blue.butler.edu", "creativecommons.org"}, name, "unexpected hosts %s" % ext)
            expect(mk.get("data.html", 0) >= len(hosted) + 3, name, "data page lacks download markers")
            expect(rmk == 5, name, "README markers %d, expected 5 (repository, website, data, BibTeX, "
                   "website pointer)" % rmk)
            expect(not dfolder.exists(), name, "data folder files written without data URLs")
        if name in ("b", "c"):
            repo = "https://github.com/%s/%s" % (owner, base["repo_name"])
            site = "https://%s.github.io/%s/" % (owner.lower(), base["repo_name"])
            expect(repo in links["index.html"], name, "index lacks the repository link")
            expect(site in links["index.html"], name, "index lacks the website link")
            expect(repo in rlinks and site in rlinks, name, "README lacks repository or website link")
            for s in specs:
                expect("%s/blob/%s/%s" % (repo, base["branch"], s) in links["formats.html"], name, "no link to spec %s" % s)
            for c in model["cells"]:
                page = "cells/%s.html" % c["tag"]
                expect("%s/tree/%s/cells/%s" % (repo, base["branch"], c["dir"]) in links[page], name,
                       "%s lacks its folder link" % page)
            expect("git clone %s.git" % repo in (docs / "index.html").read_text(encoding="utf-8"), name, "no clone command")
        if name == "b":
            for f in hosted + index_files:
                expect("%s/%s" % (data, f) in links["data.html"], name, "data page lacks the link to %s" % f)
            for c in model["cells"]:
                if c["hosted"]:
                    expect("%s/%s" % (data, c["file"]) in links["cells/%s.html" % c["tag"]], name,
                           "cell %s lacks its download link" % c["tag"])
            expect(folder in rlinks and folder in links["index.html"], name, "README or index lacks the data directory link")
            expect(nmk == 0, name, "markers left %s" % {k: v for k, v in mk.items() if v})
            expect(rmk == 0, name, "README markers %d, expected 0" % rmk)
            want = ["index.html", "v1/index.html", "v1/README.txt", "v1/MANIFEST.tsv", "v1/LICENSE-DATA.txt", "v1/SHA256SUMS"]
            got = sorted(p.relative_to(dfolder).as_posix() for p in dfolder.rglob("*") if p.is_file())
            expect(got == sorted(want), name, "data folder files %s, expected %s" % (got, sorted(want)))
            sums = dfolder / "v1" / "SHA256SUMS"
            expect(sums.exists() and sums.read_bytes() == (ROOT / "SHA256SUMS").read_bytes(), name,
                   "data folder SHA256SUMS differs from the repository's")
            man = (dfolder / "v1" / "MANIFEST.tsv").read_text(encoding="utf-8").splitlines()
            expect(len(man) == len(hosted) + 1 and all(f in "\n".join(man) for f in hosted), name,
                   "MANIFEST.tsv does not list every hosted file once")
            expect(data_folder_ok(dfolder, repo, base["repo_public"]), name,
                   "data folder wording does not match repo_public=%s" % base["repo_public"])
            rd = (ROOT / "docs").resolve()
            same = all((rd / p.relative_to(docs)).exists() and (rd / p.relative_to(docs)).read_bytes() == p.read_bytes()
                       for p in docs.rglob("*") if p.is_file())
            summary.append("b  output %s the repository's docs/" % ("equals" if same else "DIFFERS FROM (run make site)"))
        if name == "c":
            expect(not [h for h in allh + rlinks if h.startswith(folder)], name, "links into the data directory")
            expect(mk.get("data.html", 0) >= len(hosted) + 3, name, "data page lacks download markers")
            expect(rmk == 1, name, "README markers %d, expected 1 (data)" % rmk)
        summary.append("%s  exit %d  pages %d  links %d  external hosts %s  markers %d (README %d)"
                       % (name, rc, len(links), len(allh), ",".join(ext), nmk, rmk))
        for line in log.splitlines():
            if line.startswith("check:"):
                summary.append("   " + line)

    # negative tests: the checker must fail
    neg = work / "neg"
    if neg.exists():
        shutil.rmtree(neg)
    shutil.copytree(work / "b", neg)
    idx = neg / "docs" / "index.html"
    t = idx.read_text(encoding="utf-8")

    def check_neg():
        return run(outs(neg) + ["--check"])

    idx.write_text(t.replace('href="data.html"', 'href="dta.html"', 1), encoding="utf-8")
    rc, log = check_neg()
    expect(rc == 1 and "missing local page dta.html" in log, "negative", "broken local link not caught (exit %d)" % rc)
    summary.append("negative 1: broken local link -> exit %d (%s)" % (rc, "caught" if rc == 1 else "MISSED"))
    idx.write_text(t.replace('href="https://arxiv.org/abs/', 'href="https://arxiv.example.org/abs/', 1), encoding="utf-8")
    rc, log = check_neg()
    expect(rc == 1 and "is not one the config names" in log, "negative", "unknown host not caught (exit %d)" % rc)
    summary.append("negative 2: link to a host the config does not name -> exit %d (%s)" % (rc, "caught" if rc == 1 else "MISSED"))
    idx.write_text(t, encoding="utf-8")

    def bad_config(fname, **kw):
        cfg = json.loads((neg / "site.json").read_text(encoding="utf-8"))
        cfg.update(kw)
        (neg / fname).write_text(json.dumps(cfg), encoding="utf-8")
        o = neg / ("out_" + fname.replace(".json", ""))
        return run(["--config", str(neg / fname), "--out", str(o / "docs"), "--readme-out", str(o / "R.md"),
                    "--readme", "none", "--data-folder", str(o / "data-folder")]), o

    (rc, log), _ = bad_config("owner.json", github_owner="bad owner")
    expect(rc == 2, "negative", "bad owner not rejected (exit %d)" % rc)
    summary.append("negative 3: github_owner 'bad owner' -> exit %d (%s)" % (rc, "rejected" if rc == 2 else "MISSED"))

    # 4: a plain-http data URL is refused (browsers block http downloads from an https page)
    (rc, log), _ = bad_config("http.json", data_base_url=data.replace("https://", "http://"))
    expect(rc == 2 and "https" in log, "negative", "http data URL not rejected (exit %d)" % rc)
    summary.append("negative 4: http data_base_url -> exit %d (%s)" % (rc, "rejected" if rc == 2 else "MISSED"))

    # 5: right host, non-canonical form (the website link without its trailing slash)
    site = "https://%s.github.io/%s/" % (owner.lower(), base["repo_name"])
    idx.write_text(t.replace('href="%s"' % site, 'href="%s"' % site.rstrip("/"), 1), encoding="utf-8")
    rc, log = check_neg()
    expect(rc == 1 and "canonical link forms" in log and "the website link" in log, "negative",
           "non-canonical website link not caught (exit %d)" % rc)
    summary.append("negative 5: website link without its trailing slash -> exit %d (%s)" % (rc, "caught" if rc == 1 else "MISSED"))
    idx.write_text(t, encoding="utf-8")

    # 6: the README block loses its data directory link (points at a file in v1/ instead)
    rd = neg / "README_header.md"
    r0 = rd.read_text(encoding="utf-8")
    rd.write_text(r0.replace("(%s)" % folder, "(%s/%s)" % (data, hosted[0]), 1), encoding="utf-8")
    rc, log = check_neg()
    expect(rc == 1 and "README block: the data link" in log, "negative",
           "missing README data link not caught (exit %d)" % rc)
    summary.append("negative 6: README block without the data directory link -> exit %d (%s)" % (rc, "caught" if rc == 1 else "MISSED"))
    rd.write_text(r0, encoding="utf-8")

    # 7: data_base_url and the url column of results/cells.tsv disagree
    (rc, log), _ = bad_config("v2.json", data_version="v2", data_base_url=folder + "v2")
    expect(rc == 2 and "not data_base_url + file name" in log, "negative",
           "a data_base_url that disagrees with results/cells.tsv is not rejected (exit %d)" % rc)
    summary.append("negative 7: data_base_url .../v2 against the v1 urls of results/cells.tsv -> exit %d (%s)"
                   % (rc, "rejected" if rc == 2 else "MISSED"))

    # 8: a page that no longer matches its sources
    page = neg / "docs" / "cells" / "111.html"
    p0 = page.read_text(encoding="utf-8")
    page.write_text(p0.replace("</main>", "<p>edited by hand</p>\n</main>", 1), encoding="utf-8")
    rc, log = check_neg()
    expect(rc == 1 and "out of date" in log and "111.html" in log, "negative", "stale page not caught (exit %d)" % rc)
    summary.append("negative 8: a hand-edited page -> exit %d (%s)" % (rc, "caught" if rc == 1 else "MISSED"))
    page.write_text(p0, encoding="utf-8")

    # 9: a null arXiv id builds cleanly: no link to the paper, never "arXiv:None"
    (rc, log), o = bad_config("noarxiv.json", arxiv_id=None)
    p4 = o / "docs" / "index.html"
    text = p4.read_text(encoding="utf-8") if p4.exists() else ""
    expect(rc == 0 and text and "arxiv.org/abs/%s" % base["arxiv_id"] not in text and "None" not in text, "check 9",
           "a null arxiv_id does not build cleanly (exit %d)" % rc)
    summary.append("check 9: arxiv_id null -> exit %d, no link to the paper, no 'None' (%s)" % (rc, "ok" if rc == 0 else "FAILED"))

    # 10: the other value of repo_public switches the data folder wording, and nothing else
    flip = not base["repo_public"]
    (rc, log), o = bad_config("public.json", repo_public=flip)
    repo = "https://github.com/%s/%s" % (owner, base["repo_name"])
    ok10 = rc == 0 and data_folder_ok(o / "data-folder", repo, flip)
    same_docs = all((work / "b" / "docs" / p.relative_to(o / "docs")).read_bytes() == p.read_bytes()
                    for p in (o / "docs").rglob("*") if p.is_file())
    expect(ok10 and same_docs, "check 10", "repo_public=%s: exit %d, data folder wording %s, docs/ %s"
           % (flip, rc, "ok" if ok10 else "WRONG", "unchanged" if same_docs else "CHANGED"))
    summary.append("check 10: repo_public %s -> exit %d, data folder wording %s, docs/ %s"
                   % (str(flip).lower(), rc, "ok" if ok10 else "WRONG", "unchanged" if same_docs else "CHANGED"))

    print("\n".join(summary))
    if failures:
        print("test builds kept in %s" % work)
        for f in failures:
            print("FAIL: " + f)
        return 1
    if args.keep:
        print("test builds kept in %s" % work)
    else:
        shutil.rmtree(work, ignore_errors=True)
    print("all config tests passed")
    return 0


if __name__ == "__main__":
    sys.exit(main())
