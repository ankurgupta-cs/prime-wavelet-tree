# Makefile -- tools, tests, and the per-cell checks of the prime wavelet tree
# repository. GNU make; a C++20 compiler (tested: GCC 16.1.0); bash for the scripts.
#
#   make                 build the tools into build/
#   make test            unit tests (tests/unit, fixtures in tests/fixtures)
#   make check           the prefix self-test: every encoder and decoder on the
#                        8,241 Carmichael numbers below 10^12 against frozen goldens
#   make test-release-tools   fixture checks of list_decode, pwt2_decode --emit-paths, paths_cmp
#
#   make fetch-CELL      download a hosted file into data/ and check its SHA-256
#   make verify-CELL     size + SHA-256 + structural decode of the cell's file
#   make certify-CELL    decode and compare every Carmichael number with the oracle
#   make rebuild-CELL    rebuild the file from the oracle; require the published bytes
#   make decode-CELL     write the decoded values to work/
#   CELL = 000 001 100 101 010 011 110 111 111x
#   make fetch-oracle | verify-oracle | rebuild-oracle | fetch-text | verify-text
#   make search-dstar    rerun the d* search as logged (about 33 h, 16 GB RAM)
#
#   make table           render the README table, the CELL manifests, SHA256SUMS
#                        and results/comparison.tex from results/cells.tsv
#   make pdf             build COMPARISON.pdf from results/comparison.tex
#   make site            build the website (docs/), the README links block and the
#                        data folder files (dist/data-folder/) from results/cells.tsv,
#                        the cell scripts and site/
#   make site-check      check that docs/ and the README block are up to date, check
#                        every page and link, and run the site tests
#   make table-check     every rendered file (table and site) up to date; runs site-check
#
# NATIVE=1 (default) compiles with -march=native. The published files were
# written by an x86-64 build with FMA3 (GCC 16.1.0, -O2 -march=native); the
# PWT2 encoder's model choice and the d* search use floating point, so byte
# identity of a rebuild is established for that build only. NATIVE=0 gives a
# portable build; every decoder is integer-only and portable either way.
# Variables: CXX, CXXFLAGS, NATIVE, BIN (build dir), DATA, WORK, ORC1, THREADS.

CXX      ?= g++
NATIVE   ?= 1
CXXFLAGS ?= -std=c++20 -O2 -Wall -Wextra
ifeq ($(NATIVE),1)
ARCHFLAGS = -march=native
else
ARCHFLAGS =
endif
ALL_CXXFLAGS = $(CXXFLAGS) $(ARCHFLAGS) -Isrc

ifeq ($(OS),Windows_NT)
EXE = .exe
else
EXE =
endif

B      ?= build
BASH   ?= bash
PYTHON ?= $(shell for p in python3 py python; do if $$p -c "import sys; sys.exit(sys.version_info < (3, 8))" >/dev/null 2>&1; then echo $$p; break; fi; done)
PDFLATEX ?= pdflatex

TOOLS = oracle_encode oracle_decode nstats_orc1 trie_paths dlist_from_paths \
        list_build list_decode cnd_encode cnd_decode pwt2_encode pwt2_decode \
        dsel_set2 paths_cmp factor_bench cofactor_bench flat_bench
TESTS = test_sieve test_oracle test_splitter test_cnd test_pwt test_pwt2 test_cnfactor
CELLS = 000 001 100 101 010 011 110 111 111x

HDRS      = $(wildcard src/*.hpp)
TOOL_BINS = $(TOOLS:%=$(B)/%$(EXE))
TEST_BINS = $(TESTS:%=$(B)/%$(EXE))

all: $(TOOL_BINS)

$(B):
	mkdir -p $(B)

$(B)/test_%$(EXE): tests/unit/test_%.cpp $(HDRS) | $(B)
	$(CXX) $(ALL_CXXFLAGS) $< -o $@

$(B)/%$(EXE): src/%.cpp $(HDRS) | $(B)
	$(CXX) $(ALL_CXXFLAGS) $< -o $@

# ------------------------------------------------------------------ tests
test: $(TEST_BINS)
	$(B)/test_sieve$(EXE)
	cd $(B) && ./test_oracle$(EXE)
	$(B)/test_splitter$(EXE) tests/fixtures
	$(B)/test_cnd$(EXE) tests/fixtures
	$(B)/test_pwt$(EXE) tests/fixtures
	$(B)/test_cnfactor$(EXE) tests/fixtures
	$(B)/test_pwt2$(EXE) tests/fixtures

check: $(TOOL_BINS)
	$(BASH) tests/prefix/selftest_prefix.sh --bin $(B)

test-release-tools: $(TOOL_BINS)
	$(BASH) tests/test_release_tools.sh $(B)

# ------------------------------------------------------------------ cells
fetch-%:
	$(BASH) scripts/cell.sh fetch $*

verify-%: $(TOOL_BINS)
	$(BASH) scripts/cell.sh verify $*

certify-%: $(TOOL_BINS)
	$(BASH) scripts/cell.sh certify $*

rebuild-%: $(TOOL_BINS)
	$(BASH) scripts/cell.sh rebuild $*

decode-%: $(TOOL_BINS)
	$(BASH) scripts/cell.sh decode $*

paths-%: $(TOOL_BINS)
	$(BASH) scripts/cell.sh paths $*

search-dstar: $(TOOL_BINS)
	ROUTE=search $(BASH) cells/111x-pwt-dstar/build.sh

verify-hosted: $(TOOL_BINS)
	for c in 001 101 011 111 111x oracle; do $(BASH) scripts/cell.sh verify $$c || exit 1; done

# ------------------------------------------------------------------ table and PDF
table:
	@test -n "$(PYTHON)" || { echo "Python 3.8+ not found (set PYTHON=...)"; exit 2; }
	$(PYTHON) results/render.py

table-check: site-check
	@test -n "$(PYTHON)" || { echo "Python 3.8+ not found (set PYTHON=...)"; exit 2; }
	$(PYTHON) results/render.py --check

site:
	@test -n "$(PYTHON)" || { echo "Python 3.8+ not found (set PYTHON=...)"; exit 2; }
	$(PYTHON) site/build_site.py

site-check:
	@test -n "$(PYTHON)" || { echo "Python 3.8+ not found (set PYTHON=...)"; exit 2; }
	$(PYTHON) site/build_site.py --check
	$(PYTHON) site/test_site.py

pdf: results/comparison.tex
	cd results && $(PDFLATEX) -interaction=nonstopmode -halt-on-error -jobname=COMPARISON comparison.tex >/dev/null && \
	  $(PDFLATEX) -interaction=nonstopmode -halt-on-error -jobname=COMPARISON comparison.tex >/dev/null
	mv results/COMPARISON.pdf COMPARISON.pdf
	rm -f results/COMPARISON.aux results/COMPARISON.log results/COMPARISON.out

info:
	@echo "CXX       $(CXX)"; $(CXX) --version | head -n 1
	@echo "flags     $(ALL_CXXFLAGS)"
	@echo "python    $(PYTHON)"

clean:
	rm -rf $(B)

.PHONY: all test check test-release-tools search-dstar verify-hosted table table-check site site-check pdf info clean
.SECONDARY:
