# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line, and also
# from the environment for the first two.
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

# TikZ diagrams -> SVG
TIKZDIR   ?= _tikz
TIKZOUT   ?= _static/svg
TIKZBUILD ?= $(patsubst %/,%,$(or $(TMPDIR),/tmp))/sphinx
TIKZSVG    = $(patsubst $(TIKZDIR)/%.tex,$(TIKZOUT)/%.svg,$(wildcard $(TIKZDIR)/*.tex))

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile 

clean:
	-rm -rf $(BUILDDIR)
	-rm -rf ref/generated

.PHONY: tikz
tikz: $(TIKZSVG)

$(TIKZOUT)/%.svg: $(TIKZDIR)/%.tex
		@mkdir -p $(TIKZBUILD) $(@D)
		@cd $(TIKZDIR) && latex -interaction=nonstopmode -halt-on-error \
				-output-directory=$(abspath $(TIKZBUILD)) -jobname=$* \
				"\PassOptionsToClass{dvisvgm}{standalone}\def\pgfsysdriver{pgfsys-dvisvgm.def}\input{$*.tex}" \
				>/dev/null || { echo "latex failed: see $(TIKZBUILD)/$*.log"; exit 1; }
		dvisvgm --no-fonts -v1 -o $@ $(TIKZBUILD)/$*.dvi
		@rm -f $(TIKZBUILD)/$*.*

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
