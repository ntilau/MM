# ============================================================================
# Build configuration for Master's thesis manuscript and defense presentation
# ============================================================================

# --- Manuscript Build Configuration ---
# MS_TEX: LaTeX engine (default: pdflatex)
# MS_TEXFLAGS: Flags passed to pdflatex
#   -interaction=nonstopmode: non-interactive mode (skip prompts)
#   -file-line-error: show errors in file:line format
# MS_BIB: Bibliography processor (optional, empty by default)
# MS_FILE: Main LaTeX source file (default: thesis)
MS_TEX      = pdflatex
MS_TEXFLAGS ?= 
#-interaction=nonstopmode -file-line-error
MS_BIB      ?= 
MS_FILE     ?= paper


# LaTeX auxiliary files to remove during clean
# Includes: cross-references, bibliography, font database, logs, etc.
CLEAN_EXTS = *.aux *-blx.bib *.bcf *.blg *.bbl *.brf *.fdb_latexmk *.fls *.loa *.lof *.log *.lot *.lpr *.nav *.out *.run.xml *.snm *.toc *.vrb *.synctex.gz

.SUFFIXES: .aux .pdf .tex
.PHONY: all paper clean distclean help

# Build both manuscript and presentation
all: paper

# Build manuscript PDF
# Process: LaTeX → (Bibliography if configured) → LaTeX → LaTeX
# Three LaTeX passes ensure cross-references and citations are resolved
paper:
	cd manuscript && $(MS_TEX) $(MS_TEXFLAGS) $(MS_FILE).tex
	@if [ -n "$(MS_BIB)" ]; then cd manuscript && $(MS_BIB) $(MS_FILE); fi
	cd manuscript && $(MS_TEX) $(MS_TEXFLAGS) $(MS_FILE).tex
	cd manuscript && $(MS_TEX) $(MS_TEXFLAGS) $(MS_FILE).tex

# Remove temporary LaTeX files (keep PDFs)
clean:
	cd manuscript && rm -f $(CLEAN_EXTS)

# Remove temporary files AND generated PDFs
distclean: clean
	rm -f manuscript/$(MS_FILE).pdf presentation/$(PR_FILE).pdf

# Display help message with available targets and variables
help:
	@echo "========== LaTeX Build System =========="
	@echo ""
	@echo "TARGETS:"
	@echo "  make / make all      - Build BOTH manuscript and presentation"
	@echo "  make manuscript      - Build ONLY manuscript/$(MS_FILE).pdf"
	@echo "  make presentation    - Build ONLY presentation/$(PR_FILE).pdf"
	@echo "  make clean           - Remove temporary LaTeX files (keeps PDFs)"
	@echo "  make distclean       - Remove temporary files AND PDFs"
	@echo "  make help            - Show this help message"
	@echo ""
	@echo "CONFIGURABLE VARIABLES (override via: make VAR=value):"
	@echo "  Manuscript:"
	@echo "    MS_FILE=$(MS_FILE)          - Main source file (no .tex extension)"
	@echo "    MS_TEX=$(MS_TEX)            - LaTeX engine"
	@echo "    MS_TEXFLAGS                - Compiler flags"
	@echo "    MS_BIB=$(MS_BIB)            - Bibliography processor (biber, bibtex, or empty)"
	@echo ""
	@echo "  Presentation:"
	@echo "    PR_FILE=$(PR_FILE)          - Main source file (no .tex extension)"
	@echo "    PR_TEX=$(PR_TEX)            - LaTeX engine"
	@echo "    PR_TEXFLAGS                - Compiler flags"
	@echo "    PR_BIB=$(PR_BIB)            - Bibliography processor (biber, bibtex, or empty)"
	@echo ""
	@echo "EXAMPLES:"
	@echo "  make                             # Build both documents"
	@echo "  make manuscript                  # Build thesis only"
	@echo "  make clean                       # Clean all temporary files"
	@echo "  make MS_BIB=biber manuscript     # Build with biber bibliography"
	@echo "  make help                        # Show this help message"

