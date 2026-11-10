
#!/usr/bin/env bash
set -Eeuo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
PAPER_DIR="$ROOT_DIR/papers/01_maxwell_helmholtz"
PDF_SOURCE="$PAPER_DIR/main.pdf"
PDF_DEST="$ROOT_DIR/docs/papers/maxwell-helmholtz.pdf"

echo "======================================="
echo " WaveCore Research - Local CI"
echo "======================================="

echo "[1/6] Checking dependencies..."

for tool in latexmk lualatex biber python3; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "ERROR: Missing dependency: $tool" >&2
        exit 1
    fi
done

if ! python3 -m mkdocs --version >/dev/null 2>&1; then
    echo "ERROR: MkDocs is not installed." >&2
    echo "Activate .venv and install requirements-docs.txt" >&2
    exit 1
fi

for file in \
    "$PAPER_DIR/main.tex" \
    "$ROOT_DIR/mkdocs.yml" \
    "$ROOT_DIR/docs/index.md" \
    "$ROOT_DIR/docs/publications.md" \
    "$ROOT_DIR/docs/about.md"
do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: Missing file: $file" >&2
        exit 1
    fi
done

echo "[2/6] Building PDF with LuaLaTeX..."

(
    cd "$PAPER_DIR"

    latexmk -lualatex \
        -interaction=nonstopmode \
        -file-line-error \
        -halt-on-error \
        main.tex
)

echo "[3/6] Checking generated PDF..."

if [[ ! -s "$PDF_SOURCE" ]]; then
    echo "ERROR: PDF is missing or empty." >&2
    exit 1
fi

if command -v pdfinfo >/dev/null 2>&1; then
    pdfinfo "$PDF_SOURCE" | grep -E 'Pages:|Page size:'
fi

echo "[4/6] Preparing website files..."

mkdir -p "$(dirname "$PDF_DEST")"
cp "$PDF_SOURCE" "$PDF_DEST"

echo "[5/6] Building MkDocs website..."

cd "$ROOT_DIR"
python3 -m mkdocs build --strict

echo "[6/6] Checking generated website..."

if [[ ! -s "$ROOT_DIR/site/index.html" ]]; then
    echo "ERROR: Website index.html was not generated." >&2
    exit 1
fi

echo
echo "======================================="
echo " LOCAL CI PASSED"
echo "======================================="
echo "PDF:  $PDF_DEST"
echo "Site: $ROOT_DIR/site/index.html"
