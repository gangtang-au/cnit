"""Build Sphinx documentation using notebooks from the top-level notebooks folder."""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path


ROOT_DIR = Path(__file__).resolve().parents[1]
DOCS_DIR = ROOT_DIR / "docs"
NOTEBOOKS_DIR = ROOT_DIR / "notebooks"


def _validate_notebook_references(notebooks: list[Path]) -> None:
    """Ensure all notebooks are referenced and all references have a source."""
    available = {notebook.stem for notebook in notebooks}
    referenced = set()
    for rst_file in DOCS_DIR.rglob("*.rst"):
        content = rst_file.read_text(encoding="utf-8")
        notebook_refs = re.finditer(
            r"(?m)^\s{3,}(?:[^<>\n]*<[^<>]*?(CNit_[\w-]+)>|(CNit_[\w-]+))\s*$",
            content,
        )
        referenced.update(match.group(1) or match.group(2) for match in notebook_refs)

    missing = referenced - available
    unused = available - referenced
    if missing or unused:
        details = []
        if missing:
            details.append(f"missing notebook sources: {', '.join(sorted(missing))}")
        if unused:
            details.append(f"unreferenced notebooks: {', '.join(sorted(unused))}")
        raise SystemExit("Notebook documentation check failed: " + "; ".join(details))


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DOCS_DIR / "_build" / "html",
        help="HTML output directory (default: docs/_build/html)",
    )
    args = parser.parse_args()

    notebooks = sorted(NOTEBOOKS_DIR.glob("*.ipynb"))
    if not notebooks:
        raise SystemExit(f"No Jupyter notebooks found in {NOTEBOOKS_DIR}")
    _validate_notebook_references(notebooks)

    args.output_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="cnit-docs-") as temporary_dir:
        staged_docs = Path(temporary_dir) / "docs"
        shutil.copytree(
            DOCS_DIR,
            staged_docs,
            ignore=shutil.ignore_patterns(
                "_build", ".doctrees", ".ipynb_checkpoints", "__pycache__"
            ),
        )
        for notebook in notebooks:
            shutil.copy2(notebook, staged_docs / notebook.name)

        command = [
            sys.executable,
            "-m",
            "sphinx",
            "-b",
            "html",
            "-c",
            str(DOCS_DIR),
            str(staged_docs),
            str(args.output_dir.resolve()),
        ]
        result = subprocess.run(command, cwd=ROOT_DIR, check=False)
        if result.returncode == 0:
            # nbsphinx copies the source notebooks into the HTML output. The
            # workflow pages link to the canonical notebooks on GitHub instead.
            for notebook in notebooks:
                (args.output_dir.resolve() / notebook.name).unlink(missing_ok=True)
        return result.returncode


if __name__ == "__main__":
    raise SystemExit(main())
