"""Paths and helpers shared by the Nastran workflow scripts.

Layout: nastran/simple_beams/<example>/mesh<N>/ holds the .bdf, the Nastran outputs (.op2, .f06, ...), the extracted response and meta.json of each case.
"""
import json
from pathlib import Path

NASTRAN_DIR = Path(__file__).resolve().parent
META_FILE = "meta.json"


def case_dir(example, mesh):
    return NASTRAN_DIR / "simple_beams" / example / f"mesh{mesh}"


def meta_path(directory, job=None):
    """Metadata file of a case: meta.json, or <job>.meta.json for a job with its own metadata (e.g. one per reference axis)."""
    return Path(directory) / (META_FILE if job is None else f"{job}.{META_FILE}")


def write_meta(directory, job=None, **data):
    meta_path(directory, job).write_text(json.dumps(data, indent=2))


def read_meta(directory, job=None):
    """Metadata of the job if it has its own file, otherwise the one of the case (meta.json)."""
    path = meta_path(directory, job)
    return json.loads((path if path.exists() else meta_path(directory)).read_text())


def write_bdf_if_changed(model, path, size=8):
    """Writes the BDF only if its content differs from the existing file, so an unchanged model keeps its Nastran results. When the content of an existing file changes, the outputs of the previous run (op2, f06, extracted response, ...) are stale and are removed (not when the file is new). Returns True if written. size: field width of the cards (8 or 16)."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.stem + ".tmp.bdf")
    model.write_bdf(str(tmp), size=size)
    if path.exists() and path.read_bytes() == tmp.read_bytes():
        tmp.unlink()
        return False
    existed = path.exists()
    tmp.replace(path)
    # The outputs of a previous deck are stale when it changes. For a deck that did not exist, any outputs with its name (e.g. the results of the same job run elsewhere, and brought here before the deck was generated) are kept
    if existed:
        for stale in path.parent.glob(path.stem + ".*"):
            if stale != path:
                stale.unlink()
    return True
