"""Execution smoke-test for the planet teaching notebook (Planet plan §9 / ADR 0002).

Per ADR 0002 the notebook (`planet.ipynb`) is a **reach** layer, not a correctness
one: its physics is already validated behind the `ebm` / `albedo` / `precip` /
`biomes` triads, and it is a *thin skin* — each compute cell calls those module
functions directly. So this test asks the one thing that matters for a notebook:
**does it execute clean, top to bottom**, no cell raising. It is *not* a physics
check.

Why the discipline matters here (and why the load-bearing compute lives in plain
cells, not `interact` callbacks): `ipywidgets.interact` runs its callback inside an
`Output` context manager that *captures* exceptions and paints them as output
instead of re-raising — so a broken module call inside an interact callback would
leave the cell "successful" and this test green. The notebook therefore puts every
validated call in a direct cell (`allow_errors=False` below makes any such cell's
exception fail this test); the interact cells are sugar on top. (Same rule, and the
same rationale, as `projects/chip/tests/test_chip_notebook.py` and
`projects/steel/tests/test_steel_notebook.py`.)

Three gates keep a headless / clean checkout *skipping* rather than *erroring*,
like the viz tests: the optional `[notebook]` execution stack, matplotlib
(`[viz]`), and a **registered Jupyter kernelspec** — separate from merely having
`ipykernel` importable, so it is checked explicitly.

**Why a subprocess.** The notebook is executed in a *fresh* child process rather than
in-process. On Windows the kernel client polls zmq over asyncio, and if this process
already has a running/cached event loop (left by another test in the suite), nbclient
takes a thread-runner path that can deadlock pyzmq on the Proactor loop. A clean child
process always gets the fast path, and `subprocess.run(timeout=…)` wall-clocks it so a
pathological hang fails *this test* fast instead of wedging the whole suite. The child
entry point is the ``__main__`` block at the bottom of this file.
"""
import subprocess
import sys
from pathlib import Path

import pytest

NOTEBOOK = Path(__file__).resolve().parents[1] / "planet.ipynb"
REPO_ROOT = Path(__file__).resolve().parents[3]


@pytest.mark.slow
def test_planet_notebook_executes_clean():
    # @slow (ADR 0003): spawns a fresh kernel in a child process — deselected from the
    # fast inner loop (`pytest -m "not slow"`), always run in the full commit gate.
    # Gate on the optional execution stack (the [notebook] extra) + the render dep.
    pytest.importorskip("nbformat")
    pytest.importorskip("nbclient")
    pytest.importorskip("ipykernel")
    pytest.importorskip("ipywidgets")
    pytest.importorskip("matplotlib")

    # importorskip checks the *packages*; executing also needs a registered Jupyter
    # kernelspec (NOT guaranteed by `pip install ipykernel` alone). Skip — don't
    # error — if none is available, mirroring the importorskip philosophy.
    from jupyter_client.kernelspec import KernelSpecManager

    specs = KernelSpecManager().find_kernel_specs()
    kernel = "python3" if "python3" in specs else next(iter(specs), None)
    if kernel is None:
        pytest.skip("no registered Jupyter kernelspec to execute the notebook")

    assert NOTEBOOK.exists(), f"missing teaching notebook: {NOTEBOOK}"

    # Run the executor (this file's __main__) in a clean child process. It exits
    # non-zero with a traceback on any cell error (allow_errors=False); a clean run
    # exits 0. errors="replace" keeps a unicode traceback (°C, →, φ) from masking the
    # real failure behind a decode error on a legacy codepage.
    proc = subprocess.run(
        [sys.executable, str(Path(__file__).resolve()), str(NOTEBOOK), str(REPO_ROOT), kernel],
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=300,
    )
    assert proc.returncode == 0, (
        "planet.ipynb did not execute clean:\n"
        f"--- stdout ---\n{proc.stdout[-2000:]}\n--- stderr ---\n{proc.stderr[-3000:]}"
    )


if __name__ == "__main__":
    # Child entry: execute the notebook headless; raise (→ non-zero exit) on any cell
    # error. Invoked as `python <thisfile> <notebook> <repo_root> <kernel_name>`.
    import nbformat as nbf
    from nbclient import NotebookClient

    nb_path, repo_root, kernel_name = sys.argv[1], sys.argv[2], sys.argv[3]
    notebook = nbf.read(nb_path, as_version=4)
    NotebookClient(
        notebook,
        timeout=180,
        kernel_name=kernel_name,
        resources={"metadata": {"path": repo_root}},
    ).execute()
    print("planet.ipynb executed clean")
