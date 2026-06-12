"""``python -m planet`` — the one front door to the simulator.

Every leg of planet-sim (the EBM/Snowball climate demos, the shallow-water circulation,
the emergent jet & eddy animations, the interactive globes, and the teaching notebook) is
reachable from here, so a newcomer following the README has *one* command to remember and a
menu to explore from — not a scavenger hunt across nine ``python -m planet.demo_*`` modules.

    python -m planet                # interactive menu (pick a number or a name)
    python -m planet snowball       # run one demo straight off
    python -m planet notebook       # open the teaching notebook in JupyterLab (opens the browser)
    python -m planet globes         # just open a saved interactive globe (no compute)
    python -m planet list           # print the catalogue and exit
    python -m planet all            # run every demo and bank every figure

Each demo prints its validation table and banks a figure under ``docs/figures/`` (and
``outputs/``); the optional render stacks are opt-in (``pip install -e .[viz]`` /
``.[webviz]`` / ``.[notebook]``) — a demo whose stack is missing still prints its physics
summary and tells you which extra to install rather than crashing.
"""
from __future__ import annotations

import importlib
import re
import shutil
import subprocess
import sys
import time
import webbrowser
from pathlib import Path

from planet.catalog import DEMOS, Demo, _REPO_ROOT  # the single source of truth (re-exported)

_HERE = Path(__file__).resolve().parent
_NOTEBOOK = _HERE / "planet.ipynb"
_SERVER_LOG = _REPO_ROOT / "outputs" / "jupyter-server.log"  # gitignored; the server logs here

_BY_KEY = {d.key: d for d in DEMOS}

# The standalone globes a fresh clone already ships (so "just open one" needs no compute).
_BANKED_GLOBES: tuple[tuple[str, str], ...] = (
    ("the biome-map globe", "docs/figures/planet-map.html"),
    ("the coupler jet-over-temperature globe", "docs/figures/planet-coupler-map.html"),
    ("the eddy life-cycle globe animation", "docs/figures/planet-eddy-globe.html"),
)


def _utf8_stdout() -> None:
    """Make °C / ₂ / → printable on a legacy Windows code page (mirrors each demo's main())."""
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")


def _open_in_browser(path: Path, label: str = "") -> bool:
    """Open a saved figure / HTML globe in the default browser (as a proper file:// URL)."""
    if not path.exists():
        print(f"  (not found: {path.relative_to(_REPO_ROOT)} — run the demo that banks it first)")
        return False
    print(f"  opening {label or path.name} → {path.relative_to(_REPO_ROOT)}")
    webbrowser.open(path.as_uri())
    return True


def _run_demo(demo: Demo, offer_open: bool) -> None:
    """Dispatch one demo by calling the module's existing ``main()`` (graceful if its viz is absent)."""
    print(f"\n── {demo.title} " + "─" * max(0, 60 - len(demo.title)))
    print(f"   {demo.blurb}")
    if demo.sim:
        print("   (this runs a live shallow-water simulation — give it a minute)")
    print()
    try:
        importlib.import_module(demo.module).main()
    except KeyboardInterrupt:
        print("\n   interrupted.")
        return
    except Exception as exc:  # one demo blowing up must not take the menu down with it
        print(f"\n   {demo.title} failed: {type(exc).__name__}: {exc}")
        return

    artifact = _REPO_ROOT / demo.artifact
    if offer_open and artifact.exists():
        if _prompt_yes_no(f"\n   Open {artifact.name} now?", default=True):
            _open_in_browser(artifact)


# The server prints lines like `http://localhost:8889/lab?token=abc…`; we turn that into the
# /lab/tree/ URL that opens the notebook file itself. Pulled out as a pure helper so a test can
# pin the construction without launching a server.
_SERVER_URL_RE = re.compile(r"(https?://(?:localhost|127\.0\.0\.1):\d+)/lab\?token=([\w-]+)")


def _notebook_url_from_log(line: str) -> str | None:
    """If `line` is a Jupyter 'server is running at' line, return the URL that opens the notebook."""
    m = _SERVER_URL_RE.search(line)
    if not m:
        return None
    base, token = m.group(1), m.group(2)
    return f"{base}/lab/tree/planet/planet.ipynb?token={token}"


def _wait_for_notebook_url(log_path: Path, proc: subprocess.Popen, timeout: float = 90.0) -> str | None:
    """Tail the server's log file until it announces its URL — or it dies / we time out.

    We read the log *file* rather than the live stdout pipe on purpose (see ``_launch_notebook``):
    iterate with ``readline()`` and short-sleep on EOF, because a plain ``for line in f`` would stop
    at the first end-of-file — which we reach before JupyterLab has finished booting and printed its
    URL. Bail early if the server process exits (a startup crash) so we don't wait the full timeout.
    """
    deadline = time.monotonic() + timeout
    with open(log_path, "r", encoding="utf-8", errors="replace") as log:
        while time.monotonic() < deadline:
            line = log.readline()
            if not line:                              # caught up to the writer — wait for more
                if proc.poll() is not None:           # …unless the server already exited (crash)
                    return None
                time.sleep(0.2)
                continue
            sys.stdout.write(line)                    # echo the boot banner (bounded: we stop at the URL)
            sys.stdout.flush()
            nb_url = _notebook_url_from_log(line)
            if nb_url:
                return nb_url
    return None


def _launch_notebook() -> None:
    """Open the teaching notebook in JupyterLab — *we* open the browser, so it never makes you
    copy a URL.

    Two Windows traps are designed out here:

    * **The redirect file.** JupyterLab's own auto-open writes a temporary redirect ``.html`` and
      points the browser at it; on Windows the browser often can't read that file
      (``ERR_ACCESS_DENIED``). So we launch ``--no-browser`` and open the ``/lab/tree/`` URL ourselves.
    * **The frozen server.** If the server's output goes through a pipe we echo to the console, a
      *paused* console — Windows QuickEdit / selecting text in the window pauses its output — blocks
      our echo, stops us draining the pipe, fills the OS pipe buffer, and then JupyterLab blocks on
      its own log writes and goes unreachable (the browser shows "error connecting to server" and
      "File Save Error … Failed to fetch"). So the server logs to a **file** instead: a file never
      back-pressures, so no console habit can starve it. We detect the URL by tailing that file.
    """
    if not _NOTEBOOK.exists():
        print(f"Notebook not found: {_NOTEBOOK}")
        return
    jupyter = shutil.which("jupyter") or shutil.which("jupyter-lab")
    if jupyter is None:
        print("JupyterLab isn't installed. Install the notebook stack:\n"
              "  pip install -e \".[viz,notebook]\"")
        return

    print(f"Launching JupyterLab on {_NOTEBOOK.relative_to(_REPO_ROOT)} …")
    print("(leave this running; press Ctrl-C here to stop the server)\n")
    _SERVER_LOG.parent.mkdir(parents=True, exist_ok=True)
    # Truncate per launch ('w') so a stale URL/token from a previous run can't be picked up.
    log = open(_SERVER_LOG, "w", encoding="utf-8")
    try:
        # Serve from the repo root so the notebook opens at a stable /lab/tree/ path. The child gets
        # its own inherited handle to the log, so we can close ours straight after it starts.
        proc = subprocess.Popen(
            [jupyter, "lab", "--no-browser"],
            cwd=str(_REPO_ROOT),
            stdout=log, stderr=subprocess.STDOUT,
            text=True,
        )
    finally:
        log.close()

    log_rel = _SERVER_LOG.relative_to(_REPO_ROOT)
    nb_url = _wait_for_notebook_url(_SERVER_LOG, proc)
    if nb_url:
        print(f"\n→ Opening the notebook in your browser:\n  {nb_url}\n"
              "  (if it doesn't open, paste that URL into your browser)\n")
        webbrowser.open(nb_url)
    elif proc.poll() is not None:
        print(f"\nJupyterLab exited before announcing a URL — see {log_rel} for why.")
        return
    else:
        print(f"\nCouldn't detect the server URL yet — open {log_rel}, or run "
              "`jupyter lab list`, to get it.\n")

    print(f"JupyterLab is running (its log → {log_rel}). Press Ctrl-C here to stop the server.\n")
    try:
        proc.wait()
    except KeyboardInterrupt:
        print("\nStopping JupyterLab …")
        proc.terminate()
        try:
            proc.wait(timeout=10)
        except subprocess.TimeoutExpired:
            proc.kill()


def _open_banked_globe(choice: str | None = None) -> None:
    """Open one of the committed interactive globes — no compute, just the saved HTML."""
    if choice is None:
        print("\nWhich interactive globe?")
        for i, (label, rel) in enumerate(_BANKED_GLOBES, 1):
            print(f"  {i}  {label}")
        choice = _read("  pick a number (Enter = 1): ") or "1"
    try:
        idx = int(choice) - 1
        label, rel = _BANKED_GLOBES[idx]
    except (ValueError, IndexError):
        print(f"  no such globe: {choice!r}")
        return
    _open_in_browser(_REPO_ROOT / rel, label)


def _build_and_open_site() -> None:
    """(Re)generate the landing page from the catalogue and open it in the browser."""
    from planet import site                      # local import: keeps the launcher start-up lean
    path = site.write_site()
    print(f"  landing page generated → {path.relative_to(_REPO_ROOT)}")
    _open_in_browser(path, "the planet-sim landing page")


def _run_all() -> None:
    """Run every demo in catalogue order (banks every figure). Never auto-opens — that'd be 9 tabs."""
    for demo in DEMOS:
        _run_demo(demo, offer_open=False)
    print("\nAll demos run. Figures are under docs/figures/ (and outputs/).")


# --- the interactive menu --------------------------------------------------- #
def _extras_tag(demo: Demo) -> str:
    return f"[{','.join(demo.extras)}]" if demo.extras else ""


def print_catalog() -> None:
    """Print the grouped catalogue — the menu body, reused by ``list`` and the non-TTY fallback."""
    print("\n  planet-sim — a planetary climate simulator\n")
    last_section = None
    for i, demo in enumerate(DEMOS, 1):
        if demo.section != last_section:
            print(f"  {demo.section}:")
            last_section = demo.section
        print(f"    {i:>2}  {demo.key:<13} {demo.title}  {_extras_tag(demo)}")
        print(f"        {demo.blurb}")
    print("\n  Open a saved result (no compute):")
    print("     g  globes        open a banked interactive globe (biome / coupler / eddy)")
    print("     s  site          build & open the landing page (links to every demo & globe)")
    print("     n  notebook      open the teaching notebook in JupyterLab          [notebook]")
    print("\n     a  all           run every demo and bank every figure")
    print("     q  quit\n")


def _menu_loop() -> None:
    print_catalog()
    while True:
        try:
            raw = _read("planet> ")
        except KeyboardInterrupt:             # Ctrl-C at the prompt = quit cleanly, not a traceback
            print("\n(quit)")
            return
        if raw is None:                       # EOF (Ctrl-D / closed stdin)
            print()
            return
        cmd = raw.strip().lower()
        if cmd in ("", "q", "quit", "exit"):
            return
        if cmd in ("l", "list", "?", "help", "h"):
            print_catalog()
        elif cmd in ("a", "all"):
            _run_all()
        elif cmd in ("n", "notebook"):
            _launch_notebook()
        elif cmd in ("g", "globe", "globes"):
            _open_banked_globe()
        elif cmd in ("s", "site", "page", "web"):
            _build_and_open_site()
        elif cmd.isdigit() and 1 <= int(cmd) <= len(DEMOS):
            _run_demo(DEMOS[int(cmd) - 1], offer_open=True)
        elif cmd in _BY_KEY:
            _run_demo(_BY_KEY[cmd], offer_open=True)
        else:
            print(f"  not a choice: {raw!r}  (type a number, a name, or 'q')")


# --- small IO helpers (so tests can drive a non-interactive run) ------------ #
def _read(prompt: str) -> str | None:
    """input() that returns None on EOF instead of raising — keeps the menu loop simple."""
    try:
        return input(prompt)
    except EOFError:
        return None


def _prompt_yes_no(prompt: str, default: bool = True) -> bool:
    suffix = " [Y/n] " if default else " [y/N] "
    ans = _read(prompt + suffix)
    if ans is None or ans.strip() == "":
        return default
    return ans.strip().lower().startswith("y")


def _dispatch_one(token: str) -> int:
    """Run a single CLI token (a demo key, or one of the verbs). Returns a process exit code."""
    tok = token.lower()
    if tok in ("list", "--list"):
        print_catalog()
    elif tok in ("all",):
        _run_all()
    elif tok in ("notebook", "nb"):
        _launch_notebook()
    elif tok in ("globes", "globe"):
        _open_banked_globe()
    elif tok in ("site", "page", "web"):
        _build_and_open_site()
    elif tok in _BY_KEY:
        _run_demo(_BY_KEY[tok], offer_open=sys.stdout.isatty())
    else:
        print(f"Unknown demo: {token!r}\nValid names: {', '.join(_BY_KEY)}")
        print("Run `python -m planet list` for the full catalogue.")
        return 2
    return 0


def main(argv: list[str] | None = None) -> int:
    """Entry point. No args ⇒ interactive menu on a TTY, else the catalogue + usage."""
    _utf8_stdout()
    args = sys.argv[1:] if argv is None else argv
    if args:
        if args[0] in ("-h", "--help"):
            print(__doc__)
            return 0
        return _dispatch_one(args[0])

    if sys.stdin.isatty():
        _menu_loop()
    else:
        # Piped / non-interactive with no args: print the catalogue rather than hang on input().
        print_catalog()
        print("Run e.g. `python -m planet snowball`, or `python -m planet` in a terminal for the menu.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
