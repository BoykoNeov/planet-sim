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
import webbrowser
from dataclasses import dataclass
from pathlib import Path

_HERE = Path(__file__).resolve().parent
_REPO_ROOT = _HERE.parent
_NOTEBOOK = _HERE / "planet.ipynb"


@dataclass(frozen=True)
class Demo:
    """One catalogue entry — a runnable ``planet.<module>`` exposing ``main()``."""

    key: str                 # the short name the user types
    module: str              # the importable module whose main() we call
    title: str               # one-line headline for the menu
    blurb: str               # a sentence of what it shows
    extras: tuple[str, ...]  # the pip extras its figure needs ("viz" / "webviz")
    artifact: str            # repo-relative primary artifact (offered to open afterwards)
    section: str             # menu grouping header
    sim: bool = False        # True ⇒ runs a multi-second fluid simulation (warn first)


# The single source of truth. The interactive menu, the CLI dispatch, and the catalogue
# test all read this one list — add a demo here and it shows up in every surface at once.
DEMOS: tuple[Demo, ...] = (
    Demo("snowball", "planet.demo_snowball", "Snowball-Earth hysteresis",
         "one knob (the solar constant), two stable climates, a catastrophic freeze",
         ("viz",), "docs/figures/planet-snowball.png", "Climate — energy balance"),
    Demo("biomes", "planet.demo_biomes", "Climate → biome map",
         "the Whittaker (temperature, rainfall) classifier; warming migrates the bands poleward",
         ("viz",), "docs/figures/planet-biomes.png", "Climate — energy balance"),
    Demo("exoplanet", "planet.demo_exoplanet", "Exoplanet knobs",
         "a redder star and a bigger planet reshape the climate and the ice line",
         ("viz",), "docs/figures/planet-exoplanet.png", "Climate — energy balance"),
    Demo("obliquity", "planet.demo_obliquity", "Axial tilt (obliquity)",
         "how the planet's tilt reshapes the pole-to-equator sunlight",
         ("viz",), "docs/figures/planet-obliquity.png", "Climate — energy balance"),
    Demo("shallowwater", "planet.demo_shallowwater", "Rotating shallow-water atmosphere",
         "geostrophic adjustment on the sphere — the circulation engine on its own",
         ("viz",), "docs/figures/planet-shallowwater.png",
         "Circulation — shallow-water (runs a short sim)", sim=True),
    Demo("coupler", "planet.demo_coupler", "EBM → circulation coupler",
         "an emergent jet grows from the pole-to-equator temperature gradient (one-way coupling)",
         ("viz", "webviz"), "docs/figures/planet-coupler.png",
         "Circulation — shallow-water (runs a short sim)", sim=True),
    Demo("eddy_life", "planet.demo_eddy_life", "Eddy life cycle — GIF",
         "the emergent eddy stirring the temperature, animated as a two-panel GIF",
         ("viz",), "docs/figures/planet-eddy-life.gif",
         "Circulation — shallow-water (runs a short sim)", sim=True),
    Demo("eddy_globe", "planet.demo_eddy_globe", "Eddy life cycle — globe",
         "the same eddy life cycle, animated on the interactive globe",
         ("webviz",), "docs/figures/planet-eddy-globe.html",
         "Circulation — shallow-water (runs a short sim)", sim=True),
    Demo("map", "planet.planetmap", "Interactive biome-map globe",
         "the present-day globe — rotate / zoom / hover (the live sliders run in the notebook)",
         ("webviz",), "docs/figures/planet-map.html", "Interactive globes"),
)

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


def _launch_notebook() -> None:
    """Open the teaching notebook in JupyterLab — *we* open the browser, so it never makes you
    copy a URL.

    JupyterLab's own auto-open writes a temporary redirect ``.html`` and points the browser at
    that file; on Windows the browser often can't read it (``ERR_ACCESS_DENIED``). So we launch
    with ``--no-browser``, read the ``http://localhost:PORT/lab?token=…`` line off its log, and
    open the notebook URL ourselves — the one mechanism that actually broke is removed entirely.
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
    # Serve from the repo root so the notebook opens at a stable /lab/tree/ path.
    proc = subprocess.Popen(
        [jupyter, "lab", "--no-browser"],
        cwd=str(_REPO_ROOT),
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
        text=True, bufsize=1,
    )
    opened = False
    try:
        assert proc.stdout is not None
        for line in proc.stdout:
            sys.stdout.write(line)
            sys.stdout.flush()
            if not opened:
                nb_url = _notebook_url_from_log(line)
                if nb_url:
                    print(f"\n→ Opening the notebook in your browser:\n  {nb_url}\n"
                          "  (if it doesn't open, paste that URL into your browser)\n")
                    webbrowser.open(nb_url)
                    opened = True
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
