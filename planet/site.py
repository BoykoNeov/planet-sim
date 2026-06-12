"""Generate ``docs/index.html`` — the clickable landing page for planet-sim.

A self-contained gallery built **from** the demo catalogue (:data:`planet.catalog.DEMOS`), so it
can never drift: add a demo to the catalogue, run ``python -m planet site`` to regenerate, and a
golden test (``planet/tests/test_site.py``) fails the build if the committed page is out of date.

The page links to the three interactive globes and the still figures under ``docs/figures/`` with
**relative** hrefs (so it works both opened straight off disk and served by GitHub Pages from the
``/docs`` folder), and to the notebook / build plan / ADRs on GitHub (so the Markdown renders). It
is intentionally dependency-free — inline CSS, no CDN or web font — and **deterministic** (no build
timestamp, links emitted straight from the catalogue, never branched on which figures are on disk)
so the golden comparison is stable.

Run it (writes the page, prints the path)::

    python -m planet site
"""
from __future__ import annotations

import html
import sys
from pathlib import Path

from planet.catalog import DEMOS, Demo, _REPO_ROOT, globe_href

SITE_PATH = _REPO_ROOT / "docs" / "index.html"

_REPO_URL = "https://github.com/BoykoNeov/planet-sim"
_BLOB = f"{_REPO_URL}/blob/main"      # absolute, for the Markdown docs (render on GitHub)


def _rel_href(repo_rel: str) -> str:
    """A docs-relative href to an artifact (the page lives in docs/, the artifacts in docs/figures/)."""
    prefix = "docs/"
    assert repo_rel.startswith(prefix), repo_rel
    return Path(repo_rel[len(prefix):]).as_posix()


def _esc(text: str) -> str:
    return html.escape(text, quote=True)


# --- page pieces (each returns a chunk of HTML; all dynamic text is escaped) ------------------ #
_CSS = """\
:root { color-scheme: dark; }
* { box-sizing: border-box; }
body { margin: 0; font: 16px/1.55 system-ui, -apple-system, Segoe UI, Roboto, sans-serif;
       color: #e8ecf4; background: #0b1020; }
a { color: #8ab4ff; text-decoration: none; }
a:hover { text-decoration: underline; }
header { padding: 3rem 1.5rem 1.5rem; text-align: center;
         background: radial-gradient(1200px 400px at 50% -10%, #1b2647 0%, #0b1020 70%); }
header h1 { margin: 0; font-size: 2.6rem; letter-spacing: -.02em; }
.tagline { max-width: 46rem; margin: .6rem auto 0; color: #aeb7cc; }
.links { margin-top: .8rem; color: #6b7796; }
main { max-width: 72rem; margin: 0 auto; padding: 1rem 1.5rem 4rem; }
h2 { margin: 2.4rem 0 .3rem; font-size: 1.45rem; }
h2 + .sub { margin: 0 0 1rem; color: #8b95ad; }
h3.section { margin: 1.6rem 0 .6rem; font-size: 1rem; text-transform: uppercase;
             letter-spacing: .08em; color: #8b95ad; font-weight: 600; }
.grid { display: grid; gap: 1.1rem;
        grid-template-columns: repeat(auto-fill, minmax(15rem, 1fr)); }
.card, .globe { background: #131a30; border: 1px solid #232c49; border-radius: 12px;
                overflow: hidden; display: flex; flex-direction: column; }
.card .thumb { display: block; background: #0c1226; aspect-ratio: 16 / 10; overflow: hidden; }
.card .thumb img { width: 100%; height: 100%; object-fit: cover; display: block; }
.card .body { padding: .8rem .9rem 1rem; display: flex; flex-direction: column; gap: .35rem; }
.card h4 { margin: 0; font-size: 1.02rem; }
.card p { margin: 0; color: #b9c1d6; font-size: .92rem; }
.badge { font-size: .72rem; color: #9fb0d6; }
.cmd code, .card code { background: #0c1226; border: 1px solid #232c49; border-radius: 6px;
                        padding: .12rem .4rem; font-size: .82rem; color: #cfe0ff; }
.cta { margin-top: .2rem; font-weight: 600; }
.globe { padding: 1.4rem 1.1rem; text-align: center; align-items: center; gap: .4rem;
         background: radial-gradient(360px 160px at 50% 0%, #1c2a55 0%, #131a30 75%); }
.globe .orb { font-size: 2.4rem; }
.globe .t { font-weight: 600; }
.globe .cta { color: #8ab4ff; }
footer { border-top: 1px solid #232c49; margin-top: 2.5rem; padding-top: 1.2rem;
         color: #8b95ad; font-size: .9rem; }
"""


def _header() -> str:
    return f"""\
<header>
  <h1>planet-sim</h1>
  <p class="tagline">Stellar flux and planet parameters in — climate, circulation, and biomes out.
  A clickable tour of every demonstration, interactive globe, and the teaching notebook.</p>
  <p class="links"><a href="{_REPO_URL}">GitHub repo ↗</a> ·
     <a href="{_BLOB}/README.md">README ↗</a> ·
     <a href="{_BLOB}/docs/plans/planet-earth-system.md">build plan ↗</a></p>
</header>
<main>"""


def _interactive_hero() -> str:
    """The headline call-to-action: the no-install browser what-if (docs/interactive/index.html)."""
    return ('<section>\n  <h2>Turn a knob, build a climate</h2>\n'
            '  <p class="sub">No install, no notebook — drag the Sun and the greenhouse and watch the '
            'planet respond, with a plain-language explanation of <em>what changed and why</em>. '
            'Every number is the real model.</p>\n'
            '  <a class="globe" href="interactive/index.html" '
            'style="display:block;max-width:42rem;text-decoration:none">\n'
            '    <span class="orb">🛠</span>\n'
            '    <span class="t">Build a climate — the interactive what-if</span>\n'
            '    <span class="cta">Open it in your browser ↗</span>\n'
            '  </a>\n</section>')


def _hero_globes(demos: tuple[Demo, ...]) -> str:
    cards = []
    for d in demos:
        href = globe_href(d)
        if href is None:
            continue
        cards.append(f"""\
    <a class="globe" href="{_rel_href(href)}">
      <span class="orb">🌍</span>
      <span class="t">{_esc(d.title)}</span>
      <span class="cta">Open interactive globe ↗</span>
    </a>""")
    return ('<section>\n  <h2>Explore the globes in 3D</h2>\n'
            '  <p class="sub">Drag to rotate, scroll to zoom, hover for the numbers — no install, '
            'they open straight in your browser.</p>\n'
            f'  <div class="grid">\n' + "\n".join(cards) + "\n  </div>\n</section>")


def _card(d: Demo) -> str:
    if d.artifact.endswith(".html"):                     # the demo's primary view *is* a globe
        thumb = (f'<a class="thumb" href="{_rel_href(d.artifact)}" '
                 f'style="display:flex;align-items:center;justify-content:center;font-size:2.4rem">🌍</a>')
    else:
        src = _rel_href(d.artifact)
        thumb = f'<a class="thumb" href="{src}"><img src="{src}" alt="{_esc(d.title)}" loading="lazy"></a>'

    extras = ("needs " + " ".join(f"<code>.[{e}]</code>" for e in d.extras)) if d.extras else ""
    sim = " · runs a short simulation" if d.sim else ""
    links = []
    g = globe_href(d)
    if g and not d.artifact.endswith(".html"):           # a still primary that also has a live globe
        links.append(f'<a class="cta" href="{_rel_href(g)}">Open interactive globe ↗</a>')
    return f"""\
    <article class="card">
      {thumb}
      <div class="body">
        <h4>{_esc(d.title)}</h4>
        <p>{_esc(d.blurb)}</p>
        <p class="badge">{extras}{sim}</p>
        <p class="cmd"><code>python -m planet {d.key}</code></p>
        {''.join(links)}
      </div>
    </article>"""


def _gallery(demos: tuple[Demo, ...]) -> str:
    out = ['<section>\n  <h2>Demonstrations</h2>\n'
           '  <p class="sub">Each prints its validation table and banks a figure; click a figure for '
           'the full image, or copy the command to reproduce it.</p>']
    last_section = None
    for d in demos:
        if d.section != last_section:
            if last_section is not None:
                out.append("  </div>")
            out.append(f'  <h3 class="section">{_esc(d.section)}</h3>')
            out.append('  <div class="grid">')
            last_section = d.section
        out.append(_card(d))
    out.append("  </div>\n</section>")
    return "\n".join(out)


def _notebook_section() -> str:
    return f"""\
<section>
  <h2>Experiment in the notebook</h2>
  <p class="sub">The notebook is the sandbox — this is where you run your own experiments.</p>
  <article class="card" style="max-width:42rem">
    <div class="body">
      <p><code>planet.ipynb</code> is the teaching narrative and the bench: build a world in the
      §7 <em>design-a-world</em> sandbox, drag the exoplanet / obliquity / size knobs, and
      predict-then-check each climate. It needs a running kernel, so launch it locally:</p>
      <p class="cmd"><code>python -m planet notebook</code></p>
      <p>Prefer to read it first?
         <a href="{_BLOB}/planet/planet.ipynb">View the notebook rendered on GitHub ↗</a>.
         The live-slider globe lives here too (<code>planet.planetmap.interactive_map()</code>).</p>
    </div>
  </article>
</section>"""


def _footer() -> str:
    adrs = " · ".join(
        f'<a href="{_BLOB}/docs/decisions/{n}">{n.split("-")[0]}</a>'
        for n in (
            "0001-language-and-performance.md",
            "0002-visualization-and-ux.md",
            "0003-test-execution-policy.md",
            "0004-interactive-maps-and-state-interchange.md",
            "0005-engines-are-living-contracts.md",
        )
    )
    return f"""\
<footer>
  <p>The science &amp; the decisions behind each leg:
     <a href="{_BLOB}/docs/plans/planet-earth-system.md">the build plan</a> ·
     ADRs {adrs}.</p>
  <p>This page is generated from the demo catalogue (<code>planet/catalog.py</code>) —
     regenerate it with <code>python -m planet site</code> after adding a demo; a test fails the
     build if it drifts.</p>
</footer>
</main>
</body>
</html>"""


def build_index_html(demos: tuple[Demo, ...] = DEMOS) -> str:
    """Render the whole landing page as a deterministic, self-contained HTML string."""
    head = (
        "<!doctype html>\n<html lang=\"en\">\n<head>\n"
        "<meta charset=\"utf-8\">\n"
        "<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
        "<title>planet-sim — a planetary climate simulator</title>\n"
        f"<style>\n{_CSS}</style>\n</head>\n<body>"
    )
    return "\n".join([
        head,
        _header(),
        _interactive_hero(),
        _hero_globes(demos),
        _gallery(demos),
        _notebook_section(),
        _footer(),
    ]) + "\n"


def write_site(path: Path = SITE_PATH) -> Path:
    """Write the landing page (LF newlines, UTF-8) and return its path."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(build_index_html(), encoding="utf-8", newline="\n")
    return path


def main() -> None:
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")
    saved = write_site()
    print(f"Landing page written → {saved.relative_to(_REPO_ROOT)}")
    print("  open it in a browser, or publish docs/ via GitHub Pages.")


if __name__ == "__main__":
    main()
