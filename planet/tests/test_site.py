"""Guards for the generated landing page (:mod:`planet.site`).

Two jobs: prove the page is built *from* the catalogue (so it can't omit a demo), and — the one
that actually enforces "keep it up to date" — prove the committed ``docs/index.html`` equals a fresh
generation, so adding a demo without re-running ``python -m planet site`` fails the build.
"""
from __future__ import annotations

import re

from planet import site
from planet.catalog import DEMOS, globe_href


def test_committed_index_is_up_to_date():
    """docs/index.html must equal a fresh build — regenerate with `python -m planet site` if this fails."""
    expected = site.build_index_html()
    actual = site.SITE_PATH.read_text(encoding="utf-8")     # universal newlines neutralize CRLF
    assert actual == expected, "docs/index.html is stale — run `python -m planet site` and commit it"


def test_site_covers_every_catalogue_entry():
    page = site.build_index_html()
    for d in DEMOS:
        assert site._esc(d.title) in page, f"{d.key}: title missing from page"
        assert site._esc(d.blurb) in page, f"{d.key}: blurb missing from page"
        assert f"python -m planet {d.key}" in page, f"{d.key}: reproduce command missing"
        assert site._rel_href(d.artifact) in page, f"{d.key}: artifact link missing"
        g = globe_href(d)
        if g:
            assert site._rel_href(g) in page, f"{d.key}: interactive globe link missing"


def test_all_three_interactive_globes_are_linked():
    page = site.build_index_html()
    for href in ("figures/planet-map.html",
                 "figures/planet-coupler-map.html",
                 "figures/planet-eddy-globe.html"):
        assert href in page


def test_media_hrefs_are_relative_and_posix():
    """Figure/globe links must be docs-relative with forward slashes (work offline + on GitHub Pages)."""
    page = site.build_index_html()
    media = [h for h in re.findall(r'(?:href|src)="([^"]+)"', page) if h.startswith("figures/")]
    assert media, "expected relative figures/ links"
    assert all("\\" not in h for h in media), "hrefs must use forward slashes"


def test_rel_href_strips_docs_prefix():
    assert site._rel_href("docs/figures/planet-map.html") == "figures/planet-map.html"


def test_interactive_what_if_is_linked():
    """The headline no-install what-if must be linked with a relative (offline + Pages) href."""
    page = site.build_index_html()
    assert 'href="interactive/index.html"' in page
    assert "what changed and why" in page.lower()


def test_notebook_is_linked_not_run():
    """The notebook gets a GitHub-render link + the local command — never a fake in-page 'run'."""
    page = site.build_index_html()
    assert "planet/planet.ipynb" in page
    assert "python -m planet notebook" in page
