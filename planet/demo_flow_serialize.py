"""R1 banked artifact: the vector-flow **interchange** seam — one schema + one renderer, any producer (§9.3 / §11).

The §9.5 flow renderers consume a :class:`~planet.flow_globe.FlowField`. R1 (:mod:`planet.flow_serialize`)
makes that field **serializable** through the planet-spec schema and proves the seam is
**producer-agnostic** — the serialization and the renderer do not care *what* produced the field. This
demo shows it with **two producers run through the identical path**:

1. a **synthetic global** ``(u, v)`` (``is_global=True``) — a fabricated field standing in for a *future*
   ocean engine's output (the spin-out preview, plan §11): the banked, committed globe below; and
2. the model's **real emergent eddy band** (``is_global=False``) — embedded NH-sector-only, zeros
   elsewhere (the honesty edge), banked to ``outputs/`` (gitignored).

Both are serialized into **one** schema (a ``VECTOR_OVERLAY`` layer carrying the coverage-extent +
provenance), **round-tripped identically** (``load(save(spec)) == spec`` — the real proof), and painted
by the **one** generic :func:`planet.planetmap.render` (the kind-dispatching cone overlay) — no
per-producer special-casing. *Producer-agnosticism is the exact property the ClimaOcean spin-out relies
on:* an ocean engine's ``(u, v)`` + SST will flow through the same seam as the eddy band.

Run headless (saves the globe + the interchange artifacts, prints the round-trip proof):

    python -m planet.demo_flow_serialize
"""
from __future__ import annotations

from pathlib import Path

from . import flow_serialize as fs
from . import planet_spec as ps

_REPO_ROOT = Path(__file__).resolve().parents[1]
DOCS_FIGURE = _REPO_ROOT / "docs" / "figures" / "planet-flow-serialize.html"
OUTPUTS = _REPO_ROOT / "outputs"

PROV_SYNTHETIC = "synthetic analytic global flow"
PROV_EDDY = "planet.eddy_flux (saturated eddy frame)"

NX = NY = 64              # the eddy-channel resolution for the real-producer leg (modest → fast)
N_FRAMES = 4              # only the saturated snapshot is used — a few frames suffice


def _round_trip(field, provenance: str, stem: str) -> ps.PlanetSpec:
    """Serialize ``field`` into the planet-spec schema, write it, reload, and **assert** ``load(save) == spec``.

    The interchange artifact (``<stem>.json`` + ``<stem>.npz``) is written to ``outputs/`` so it can be
    inspected; the equality assertion *is* the producer-agnosticism proof, run live in the demo."""
    spec = fs.vector_spec_from_flow_field(field, provenance=provenance)
    OUTPUTS.mkdir(parents=True, exist_ok=True)
    ps.save(spec, OUTPUTS / stem)
    reloaded = ps.load(OUTPUTS / stem)
    assert reloaded == spec, f"round-trip identity FAILED for {provenance!r}"
    return spec


def _render_globe(spec: ps.PlanetSpec, field, path: Path) -> Path:
    """Paint a serialized vector spec as a speed surface + flow cones, with the field's honest caption."""
    from .planetmap import save_html
    caption = ("<b>How to read it.</b> Colour is <b>flow speed</b> and the arrows trace the "
               "<b>(u, v) velocity</b>. " + field.honesty)
    path.parent.mkdir(parents=True, exist_ok=True)
    return save_html(spec.view(), path, active=fs.SPEED_LAYER, caption=caption)


def main() -> None:
    import sys
    if hasattr(sys.stdout, "reconfigure"):
        sys.stdout.reconfigure(encoding="utf-8")     # °C, ², ⁻, → on legacy codepages

    print("\nThe vector-flow interchange seam — one schema, one renderer, any producer\n")

    # --- producer 1: a synthetic GLOBAL field (the committed globe; needs no sim) ------------------ #
    synthetic = fs.synthetic_flow_field()
    syn_spec = _round_trip(synthetic, PROV_SYNTHETIC, "flow-synthetic")   # asserts the round-trip (loud)
    print(f"  synthetic global : is_global=True   round-trip load(save)==spec  ✓   "
          f"({len(syn_spec.layers)} layers)")

    # --- producer 2: the model's REAL emergent eddy band (runs the short sim) ---------------------- #
    # Unguarded on purpose: _round_trip ASSERTS load(save)==spec, and a real round-trip regression must
    # surface, not be swallowed as "skipped" — the eddy sim runs on core (no optional dep to guard).
    from . import demo_eddy_life
    from .flow_globe import flow_field_from_eddy
    eddy_field = flow_field_from_eddy(demo_eddy_life.compute(nx=NX, ny=NY, n_frames=N_FRAMES).eddy)
    eddy_spec = _round_trip(eddy_field, PROV_EDDY, "flow-eddy")
    cov = eddy_spec.view().layer(fs.VECTOR_LAYER).style["coverage"]
    print(f"  real eddy band   : is_global={cov['is_global']}  round-trip load(save)==spec  ✓   "
          f"(NH sector {cov['lat_min']:.0f}–{cov['lat_max']:.0f}° lat, zeros elsewhere)")

    print("\n  → both producers serialize into the SAME schema and round-trip identically;")
    print("    coverage-extent (is_global) is only exercised because one producer is non-global.\n")

    # --- render both through the ONE generic renderer (cones); commit the synthetic globe ---------- #
    # Only the render needs Plotly ([webviz]); the round-trip proof above already ran on bare core.
    try:
        saved = _render_globe(syn_spec, synthetic, DOCS_FIGURE)
        print(f"Synthetic-global globe saved → {saved.relative_to(_REPO_ROOT)}")
        eddy_html = _render_globe(eddy_spec, eddy_field, OUTPUTS / "planet-flow-eddy.html")
        print(f"Real eddy-band globe saved   → {eddy_html.relative_to(_REPO_ROOT)} (gitignored)")
    except ImportError:
        print("(plotly not installed — install the webviz extra to render the globes: "
              "pip install -e .[webviz])")


if __name__ == "__main__":
    main()
