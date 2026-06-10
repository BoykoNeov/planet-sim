"""BigSim shared solver toolkit (ARCHITECTURE.md §5).

Each engine is a standalone, separately-tested library guarded by a passing
validation suite — a *living, versioned* contract (its ``CONTRACT.md``), extended
directly by consumers as needs grow, with the suite (and the full-repo gate on any
engine edit) as the guardrail rather than a freeze (ADR 0005). The first engine is
``engines.diffusion`` — the 1-D conservative parabolic (diffusion/heat) spine;
``engines.fluid`` (rotating shallow water) is the second.
"""
