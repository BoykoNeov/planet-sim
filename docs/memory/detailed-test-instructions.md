---
name: detailed-test-instructions
description: "When asking the user to test/verify something manually, always give step-by-step instructions (where, how to start, what to look at, expected result)"
metadata: 
  node_type: memory
  type: feedback
  originSessionId: e70e5026-9f6c-4a79-be0f-d41a9982163b
---

When I need the user to manually test or verify something (the repo's known
hand-off boundary is browser/WebGL play-through — no browser here), I must give
**detailed, concrete instructions**, not a vague "please eyeball it."

**Why:** The user can't act on "verify the page works." They need to know
exactly where to start, how to start it, where to look, and what a pass vs.
fail looks like — otherwise the hand-off stalls or they test the wrong thing.

**How to apply:** Every manual-test ask must spell out: (1) the exact command
or file to open and how (e.g. `python -m planet interactive`, or open
`docs/interactive/index.html` off disk); (2) what to do step by step (which
knob to drag, to what value); (3) what to look at (which element on screen);
(4) the expected result for each step + what failure looks like; (5) any
diagnostic to copy back (e.g. the browser console line) if it fails. Numbered
steps, not prose. Applies to [[interactive-what-if]] play-throughs and the
[[planet-viz-animation-rungs]] globe hand-offs.
