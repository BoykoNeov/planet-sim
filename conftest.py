"""Root pytest config: cap xdist parallelism at half the logical cores.

`addopts` runs the whole suite under `-n auto` (pytest-xdist) so every gate —
the local inner loop and CI alike — parallelises by default. Left to itself
`auto` would grab a worker per (physical) core; this hook instead pins the
count to half the machine's *logical* CPUs, so a test run leaves the other half
free for the rest of the system and the slow kernel/subprocess tests don't pile
onto a saturated box.

The hook fires only while xdist is resolving `-n auto` / `-n logical`; a serial
`pytest -n0` run (the documented escape hatch — it keeps xdist loaded, so the
inherited `-n auto` still parses) sets numprocesses to 0 and never calls it.
`optionalhook=True` is belt-and-suspenders: it lets this conftest load without
error even in an environment where xdist is absent entirely (e.g. a minimal env
that also clears addopts), where the hookspec would otherwise be unregistered.
"""

import os

import pytest


@pytest.hookimpl(optionalhook=True)
def pytest_xdist_auto_num_workers(config):
    return max(1, (os.cpu_count() or 2) // 2)
