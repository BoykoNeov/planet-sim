"""Guards for the ``python -m planet`` front door (:mod:`planet.__main__`).

The launcher carries a hand-written catalogue (:data:`planet.__main__.DEMOS`) that must not drift
from the actual demo modules, and two paths a script/CI runs head-on: the no-argument *non-TTY*
fallback and ``list`` must print and return — never block on ``input()``. These pin both, plus the
one bit of real logic in the notebook launcher: turning a Jupyter server log line into the URL that
opens the notebook (the fix for the README's copy-a-URL wart).
"""
from __future__ import annotations

import importlib

import pytest

from planet import __main__ as launcher


def test_catalog_keys_are_unique():
    keys = [d.key for d in launcher.DEMOS]
    assert len(keys) == len(set(keys)), f"duplicate catalogue keys: {keys}"


@pytest.mark.parametrize("demo", launcher.DEMOS, ids=lambda d: d.key)
def test_every_catalog_entry_resolves_to_a_module_with_main(demo):
    """Each entry must point at a real importable module exposing a callable ``main()``.

    The demo modules import their render stacks (matplotlib / plotly) lazily inside ``main()``, so
    importing the module itself is headless-safe — a missing ``[viz]``/``[webviz]`` can't mask a
    typo'd module name here.
    """
    mod = importlib.import_module(demo.module)
    assert callable(getattr(mod, "main", None)), f"{demo.module} has no callable main()"


@pytest.mark.parametrize("demo", launcher.DEMOS, ids=lambda d: d.key)
def test_catalog_artifacts_are_well_formed(demo):
    artifact = launcher._REPO_ROOT / demo.artifact
    assert artifact.parent.name == "figures"
    assert artifact.suffix in {".png", ".gif", ".html"}


def test_list_path_returns_without_input(capsys):
    """`python -m planet list` prints the catalogue and exits 0 (no prompt)."""
    assert launcher.main(["list"]) == 0
    out = capsys.readouterr().out
    for demo in launcher.DEMOS:
        assert demo.key in out


def test_no_args_non_tty_prints_and_returns(monkeypatch, capsys):
    """No args on a non-TTY (piped/CI) must fall back to the catalogue, not hang on input()."""
    monkeypatch.setattr(launcher.sys.stdin, "isatty", lambda: False)
    # If this path ever called input(), reading captured stdin would raise — so reaching 0 proves it didn't.
    assert launcher.main([]) == 0
    assert "planet-sim" in capsys.readouterr().out


def test_unknown_key_is_exit_2(capsys):
    assert launcher.main(["definitely-not-a-demo"]) == 2
    assert "Unknown demo" in capsys.readouterr().out


def test_help_flag(capsys):
    assert launcher.main(["--help"]) == 0
    assert "python -m planet" in capsys.readouterr().out


def test_notebook_url_from_a_real_server_log_line():
    """The exact line the user saw must yield the /lab/tree/ URL that opens the notebook file."""
    line = ("[I 2026-06-12 12:35:38.584 ServerApp]     "
            "http://localhost:8889/lab?token=b61722fe27ae41c7999aa617bab23c592e37897371ddd8b5")
    assert launcher._notebook_url_from_log(line) == (
        "http://localhost:8889/lab/tree/planet/planet.ipynb"
        "?token=b61722fe27ae41c7999aa617bab23c592e37897371ddd8b5"
    )


def test_notebook_url_ignores_non_server_lines():
    assert launcher._notebook_url_from_log("[I jupyterlab | extension was successfully linked.") is None
