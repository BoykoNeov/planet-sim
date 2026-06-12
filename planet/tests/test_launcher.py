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


# --- the interactive menu (the headline `python -m planet` surface) --------- #
def test_menu_loop_dispatches_by_number_name_and_quits(monkeypatch):
    """Drive _menu_loop with scripted input: a number, a name, the globe & notebook verbs, then quit."""
    calls: list = []
    inputs = iter(["1", "biomes", "g", "n", "q"])
    monkeypatch.setattr(launcher, "_read", lambda prompt="": next(inputs))
    monkeypatch.setattr(launcher, "_run_demo", lambda demo, offer_open: calls.append(("demo", demo.key)))
    monkeypatch.setattr(launcher, "_open_banked_globe", lambda *a, **k: calls.append(("globes",)))
    monkeypatch.setattr(launcher, "_launch_notebook", lambda: calls.append(("notebook",)))
    launcher._menu_loop()
    assert calls == [("demo", "snowball"), ("demo", "biomes"), ("globes",), ("notebook",)]


def test_menu_loop_reports_bad_choice_then_continues(monkeypatch, capsys):
    inputs = iter(["nope", "q"])
    monkeypatch.setattr(launcher, "_read", lambda prompt="": next(inputs))
    launcher._menu_loop()
    assert "not a choice" in capsys.readouterr().out


def test_menu_loop_survives_ctrl_c_at_prompt(monkeypatch):
    """Ctrl-C at `planet>` must quit cleanly, not blow up the loop with a traceback."""
    def interrupt(prompt=""):
        raise KeyboardInterrupt
    monkeypatch.setattr(launcher, "_read", interrupt)
    launcher._menu_loop()  # returns rather than raising


@pytest.mark.parametrize("answer,default,expected", [
    ("", True, True), ("", False, False),
    ("y", False, True), ("yes", False, True),
    ("n", True, False), ("no", True, False),
])
def test_prompt_yes_no(monkeypatch, answer, default, expected):
    monkeypatch.setattr(launcher, "_read", lambda prompt="": answer)
    assert launcher._prompt_yes_no("?", default=default) is expected


def test_open_banked_globe_by_choice(monkeypatch):
    opened: list = []
    monkeypatch.setattr(launcher, "_open_in_browser", lambda path, label="": opened.append(path.name))
    launcher._open_banked_globe("2")                       # 2 = the coupler globe
    assert opened == ["planet-coupler-map.html"]


def test_open_banked_globe_rejects_out_of_range(monkeypatch, capsys):
    monkeypatch.setattr(launcher, "_open_in_browser", lambda *a, **k: pytest.fail("should not open"))
    launcher._open_banked_globe("99")
    assert "no such globe" in capsys.readouterr().out


def test_site_verb_builds_and_opens(monkeypatch):
    """`python -m planet site` wires write_site() → open. Stubbed so the test writes/opens nothing."""
    opened: list = []
    # Don't touch the shared docs/index.html (xdist-safe) and don't spawn a browser.
    monkeypatch.setattr("planet.site.write_site", lambda: launcher._REPO_ROOT / "docs" / "index.html")
    monkeypatch.setattr(launcher, "_open_in_browser", lambda path, label="": opened.append(path.name))
    assert launcher.main(["site"]) == 0
    assert opened == ["index.html"]
