"""CLI smoke tests."""
import subprocess
import sys

from astrocats import __version__


def test_module_version():
    result = subprocess.run(
        [sys.executable, "-m", "astrocats", "--version"],
        capture_output=True,
        text=True,
        check=False,
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0
    assert __version__ in output


def test_module_help():
    result = subprocess.run(
        [sys.executable, "-m", "astrocats"],
        capture_output=True,
        text=True,
        check=False,
    )
    output = result.stdout + result.stderr
    assert "Generate catalogs" in output or "usage:" in output.lower()
