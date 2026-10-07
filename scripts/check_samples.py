"""Recompute every recorded JSON/CSV value in an isolated output directory."""

import json
import tempfile
from pathlib import Path

import numpy as np

from heatstart.study import run


def compare(actual, expected, path="metrics"):
    if isinstance(expected, dict):
        assert actual.keys() == expected.keys(), path
        for key in expected:
            compare(actual[key], expected[key], f"{path}.{key}")
    elif isinstance(expected, list):
        assert len(actual) == len(expected), path
        for i, (a, e) in enumerate(zip(actual, expected, strict=True)):
            compare(a, e, f"{path}[{i}]")
    else:
        np.testing.assert_allclose(actual, expected, rtol=5e-5, atol=5e-12, err_msg=path)


def main():
    root = Path(__file__).resolve().parents[1]
    expected = root / "examples/results"
    (root / "tmp").mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(dir=root / "tmp") as folder:
        actual = run(folder)
        compare(actual, json.loads((expected / "metrics.json").read_text()))
        for path in expected.glob("*.csv"):
            np.testing.assert_allclose(
                np.loadtxt(Path(folder) / path.name, delimiter=",", skiprows=1),
                np.loadtxt(path, delimiter=",", skiprows=1),
                rtol=2e-10,
                atol=5e-12,
            )
    print("All recorded JSON metrics and complete CSV histories reproduced.")


if __name__ == "__main__":
    main()
