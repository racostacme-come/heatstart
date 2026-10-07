import json

import numpy as np

from heatstart.study import run


def test_recorded_experiment(tmp_path):
    result = run(tmp_path)
    methods = result["pulse"]["methods"]
    assert methods["cn"]["first_min"] < -0.5
    assert 0.50 < methods["cn"]["first_energy_ratio"] < 0.51
    assert methods["rannacher"]["all_min"] >= -1e-14
    assert methods["rannacher"]["final_time_error"] < methods["cn"]["final_time_error"] / 300
    for row in methods.values():
        assert row["relative_heat_drift"] < 1e-11
        assert row["max_energy_increase"] < 0
    assert result["matrix_exponential_max_error"] < 1e-12
    for method, order in result["temporal_last_orders"].items():
        assert abs(order - (1 if method == "be" else 2)) < 0.01
    assert np.all(np.array(result["spatial_orders"]) > 1.99)
    assert json.loads((tmp_path / "metrics.json").read_text()) == result
    assert len(list(tmp_path.glob("*.png"))) == 3
    assert all(b"\r" not in p.read_bytes() for p in tmp_path.glob("*.csv"))
    assert np.loadtxt(tmp_path / "history_cn.csv", delimiter=",", skiprows=1).shape == (13, 129)
