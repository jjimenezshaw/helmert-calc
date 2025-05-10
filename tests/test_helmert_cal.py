import pytest
import pyproj

import helmert_calc as hc


def basic_dhdn_to_etrs89():
    input = {
        "source": {
            "ellipsoid": {"name": "bessel"},
            "points": [[53, 8, 0], [52, 13, 0], [50, 12, 0], [47, 7, 0], [49, 10, 0]],
        },
        "target": {
            "ellipsoid": {"name": "GRS80"},
            "points": [
                [52.9985065862, 7.9990898704, 41.7584],
                [51.9986488003, 12.9983316939, 42.5789],
                [49.9988693078, 11.99852119, 48.0281],
                [46.9991910348, 6.9992606373, 55.158],
                [48.9989715037, 9.9988250827, 51.0665],
            ],
        },
    }
    return input


def test_nothing():
    assert True


def test_basic_dhdn_to_etrs89():
    # https://epsg.org/transformation_1776/DHDN-to-ETRS89-2.html

    sol = hc.helmert_calc(basic_dhdn_to_etrs89())
    params = sol["helmert"]["params"]
    assert params["x"] == pytest.approx(598.10, 0.01)
    assert params["y"] == pytest.approx(73.70, 0.01)
    assert params["z"] == pytest.approx(418.20, 0.01)
    assert params["rx"] == pytest.approx(0.202, 0.001)
    assert params["ry"] == pytest.approx(0.045, 0.001)
    assert params["rz"] == pytest.approx(-2.455, 0.001)
    assert params["s"] == pytest.approx(6.70, 0.01)


def test_basic_dhdn_to_etrs89_pyproj():
    input = basic_dhdn_to_etrs89()
    sol = hc.helmert_calc(input)

    tr = pyproj.transformer.Transformer.from_pipeline(sol["helmert"]["pipeline_3D"])
    res = []
    for llh in input["source"]["points"]:
        res.append(tr.transform(*llh))
    for computed, data in zip(res, input["target"]["points"]):
        assert computed[0] == pytest.approx(data[0], 1e-8)
        assert computed[1] == pytest.approx(data[1], 1e-8)
        assert computed[2] == pytest.approx(data[2], 1e-3)
