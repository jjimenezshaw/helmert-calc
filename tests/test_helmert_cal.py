import math
import pyproj
import pytest

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


def mgi_ferro_to_etrs89():
    input = {
        "source": {
            "pm": "ferro",
            "ellipsoid": {"name": "bessel"},
            "points": [
                [48.2200979, 33.3867950, 164.0],
                [48.2223853, 33.3510895, 162.0],
                [48.2027103, 33.3442230, 171.4],
                [48.1949297, 33.3661957, 160.5],
                [48.2104897, 33.3970947, 162.7],
            ],
        },
        "target": {
            "ellipsoid": {"name": "GRS80"},
            "points": [
                [48.2195734529, 15.7190118888, 208.9874],
                [48.2218593701, 15.6833110853, 207.0136],
                [48.2021864424, 15.6764470206, 216.4388],
                [48.1944075001, 15.6984173143, 205.5290],
                [48.2099667290, 15.7293109162, 207.6887],
            ],
        },
    }
    return input


def swap(input):
    input["source"], input["target"] = input["target"], input["source"]
    return input


def test_nothing():
    assert True


def test_mirrored_points():
    # the points are not properly paired
    input = {
        "source": {
            "ellipsoid": {"name": "GRS80"},
            "points": [[53, 8, 0], [52, 13, 0], [50, 12, 0]],
        },
        "target": {
            "ellipsoid": {"name": "GRS80"},
            "points": [[53, 8, 0], [50, 12, 0], [52, 13, 0]],
        },
    }
    sol = hc.helmert_calc(input)
    print(sol)
    assert math.fabs(sol["helmert"]["params"]["s"]) > 1000



@pytest.mark.parametrize("inverse", [False, True])
def test_dhdn_to_etrs89(inverse):
    # https://epsg.org/transformation_1776/DHDN-to-ETRS89-2.html

    input = basic_dhdn_to_etrs89()
    sign = 1
    if inverse:
        sign = -1
        input = swap(input)

    sol = hc.helmert_calc(input)
    params = sol["helmert"]["params"]
    assert sign * params["x"] == pytest.approx(598.10, abs=0.01)
    assert sign * params["y"] == pytest.approx(73.70, abs=0.01)
    assert sign * params["z"] == pytest.approx(418.20, abs=0.01)
    assert sign * params["rx"] == pytest.approx(0.202, abs=0.001)
    assert sign * params["ry"] == pytest.approx(0.045, abs=0.001)
    assert sign * params["rz"] == pytest.approx(-2.455, abs=0.001)
    assert sign * params["s"] == pytest.approx(6.70, abs=0.001)


@pytest.mark.parametrize("inverse", [False, True])
def test_dhdn_to_etrs89_pyproj(inverse):
    input = basic_dhdn_to_etrs89()
    if inverse:
        input = swap(input)

    sol = hc.helmert_calc(input)

    tr = pyproj.transformer.Transformer.from_pipeline(sol["helmert"]["pipeline_3D"])
    res = []
    for llh in input["source"]["points"]:
        res.append(tr.transform(*llh))
    for computed, data in zip(res, input["target"]["points"]):
        assert computed[0] == pytest.approx(data[0], abs=1e-7)
        assert computed[1] == pytest.approx(data[1], abs=1e-7)
        assert computed[2] == pytest.approx(data[2], abs=1e-3)


@pytest.mark.parametrize("inverse", [False, True])
def test_dhdn_to_etrs89_ellipsoid_by_params(inverse):
    input = basic_dhdn_to_etrs89()
    assert input["source"]["ellipsoid"]["name"] == "bessel"

    input["source"]["ellipsoid"]["name"] = "change name"
    input["source"]["ellipsoid"]["a"] = 6377397.155
    input["source"]["ellipsoid"]["rf"] = 299.1528128
    assert input["source"]["ellipsoid"]["name"] != "bessel"

    if inverse:
        input = swap(input)
        assert input["target"]["ellipsoid"]["name"] == "change name"

    sol = hc.helmert_calc(input)

    tr = pyproj.transformer.Transformer.from_pipeline(sol["helmert"]["pipeline_3D"])
    res = []
    for llh in input["source"]["points"]:
        res.append(tr.transform(*llh))
    for computed, data in zip(res, input["target"]["points"]):
        assert computed[0] == pytest.approx(data[0], abs=1e-7)
        assert computed[1] == pytest.approx(data[1], abs=1e-7)
        assert computed[2] == pytest.approx(data[2], abs=1e-3)


@pytest.mark.parametrize(
    "prime_meridian",
    [
        "ferro",
        -17.6666666666,
        pytest.param(-17, marks=pytest.mark.xfail),
        pytest.param("foo", marks=pytest.mark.xfail),
    ],
)
def test_mgi_ferro_to_etrs89(prime_meridian):
    # https://epsg.org/transformation_1619/MGI-to-ETRS89-1.html
    input = mgi_ferro_to_etrs89()
    input["source"]["pm"] = prime_meridian
    sol = hc.helmert_calc(input)

    tr = pyproj.transformer.Transformer.from_pipeline(sol["helmert"]["pipeline_3D"])
    res = []
    for llh in input["source"]["points"]:
        res.append(tr.transform(*llh))
    for computed, data in zip(res, input["target"]["points"]):
        assert computed[0] == pytest.approx(data[0], abs=1e-7)
        assert computed[1] == pytest.approx(data[1], abs=1e-7)
        assert computed[2] == pytest.approx(data[2], abs=2e-3)

    params = sol["helmert"]["params"]
    assert params["x"] == pytest.approx(577.326, abs=0.04)
    assert params["y"] == pytest.approx(90.129, abs=0.11)
    assert params["z"] == pytest.approx(463.919, abs=0.05)
    assert params["rx"] == pytest.approx(5.137, abs=0.003)
    assert params["ry"] == pytest.approx(1.474, abs=0.002)
    assert params["rz"] == pytest.approx(5.297, abs=0.002)
    assert params["s"] == pytest.approx(2.4232, abs=0.001)


@pytest.mark.parametrize("prime_meridian", ["ferro", -17.6666666666])
def test_mgi_ferro_to_etrs89_inverted(prime_meridian):
    # https://epsg.org/transformation_1619/MGI-to-ETRS89-1.html
    input = mgi_ferro_to_etrs89()
    input["source"]["pm"] = prime_meridian

    input = swap(input)
    sol = hc.helmert_calc(input)

    tr = pyproj.transformer.Transformer.from_pipeline(sol["helmert"]["pipeline_3D"])
    res = []
    for llh in input["source"]["points"]:
        res.append(tr.transform(*llh))
    for computed, data in zip(res, input["target"]["points"]):
        assert computed[0] == pytest.approx(data[0], abs=1e-8)
        assert computed[1] == pytest.approx(data[1], abs=1e-8)
        assert computed[2] == pytest.approx(data[2], abs=2e-3)

    params = sol["helmert"]["params"]
    assert params["x"] == pytest.approx(-577.326, abs=0.04)
    assert params["y"] == pytest.approx(-90.129, abs=0.11)
    assert params["z"] == pytest.approx(-463.919, abs=0.05)
    assert params["rx"] == pytest.approx(-5.137, abs=0.003)
    assert params["ry"] == pytest.approx(-1.474, abs=0.002)
    assert params["rz"] == pytest.approx(-5.297, abs=0.002)
    assert params["s"] == pytest.approx(-2.4232, abs=0.001)
