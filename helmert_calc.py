#!/usr/bin/env python3

import argparse
import json
import math

import numpy as np
import pyproj


def umeyama(X, Y):
    # https://github.com/clementinboittiaux/umeyama-python
    """
    Estimates the Sim(3) transformation between `X` and `Y` point sets.

    Estimates c, R and t such as c * R @ X + t ~ Y.

    Parameters
    ----------
    X : numpy.array
        (m, n) shaped numpy array. m is the dimension of the points,
        n is the number of points in the point set.
    Y : numpy.array
        (m, n) shaped numpy array. Indexes should be consistent with `X`.
        That is, Y[:, i] must be the point corresponding to X[:, i].

    Returns
    -------
    c : float
        Scale factor.
    R : numpy.array
        (3, 3) shaped rotation matrix.
    t : numpy.array
        (3, 1) shaped translation vector.
    """
    mu_x = X.mean(axis=1).reshape(-1, 1)
    mu_y = Y.mean(axis=1).reshape(-1, 1)
    var_x = np.square(X - mu_x).sum(axis=0).mean()
    cov_xy = ((Y - mu_y) @ (X - mu_x).T) / X.shape[1]
    U, D, VH = np.linalg.svd(cov_xy)
    S = np.eye(X.shape[0])
    if np.linalg.det(U) * np.linalg.det(VH) < 0:
        S[-1, -1] = -1
    c = np.trace(np.diag(D) @ S) / var_x
    R = U @ S @ VH
    t = mu_y - c * R @ mu_x
    return c, R, t


def solve_helmert(geoc_s, geoc_t):
    c, R, t = umeyama(np.transpose(geoc_s), np.transpose(geoc_t))
    # print("\nc", c, "\nR", R, "\nt", t)
    secs = 3600
    x = t[0][0]
    y = t[1][0]
    z = t[2][0]
    rx = math.degrees(R[2, 1]) * secs
    ry = math.degrees(R[0, 2]) * secs
    rz = math.degrees(R[1, 0]) * secs
    s = (c - 1) * 1e6
    return (x, y, z, rx, ry, rz, s)


def geographic_to_geocentric(ellipsoid_str, llhs):
    pipeline = (
        "+proj=pipeline "
        "+step +proj=axisswap +order=2,1 "
        "+step +proj=unitconvert +xy_in=deg +xy_out=rad "
        f"+step +proj=cart {ellipsoid_str}"
    )
    tr = pyproj.transformer.Transformer.from_pipeline(pipeline)
    res = []
    for llh in llhs:
        res.append(np.array(tr.transform(*llh)))
    return np.array(res)


def _read_ellipsoid(ellip):
    # Check if ellipsoid has +a and +rf defined, if not read from proj list
    # Include custom: boolean to determine in the pipeline if "+a and +rf"
    # or +ellps parameter is to be used
    el = None
    if "a" not in ellip or "rf" not in ellip:
        # print(f"Reading ellipsoid {ellip['name']} from pyproj")
        # https://proj.org/en/stable/usage/ellipsoids.html#ellipsoids
        ellps = pyproj.list.get_ellps_map()
        if ellip["name"] in ellps.keys():
            el_in_map = ellps[ellip["name"]]
            el = {
                "name": ellip["name"],
                "a": el_in_map["a"],
                "rf": el_in_map["rf"],
                "custom": False,
            }
            return el
    # no lookup necessary -> return user defined variables
    else:
        return {
            "name": ellip["name"],
            "a": ellip["a"],
            "rf": ellip["rf"],
            "custom": True,
        }


def _make_pipeline(pipe, ellps_s_str, ellps_t_str, pm_s_str, pm_t_str, is2D=False):
    p = (
        "+proj=pipeline "
        " +step +proj=axisswap +order=2,1 "
        " +step +proj=unitconvert +xy_in=deg +xy_out=rad "
        f" {pm_s_str} "
        f" {'+step +proj=push +v_3' if is2D else ''} "
        f" +step +proj=cart {ellps_s_str} "
        f" {pipe} "
        f" +step +inv +proj=cart {ellps_t_str} "
        f" {'+step +proj=pop +v_3' if is2D else ''} "
        f" {pm_t_str} "
        " +step +proj=unitconvert +xy_in=rad +xy_out=deg "
        " +step +proj=axisswap +order=2,1"
    )
    return p


def _build_ellipsoid_parameters(ellipsoid):
    if ellipsoid["custom"]:
        return f"+a={ellipsoid['a']} +rf={ellipsoid['rf']}"
    else:
        return f"+ellps={ellipsoid['name']}"


def _make_solution(x, y, z, rx, ry, rz, s, ellps_s_str, ellps_t_str, pm_s, pm_t):
    sol = {}
    helmert = sol["helmert"] = {}
    helmert["params"] = {"x": x, "y": y, "z": z, "rx": rx, "ry": ry, "rz": rz, "s": s}
    pipe = (
        f"+step +proj=helmert "
        f"+x={x:.3f} +y={y:.3f} +z={z:.3f} "
        f"+rx={rx:.4f} +ry={ry:.4f} +rz={rz:.4f} "
        f"+s={s:.3f} +convention=position_vector"
    )

    pm_s_str = f"+step +inv +proj=longlat +pm={pm_s}" if pm_s else ""
    pm_t_str = f"+step +proj=longlat +pm={pm_t}" if pm_t else ""

    helmert["string"] = pipe
    helmert["pipeline_2D"] = _make_pipeline(
        pipe,
        ellps_s_str,
        ellps_t_str,
        pm_s_str,
        pm_t_str,
        True,
    )
    helmert["pipeline_3D"] = _make_pipeline(
        pipe,
        ellps_s_str,
        ellps_t_str,
        pm_s_str,
        pm_t_str,
        False,
    )
    return sol


def _pm_correction(points, pm):
    if pm:
        pipeline = (
            "+proj=pipeline "
            "+step +proj=axisswap +order=2,1 "
            f"+step +inv +proj=longlat +pm={pm} "
            "+step +proj=axisswap +order=2,1 "
        )
    else:
        pipeline = "+proj=noop"

    tr = pyproj.transformer.Transformer.from_pipeline(pipeline)
    res = []
    for point in points:
        res.append(np.array(tr.transform(*point)))
    return np.array(res)


def _read_input(obj):
    ellps = _read_ellipsoid(obj["ellipsoid"])
    pm = obj.get("pm")
    points = _pm_correction(obj["points"], pm)
    ellipsoid_str = _build_ellipsoid_parameters(ellps)
    geoc = geographic_to_geocentric(ellipsoid_str, points)
    return ellipsoid_str, pm, geoc


def helmert_calc(input):
    ellps_s, pm_s, geoc_s = _read_input(input["source"])
    ellps_t, pm_t, geoc_t = _read_input(input["target"])
    hlmrt = solve_helmert(geoc_s, geoc_t)

    sol = _make_solution(*hlmrt, ellps_s, ellps_t, pm_s, pm_t)
    return sol


def main(input_filename, output_filename, inverse, quiet):
    with open(input_filename) as input_f:
        input = json.load(input_f)

    if inverse:
        input["source"], input["target"] = input["target"], input["source"]

    solution = helmert_calc(input)
    if not quiet:
        print(json.dumps(solution, indent=4))

    if output_filename:
        res = input
        res["solution"] = solution
        with open(output_filename, "w") as output_f:
            json.dump(res, output_f, indent=4)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate Helmert parameters.")
    parser.add_argument(
        "--input", "-i", type=str, help="file path for input values as JSON"
    )
    parser.add_argument(
        "--output", "-o", type=str, help="file path for output values as JSON"
    )
    parser.add_argument(
        "--inverse", default=False, action=argparse.BooleanOptionalAction, help="invert source and target"
    )
    parser.add_argument(
        "--quiet", "-q", default=False, action=argparse.BooleanOptionalAction, help="file path for output values as JSON"
    )
    args = parser.parse_args()

    main(args.input, args.output, args.inverse, args.quiet)
