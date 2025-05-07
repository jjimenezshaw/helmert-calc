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


def geographic_to_geocentric(ellpsoid, llhs):
    pipeline = (
        "+proj=pipeline "
        "+step +proj=axisswap +order=2,1 "
        "+step +proj=unitconvert +xy_in=deg +xy_out=rad "
        f"+step +proj=cart +ellps={ellpsoid["name"]}"
    )
    tr = pyproj.transformer.Transformer.from_pipeline(pipeline)
    res = []
    for llh in llhs:
        res.append(np.array(tr.transform(*llh)))
    return np.array(res)


def read_ellipsoid(ellip):
    # Check if ellipsoid has +a and +rf defined, if not read from proj list
    el = None
    if "a" not in ellip or "rf" not in ellip:
        print(f"Reading ellipsoid {ellip['name']} from pyproj")
        # https://proj.org/en/stable/usage/ellipsoids.html#ellipsoids
        ellps = pyproj.list.get_ellps_map()
        # @Javier: could be removed?
        # e = pyproj.crs.Ellipsoid.from_name("WGS 84")
        # print(e.to_json(True))
        if ellip["name"] in ellps.keys():
            el_in_map = ellps[ellip["name"]]
            el = {
                "type": "Ellipsoid",
                "name": ellip["name"],
                "semi_major_axis": el_in_map["a"],
                "inverse_flattening": el_in_map["rf"],
            }
            return el
    # no lookup necessary -> return user defined variables
    else:
        return {
            "type": "Ellipsoid",
            "name": ellip["name"],
            "semi_major_axis": ellip["a"],
            "inverse_flattening": ellip["rf"]
        }


def solve_helmert(geoc_s, geoc_t):
    def f(p):
        r = []
        for i in range(0, 3):
            r.append(np.array([x[i] for x in p]))
        return np.array(r)

    src = f(geoc_s)
    tgt = f(geoc_t)
    c, R, t = umeyama(src, tgt)
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
    # pipe = f"+step +proj=helmert +x={x:.3f} +y={y:.3f} +z={z:.3f} +rx={rx:.4f} +ry={ry:.4f} +rz={rz:.4f} +s={s:.3f} +convention=position_vector"
    # return pipe


def make_pipeline(pipe, ellps_s_str, ellps_t_str, is2D=False):
    p = (
        "+proj=pipeline "
        " +step +proj=axisswap +order=2,1 "
        " +step +proj=unitconvert +xy_in=deg +xy_out=rad "
        f" {"+step +proj=push +v_3" if is2D else ""} "
        f" +step +proj=cart {ellps_s_str} "
        f" {pipe} "
        f" +step +inv +proj=cart {ellps_t_str} "
        f" {"+step +proj=pop +v_3" if is2D else ""} "
        " +step +proj=unitconvert +xy_in=rad +xy_out=deg "
        " +step +proj=axisswap +order=2,1"
    )
    return p


def make_solution(x, y, z, rx, ry, rz, s, ellps_s, ellps_t):
    sol = {}
    helmert = sol["helmert"] = {}
    helmert["params"] = {"x": x, "y": y, "z": z, "rx": rx, "ry": ry, "rz": rz, "s": s}
    pipe = f"+step +proj=helmert +x={x:.3f} +y={y:.3f} +z={z:.3f} +rx={rx:.4f} +ry={ry:.4f} +rz={rz:.4f} +s={s:.3f} +convention=position_vector"
    helmert["string"] = pipe
    helmert["pipeline_2D"] = make_pipeline(
        pipe, f"+ellps={ellps_s["name"]}", f"+ellps={ellps_t["name"]}", True
    )
    helmert["pipeline_3D"] = make_pipeline(
        pipe, f"+ellps={ellps_s["name"]}", f"+ellps={ellps_t["name"]}", False
    )
    return sol


def main(input_filename, output_filename):
    with open(input_filename) as input_f:
        input = json.load(input_f)

    ellps_s = read_ellipsoid(input["source"]["ellipsoid"])
    ellps_t = read_ellipsoid(input["target"]["ellipsoid"])
    geoc_s = geographic_to_geocentric(ellps_s, input["source"]["points"])
    geoc_t = geographic_to_geocentric(ellps_t, input["target"]["points"])
    # print(geoc_s, geoc_t)
    hlmrt = solve_helmert(geoc_s, geoc_t)
    sol = make_solution(*hlmrt, ellps_s, ellps_t)
    print(sol)

    res = input
    res["solution"] = sol
    with open(output_filename, "w") as output_f:
        json.dump(res, output_f, indent=4)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate Helmert parameters.")
    parser.add_argument(
        "-input", "-i", type=str, help="file path for input values as JSON"
    )
    parser.add_argument(
        "-output", "-o", type=str, help="file path for output values as JSON"
    )
    args = parser.parse_args()

    main(args.input, args.output)
    """+proj=pipeline
  +step +proj=axisswap +order=2,1
  +step +proj=unitconvert +xy_in=deg +xy_out=rad
  +step +proj=push +v_3
  +step +proj=cart +ellps=bessel
  +step +proj=helmert +x=598.1 +y=73.7 +z=418.2 +rx=0.202 +ry=0.045 +rz=-2.455
        +s=6.7 +convention=position_vector
  +step +inv +proj=cart +ellps=GRS80
  +step +proj=pop +v_3
  +step +proj=unitconvert +xy_in=rad +xy_out=deg
  +step +proj=axisswap +order=2,1"""
