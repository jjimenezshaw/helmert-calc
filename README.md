# helmert-calc

Calculator for Helmert parameters


## Purpose

Some old or historical data are using not well defined references (definitely not the same as nowadays).
There isn't a well known transformation to a modern CRS like [ETRS89](https://spatialreference.org/ref/epsg/4258/) or [WGS 84](https://www.youtube.com/watch?v=M2ck3cAGvhg).

Given two sets of geographic coordinates -lat, lon, alt- (of the same real points), this software produces a [Helmert transformation](https://en.wikipedia.org/wiki/Helmert_transformation) to convert from one reference system to the other.
That transformation can be used in [PROJ](https://proj.org) and therefore in [pyproj](https://pyproj4.github.io/pyproj)
and [QGIS](https://qgis.org)


## Usage

Run the main script `helmert_calc.py` with input and output variables.
Points and parameters are defined within a `.json` file.
Use the `a` and `rf` parameters to define custom ellipsoids.
In case any set of coordinates use a primer meridian other than Greenwich (like Ferro), define it with the parameter `pm`, as in [example-MB25.json](./example-MB25.json).
You can use an name from [PROJ prime meridians](https://proj.org/en/stable/usage/projections.html#prime-meridian) or directly its longitude.

See the example file [example-mueffling.json](./example-mueffling.json). The parameters for Müffling 1821 are taken from:

Straßer, Georg: Ellipsoidische Parameter der Erdfigur (1800-1950), München 1957, Verlag der Bayerischen Akademien der Wissenschaften, S. 31 ([PDF](https://dgk.badw.de/fileadmin/user_upload/Files/DGK/docs/a-19.pdf))


## Residuals

As a check of the quality of the transformation, the residuals are computed in a projected reference system.
Those residuals, in meters, give you an idea about how far (or near) are the source points from the target points once the transformation is applied.


## Testing Data

The historical testing data using actual historical coordinates and their modern counterparts in the files `example-MB25.json` and `example-NKH-25.json`
were kindly provided by Andreas Heckmann. For details on the data collection see the
respective article (in german):

> Fletling, Rainer / Heckmann, Bernhard: Georeferenzierung historischer topografischer
 Karten des 19. Jahrhunderts mit dem Bezugsmeridian von Ferro. In DVW-Mitteilungen
 Hessen-Thüringen, Heft 1/2018, S. 2–15.
