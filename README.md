# helmert-calc

Calculator for Helmert parameters

## Usage

Run the main script `helmert_calc.py` with input and output variables.
Points and parameters are defined within a `.json` file.
Use the `a` and `rf` parameters to define custom ellipsoids.
See the example file [example-mueffling.json](./example-mueffling.json). The parameters for Müffling 1821 are taken from:

Straßer, Georg: Ellipsoidische Parameter der Erdfigur (1800-1950), München 1957, Verlag der Bayerischen Akademien der Wissenschaften, S. 31 ([PDF](https://dgk.badw.de/fileadmin/user_upload/Files/DGK/docs/a-19.pdf))

## Testing Data 

The historical testing data using actual historical coordinates and their modern counterparts in the files `example-MB25.json` and `example-NKH-25.json` 
were kindly provided by Andreas Heckmann. For details on the data collection see the 
respective article (in german): 

> Fletling, Rainer / Heckmann, Bernhard: Georeferenzierung historischer topografischer 
 Karten des 19. Jahrhunderts mit dem Bezugsmeridian von Ferro. In DVW-Mitteilungen
 Hessen-Thüringen, Heft 1/2018, S. 2–15.
