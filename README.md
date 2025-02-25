# barnesn

This repo contains code for implementing Barnes objective analysis, which is used for interpolating scattered data points onto a grid.  The code is translated from the [Barnes interpolation (Barnes objective analysis) toolbox](https://www.mathworks.com/matlabcentral/fileexchange/58937-barnes-interpolation-barnes-objective-analysis), originally written in Matlab by Lena Bartell.

**Original description (from Matlab):** <br>
   BARNESN: Barnes smoothing interpolation of unstructured data
```Vq = BARNESN(X, V, Xv)``` returns the smoothing interpolation of
D-dimensional observations ```V(X)``` at query points ```Xq```. Query points ```Xq``` are
created by meshing the vectors in the cell array Xv that define the
grid in each dimension. Smoothing interpolation is performed using
the Koch form of Barnes objective analysis (Daley 1993). Roughly, (in 2D) the
interpolated value ```(vq)``` at gridpoint ```(xq, yq)``` is determined as a
weighted-sum of the values ```(v)``` at data points ```(x, y)```, based on the
gaussian weighting function ```exp(-r^2 / s / g^j)```, where r is the
euclidian distance from ```(xq, yq)``` to ```(x, y)```, s is the Gaussian Variance,
and ```g``` is the Convergence Parameter.

**Bibliography:**
* Barnes, Stanley L. "Mesoscale objective map analysis using weighted
time-series observations." (1973)
* Bartell, Lena (2025). Barnes interpolation (Barnes objective analysis) (https://www.mathworks.com/matlabcentral/fileexchange/58937-barnes-interpolation-barnes-objective-analysis), MATLAB Central File Exchange. Retrieved February 25, 2025.
* Daley, Roger. Atmospheric data anlysis. No. 2. Cambridge University
Press, 1993.
* Koch, Steven E., Mary DesJardins, and Paul J. Kocin. "An
interactive Barnes objective map analysis scheme for use with
satellite and conventional data." Journal of Climate and Applied
Meteorology 22.9 (1983): 1487-1503.
