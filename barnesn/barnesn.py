'''
Translated from Matlab toolbox for Barnes interpolation (Barnes objective analysis)
Originally written in Matlab by Lena Bartell
Translated to Python by Haley Synan 
link to original toolbox: https://www.mathworks.com/matlabcentral/fileexchange/58937-barnes-interpolation-barnes-objective-analysis
Lena Bartell (2025). Barnes interpolation (Barnes objective analysis) (https://www.mathworks.com/matlabcentral/fileexchange/58937-barnes-interpolation-barnes-objective-analysis), MATLAB Central File Exchange. Retrieved February 25, 2025.


#ORIGINAL DESCRIPTION:
   BARNESN Barnes smoothing interpolation of unstructured data
Vq = BARNESN(X, V, Xv) returns the smoothing interpolation of
D-dimensional observations V(X) at query points Xq. Query points Xq are
created by meshing the vectors in the cell array Xv that define the
grid in each dimension. Smoothing interpolation is performed using
the Koch form of Barnes objective analysis [2]. Roughly, (in 2D) the
interpolated value (vq) at gridpoint (xq, yq) is determined as a
weighted-sum of the values (v) at data points (x, y), based on the
gaussian weighting function exp(-r^2 / s / g^j), where r is the
euclidian distance from (xq, yq) to (x, y), s is the Gaussian Variance,
and g is the Convergence Parameter.
-
Bibliography:
[1] Barnes, Stanley L. "Mesoscale objective map analysis using weighted
time-series observations." (1973)
[2] Koch, Steven E., Mary DesJardins, and Paul J. Kocin. "An
interactive Barnes objective map analysis scheme for use with
satellite and conventional data." Journal of Climate and Applied
Meteorology 22.9 (1983): 1487-1503.
[3] Daley, Roger. Atmospheric data anlysis. No. 2. Cambridge University
Press, 1993.
-
'''

from sklearn.metrics import mean_squared_error
from math import sqrt
from scipy.interpolate import LinearNDInterpolator
import pandas as pd
import numpy as np
import math
import matplotlib.pyplot as plt
import cartopy
import warnings
from scipy.spatial.distance import cdist

def calc_Verror(V, Vq, Xcell,ii, data): 
    interp = LinearNDInterpolator(list(zip(gridX, gridY)), Vq.flatten())
    interp_v = interp(data.lat,data.lon)
    Verr = V - interp_v 
    rmse = np.nanmean(Verr**2)
    outside =np.where(Verr.isna()==True)
    try:
        Verr[outside] = 0
    except:
        Verr[outside[0]] = 0
    print('Barnes iteration ' + str(ii) + ', average RMS error is ' + str(rmse))
    return rmse, Verr

def parse_inputs(X, V, Xv,n_interations=3, convergenceparam=0.2, gaussianvariance=float('nan')): 
    params={}
    params['iterations']=n_interations
    params['gaussianvariance']=gaussianvariance
    params['convergenceparameter']=convergenceparam #convergence parameter is to mitigate oversmoothing of data 
    
    if len(X) != len(V):
        raise Exception('The sizes of ''V'' and ''X'' do not match. V must have one value for each row in X')
    
    #remove data points with nan/inf
    remove = np.where(np.isfinite(X)==False)
    X[np.isfinite(X)]
    V[np.isfinite(V)]
    
    #setup parameters, store variable sizes
    params['D'] = X.ndim
    params['nData'] = len(X)
    params['grid_size'] = (len(Xv[0]),len(Xv[1]))
    params['nGrid'] = np.prod(params['grid_size'])
    
    
    A = np.prod(np.max(X, axis=0) - np.min(X,axis=0))
    M = params['nData']
    params['data_spacing'] = sqrt(A)*(1+sqrt(M))/(M-1) #put into dict params
    
    params['optimal_var']=(2*params['data_spacing']/math.pi)**2 * params['convergenceparameter']**(-params['iterations'])
    
    if math.isnan(params['gaussianvariance'])== True:
        params['gaussianvariance'] = params['optimal_var']
    else: 
        raise Exception('Gaussian variance is small. The optimal value is '+str(optimal_var))
    
    params['gaussianstd'] = sqrt(params['gaussianvariance'])


    # Warn if the grid spacing is not appropriate for the data spacing
    #Limits from Koch 1983: 1/3 <= (dn/[grid spacing]) <= 1/2, where dn is the average nearest-neighbor spacing of the data points
    min_grid_spacing = min((np.diff(Xv[0]).min(),np.diff(Xv[1]).min()))
    max_grid_spacing = max((np.diff(Xv[0]).max(),np.diff(Xv[1]).max()))
    if (min_grid_spacing/params['data_spacing']) < 0.333:
        warnings.warn('Grid spacing should be larger than ' + str(params['data_spacing']/3)+ '. Smallest grid spacing: '+ str(min_grid_spacing)+'. Data spacing: '
              + str(params['data_spacing']))
    if (max_grid_spacing/params['data_spacing']) > 0.5: 
        warnings.warn('Note that grid spacing can be smaller than ' + str(params['data_spacing']/2) + '. Largest grid spacing: ' + str(max_grid_spacing) + 
                      '. Data spacing: ' + str(params['data_spacing']))
        
    # Warn if some grid points are far from the data points
    r = np.round(cdist(Xq,X),decimals=4) 
    np.sum(r<=2*params['gaussianstd'])
    
    if any(np.sum(r<=2*params['gaussianstd'],axis=1)<3):
        warnings.warn('Some grid points are far from any data points. Consider modifying the grid.')
    return params, X, V, Xq 

def barnesn(X, V, Xv, data, n_interations=3, convergenceparam=0.2, gaussianvariance=float('nan')):
    params, X, V, Vq = parse_inputs(X,V,Xv,n_interations, convergenceparam, gaussianvariance)
    
    #set up for analysis 
    from scipy.spatial.distance import cdist
    r = np.round(cdist(Xq,X),decimals=4) #matches matlab
    outer_data = [1]*len(data)
    W=[]
    for ii in range(0,params['iterations']):
        w = np.exp(-r**2/params['gaussianvariance']/params['convergenceparameter']**ii)
        sum_w = np.repeat(np.sum(w,axis=1), len(X)).reshape(len(gridX),len(V))
        W.append(w/sum_w) #matches
    
    #first pass
    ii = 0
    outer_grid = [1]*len(gridX)
    f = np.tile(V.values, len(outer_grid)).reshape(len(outer_grid),len(V))
    Vq= np.sum(W[ii]*f,axis=1)
    
    #subsequent passes
    Xcell = (data.lat.values,data.lon.values)
    Verr = calc_Verror(V,Vq,Xcell,ii,data)[1]
    
    for ii in range(1,params['iterations']):
        f = np.tile(Verr.values, len(outer_grid)).reshape(len(outer_grid),len(V))
        Vq = Vq + np.sum(W[ii]*f,axis=1)
        Verr = calc_Verror(V,Vq,Xcell,ii,data)[1]
    Vq = Vq.reshape(grid_size[1],grid_size[0])
    #Vq = Vq.reshape(grid_size)
    return Vq, params