"""Small archive-coordinate helpers; no signal processing or legacy imports."""
import numpy as np

def recorded_power(time, times, values, max_distance_s=30.):
    """Interpolate recorded meter samples, allowing their five-second cadence.
    A missing per-second file is not missing transmit power. Distant or invalid
    samples remain unavailable rather than extrapolating an unmeasured session.
    """
    times=np.asarray(times,dtype=float);values=np.asarray(values,dtype=float)
    if not len(times):return float('nan')
    if times.shape!=values.shape or np.any(np.diff(times)<0):raise ValueError('Invalid power time series')
    i=int(np.searchsorted(times,time));near=[j for j in [i-1,i]if 0<=j<len(times)]
    if min(abs(times[j]-time)for j in near)>max_distance_s:return float('nan')
    return float(np.interp(time,times,values))
