"""FIR design tools without needing the scipy.signal toolbox.  Also provides
the Type II FIRLS filter design.

    import fir
"""

__author__ = "Translated from the scipy C code source by Eric Gibbons"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 01/30/2017 $"
__copyright__ = "(c) 2017 Eric Gibbons"
__license__ = "Python"


import matplotlib.pyplot as plt
import numpy as np

import sigtools
from scipy import fftpack

def remez(numtaps, bands, desired, weight=None, Hz=1, type='bandpass',
               maxiter=25, grid_density=16):
    """
    NOTE:  this runs on a separate C-routine that needs to be compiled to 
    go along with this.  Also, this is pretty much lifted from the Scipy 
    library, and comes with the same license.  
    
    Calculate the minimax optimal filter using the Remez exchange algorithm.
    Calculate the filter-coefficients for the finite impulse response
    (FIR) filter whose transfer function minimizes the maximum error
    between the desired gain and the realized gain in the specified
    frequency bands using the Remez exchange algorithm.
    Parameters
    ----------
    numtaps : int
        The desired number of taps in the filter. The number of taps is
            the number of terms in the filter, or the filter order plus one.
    bands : array_like
        A monotonic sequence containing the band edges in Hz.
        All elements must be non-negative and less than half the sampling
        frequency as given by `Hz`.
    desired : array_like
        A sequence half the size of bands containing the desired gain
        in each of the specified bands.
            weight : array_like, optional
        A relative weighting to give to each band region. The length of
        `weight` has to be half the length of `bands`.
            Hz : scalar, optional
        The sampling frequency in Hz. Default is 1.
            type : {'bandpass', 'differentiator', 'hilbert'}, optional
        The type of filter:
          * 'bandpass' : flat response in bands. This is the default.
          * 'differentiator' : frequency proportional response in bands.
            * 'hilbert' : filter with odd symmetry, that is, type III
                        (for even order) or type IV (for odd order)
                        linear phase filters.
            maxiter : int, optional
        Maximum number of iterations of the algorithm. Default is 25.
            grid_density : int, optional
        Grid density. The dense grid used in `remez` is of size
        ``(numtaps + 1) * grid_density``. Default is 16.
    Returns
    -------
    out : ndarray
        A rank-1 array containing the coefficients of the optimal
        (in a minimax sense) filter.
    See Also
    --------
    firls
    firwin
    firwin2
    minimum_phase
    References
    ----------
            .. [1] J. H. McClellan and T. W. Parks, "A unified approach to the
            design of optimum FIR linear phase digital filters",
            IEEE Trans. Circuit Theory, vol. CT-20, pp. 697-701, 1973.
            .. [2] J. H. McClellan, T. W. Parks and L. R. Rabiner, "A Computer
           Program for Designing Optimum FIR Linear Phase Digital
            Filters", IEEE Trans. Audio Electroacoust., vol. AU-21,
            pp. 506-525, 1973.
    Examples
    --------
            We want to construct a filter with a passband at 0.2-0.4 Hz, and
    stop bands at 0-0.1 Hz and 0.45-0.5 Hz. Note that this means that the
    behavior in the frequency ranges between those bands is unspecified and
    may overshoot.
    >>> from scipy import signal
            >>> bpass = signal.remez(72, [0, 0.1, 0.2, 0.4, 0.45, 0.5], [0, 1, 0])
            >>> freq, response = signal.freqz(bpass)
    >>> ampl = np.abs(response)
    >>> import matplotlib.pyplot as plt
    >>> fig = plt.figure()
    >>> ax1 = fig.add_subplot(111)
            >>> ax1.semilogy(freq/(2*np.pi), ampl, 'b-')  # freq in Hz
    >>> plt.show()
    """
    # Convert type
    try:
        tnum = {'bandpass': 1, 'differentiator': 2, 'hilbert': 3}[type]
    except KeyError:
        raise ValueError("Type must be 'bandpass', 'differentiator', "
                                                          "or 'hilbert'")
    # Convert weight
    if weight is None:
        weight = [1] * len(desired)

    
        
    bands = np.asarray(bands).copy()
    return sigtools.RemezScipy(numtaps, bands.astype(float), desired.astype(float),
                                weight.astype(float), tnum, maxiter, grid_density)


def ls(numtaps, bands, desired, weight):
    """
    FIR filter design using least-squares error minimization.

    This calculates the coefficients for a Type 1 or II filter. The
    filter type is automatically determined based on whether or not 
    numtaps is even or odd.

    Calculate the filter coefficients for the linear-phase finite
    impulse response (FIR) filter which has the best approximation
    to the desired frequency response described by `bands` and
    `desired` in the least squares sense (i.e., the integral of the
    weighted mean-squared error within the specified bands is
    minimized).

    Design borrows heavily from the Matlab implementation (firls.m).

    Parameters
    ----------
    numtaps : int
        The number of taps in the FIR filter.  `numtaps` must be odd.
    bands : array_like
        A monotonic nondecreasing sequence containing the band edges in
        Hz. All elements must be non-negative and less than or equal to
        the Nyquist frequency.  Assume 1Hz sampling.
    desired : array_like
        A sequence the same size as `bands` containing the desired gain
        at the start and end point of each band.
    weight : array_like, optional
        A relative weighting to give to each band region when solving
        the least squares problem. `weight` has to be half the size of
        `bands`.

    Returns
    -------
    coeffs : ndarray
        Coefficients of the optimal (in a least squares sense) FIR filter.


    """

    N = np.copy(numtaps)
    F = np.copy(bands)
    M = np.copy(desired)
    W = np.copy(weight)

    Nodd = N % 2

    W = np.sqrt(W)
    dF = np.diff(F)
    
    L = (N - 1)/2

    m = np.arange(0,L+1) + 0.5*(1 - Nodd)
    k = m
    
    I11, I22 = np.meshgrid(k,m,)

    I1 = I11 + I22
    I2 = - I11 + I22

    if Nodd:
        k = k[1:]
        b0 = 0

    b = np.zeros_like(k)
    G = np.zeros_like(I1)

    for s in np.arange(0,len(F),2):
        m = (M[s+1] - M[s])/(F[s+1] - F[s]) # slope
        b1 = M[s] - m*F[s]
        
        if Nodd:
            b0 += (b1*(F[s+1]-F[s]) + m/2*(F[s+1]**2 - F[s]**2)*abs(W[s/2])**2)
        
        b += ((m/(4*np.pi**2)*(np.cos(2*np.pi*k*F[s+1]) - np.cos(2*np.pi*k*F[s]))/(k**2))
              * abs(W[s/2])**2)
        
        b += ((F[s+1]*(m*F[s+1]+b1)*np.sinc(2*k*F[s+1])
               - F[s]*(m*F[s]+b1)*np.sinc(2*k*F[s]))
               * abs(W[s/2])**2)

        
        G += ((0.5*F[s+1]*(np.sinc(2*I1*F[s+1]) + np.sinc(2*I2*F[s+1]))
                 - 0.5*F[s]*(np.sinc(2*I1*F[s]) + np.sinc(2*I2*F[s])))
                * abs(W[s/2])**2)
        
    if Nodd:
        b = np.hstack((b0, b))

    
        
    a = np.linalg.pinv(G).dot(b)

    if Nodd:
        h = np.hstack((a[L:0:-1]/2, a[0], a[1:L+1]/2))
    else:
        h = 0.5*np.hstack((a[::-1],a))

    coeffs = h
        
    return coeffs

def _mag2mp(x):
    """
    
    Arguments:
    - `x`: magnitude of the analytic signal fft
    """

    n = len(x)
    xl = np.log(x)
    xlf = np.fft.fft(xl)
    xlfp = np.zeros((n),dtype=complex)
    xlfp[0] = xlf[0]
    xlfp[1:(n/2)] = 2*xlf[1:(n/2)]
    xlfp[n/2] = xlf[n/2]
    xlaf = np.fft.ifft(xlfp)
    a = np.exp(xlaf)

    return a
    
def MinPhase(h):
    """
    
    Arguments:
    - `h`:
    """

    l = len(h)

    if (l % 2) is 0:
        raise ValueError('Filter must be odd')

    lp = 8*np.exp(np.ceil(np.log(l)/np.log(2))*np.log(2))
    
    leftPad = np.zeros((int(np.ceil((lp-l)/2))))
    rightPad = np.zeros((int(np.floor((lp-l)/2))))
    hp = np.hstack((leftPad,h,rightPad))

    hpf = np.fft.fftshift(np.fft.fft(np.fft.fftshift(hp)))
    hpfs = hpf - np.amin(hpf.real)*1.000001

    hpfmp = _mag2mp(np.sqrt(abs(hpfs)))
    hpmp = np.fft.ifft(np.fft.fftshift(np.conjugate(hpfmp)))
    hmp = hpmp[0:(l+1)/2]

    return hmp
    
    

    
