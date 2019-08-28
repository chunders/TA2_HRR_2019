import numpy as np

import re

pupil_x_re = re.compile(r'<x>(\d+)</x>')
pupil_y_re = re.compile(r'<y>(\d+)</y>')
pupil_r_re = re.compile(r'<radius>(.+)</radius>')

def read_has(fname):
    """Read the contents of a .has slopes file.

    Returns x_slopes, y_slopes, pupil_data
    """

    with open(fname) as f:
        return read_has_from_stream(f)

def read_has_from_stream(f):
    state = 0
    x_slopes = []
    y_slopes = []
    for line in f:
        if state == 0:
            if line.strip() == '<buffer>':
                state = 1
        elif state == 1:
            if line.strip() == '</buffer>':
                state = 2
            else:
                x_slopes.append(line)
        elif state == 2:
            if line.strip() == '<buffer>':
                state = 3
        elif state == 3:
            if line.strip() == '</buffer>':
                state = 4
            else:
                y_slopes.append(line)
        elif state == 4:
            if line.strip() == '<projection_pupil>':
                state = 5
        elif state == 5:
            line = line.strip()
            x_match = pupil_x_re.match(line)
            if x_match is not None:
                pupil_x = int(x_match.group(1))
                continue
            y_match = pupil_y_re.match(line)
            if y_match is not None:
                pupil_y = int(y_match.group(1))
                continue
            r_match = pupil_r_re.match(line)
            if r_match is not None:
                pupil_r = float(r_match.group(1))

    x_slopes = np.array([np.fromstring(x.strip(), sep='\t') for x in x_slopes])
    y_slopes = np.array([np.fromstring(y.strip(), sep='\t') for y in y_slopes])

    return x_slopes, y_slopes, {'x': pupil_x, 'y': pupil_y, 'r': pupil_r}

## Zernike support code

import scipy.special

def z_slopes_cos(x, y, rho, phi, n, m):
    # The (x,y)-slopes of the phase term rho**2 * cos(m*phi)
    rhon2 = rho**(n-2)
    cosmp = np.cos(m * phi)
    sinmp = np.sin(m * phi)
    x_slope = rhon2 * (  m*y*sinmp + n*x*cosmp)
    y_slope = rhon2 * (- m*x*sinmp + n*y*cosmp)
    return (x_slope, y_slope)

def z_slopes_sin(x, y, rho, phi, n, m):
    # The (x,y)-slopes of the phase term rho**2 * sin(m*phi)
    rhon2 = rho**(n-2)
    cosmp = np.cos(m * phi)
    sinmp = np.sin(m * phi)
    x_slope = rhon2 * (- m*y*cosmp + n*x*sinmp)
    y_slope = rhon2 * (  m*x*cosmp + n*y*sinmp)
    return (x_slope, y_slope)

def zernike(x, y, n, m):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    if m < 0:
        m = -m
        F = z_slopes_sin
    else:
        F = z_slopes_cos
    assert m <= n
    assert (n - m) % 2 == 0
    D = int((n - m) / 2)
    zx = None
    zy = None
    for k in range(D + 1):
        coeff = (-1)**k * scipy.special.binom(n - k, k) * scipy.special.binom(n - 2*k, D - k)
        termx, termy = F(x, y, rho, phi, n-2*k, m)
        termx *= coeff
        termy *= coeff
        if zx is None:
            zx = termx
            zy = termy
        else:
            zx += termx
            zy += termy
    return (zx, zy)

def generate_zernike_indices(start, stop):
    """Generate the Zernike modes indices from n=start to n=stop, inclusive."""

    for n in range(start, stop+1):
        for m in range(-n, n+1, 2):
            yield n, m

class ZernikeDecomposer:
    """Decomposes slopes data into Zernike modes.
    
    Re-uses intermediate calculations where possible.
    """
    
    def __init__(self, max_n, *, limit_modes=None):
        self._max_n = max_n
        self._limit_modes = limit_modes
        self._Minv_cache = {}
    
    def _populate_cache(self, key, x, y, Z_modes, valid, into_raw_slopes):
        if key in self._Minv_cache:
            return
        
        M = np.zeros((np.sum(valid * 2), len(Z_modes)))
        for i, mode in enumerate(Z_modes):
            z_x, z_y = zernike(x, y, *mode)
            M[:, i] = into_raw_slopes(z_x, z_y)

        Minv = np.linalg.pinv(M)
        self._Minv_cache[key] = Minv
    
    def decompose(self, x_slopes, y_slopes, pupil_data):
        x = x_slopes.shape[1]
        y = x_slopes.shape[0]

        y, x = np.mgrid[0:y, 0:x]

        x -= pupil_data['x']
        y -= pupil_data['y']
        x = x / pupil_data['r'] + 1e-9 # Needed to avoid divide-by-zero errors
        y = y / pupil_data['r'] + 1e-9

        r = np.sqrt(x**2 + y**2)
        valid = r < 1 & np.logical_not(np.isnan(x_slopes))
        
        # Use the pupil data and set of valid points to compute a key
        key = (tuple(pupil_data.items()), tuple(valid.flat))
        
        def into_raw_slopes(x_s, y_s):
            length = np.sum(valid)
            wf = np.zeros(length * 2)
            wf[:length] = x_s[valid]
            wf[length:] = y_s[valid]
            return wf
        
        Z_modes = list(generate_zernike_indices(1, self._max_n))
        if self._limit_modes is not None:
            Z_modes = Z_modes[:self._limit_modes]

        self._populate_cache(key, x, y, Z_modes, valid, into_raw_slopes)
        Minv = self._Minv_cache[key]
        
        magic_haso_multiplier = 1.6 # This seems to be needed to match the HASO's results with "Standard" normalization
        
        return Z_modes, Minv @ into_raw_slopes(x_slopes, y_slopes) * magic_haso_multiplier

def do_zernike_decomposition(x_slopes, y_slopes, pupil_data, max_n, *, limit_modes=None):
    """Perform a Zernike decomposition using a matrix pseudoinverse.

    Returns zernike_modes, decomposition
    """
    
    return ZernikeDecomposer(max_n, limit_modes=limit_modes).decompose(x_slopes, y_slopes, pupil_data)
    
    # The following is the original (non-caching) code

    x = x_slopes.shape[1]
    y = x_slopes.shape[0]

    y, x = np.mgrid[0:y, 0:x]

    #x -= pupil_data['x']
    #y -= pupil_data['y']
    x = x - pupil_data['x']
    y = y - pupil_data['y']
    x = x / pupil_data['r'] + 1e-9 # Needed to avoid divide-by-zero errors
    y = y / pupil_data['r'] + 1e-9

    r = np.sqrt(x**2 + y**2)
    valid = r < 1 & np.logical_not(np.isnan(x_slopes))

    def into_raw_slopes(x_s, y_s):
        length = np.sum(valid)
        wf = np.zeros(length * 2)
        wf[:length] = x_s[valid]
        wf[length:] = y_s[valid]
        return wf

    Z_modes = list(generate_zernike_indices(1, max_n))
    if limit_modes is not None:
        Z_modes = Z_modes[:limit_modes]
    M = np.zeros((np.sum(valid * 2), len(Z_modes)))

    for i, mode in enumerate(Z_modes):
        z_x, z_y = zernike(x, y, *mode)
        M[:, i] = into_raw_slopes(z_x, z_y)

    Minv = np.linalg.pinv(M)
    
    magic_haso_multiplier = 1.6 # This seems to be needed to match the HASO's results with "Standard" normalization

    return Z_modes, Minv @ into_raw_slopes(x_slopes, y_slopes) * magic_haso_multiplier

__all__ = ['read_has', 'read_has_from_stream', 'do_zernike_decomposition', 'ZernikeDecomposer']
