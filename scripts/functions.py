#!/usr/bin/python

# Function to compute a sectional grid

def sectionalGrid(yMin, yMax, N, gridType):

    import numpy as np

    x = np.zeros(N)
    y = np.zeros(N+1)

    if gridType == 'logarithmic':

        a = np.power(yMax/yMin, 1.0/N)

        for i in range(0,N):

            x[i] = yMin*np.power(a, float(i)+0.5)
            y[i] = yMin*np.power(a, float(i))

        y[N] = yMax

    elif gridType == 'linear':

        a = (yMax-yMin)/N

        for i in range(0,N):

            x[i] = yMin+(float(i)+0.5)*a
            y[i] = yMin+float(i)*a

        y[N] = yMax

    else:

        raise ValueError('Invalid gridType specified')

    return x,y

# Function to compute df/dlogd of the log-normal distribution in length, given an array of diameters

def logNormal(d, sigmag, CMD):

    import numpy as np

    return \
        1.0/(np.sqrt(2.0*np.pi)*np.log(sigmag)) \
      * np.exp \
        (
          - np.square(np.log(d/CMD))
          / (2.0*np.log(sigmag)**2)
        )

# Function to compute df/dlogd of the log-normal distribution in length, given a range in mass

def logNormalInterval(yMin, yMax, sigmag, CMM):

    import numpy as np
    from scipy.special import erf

    return 0.5* \
    (
        erf
        (
            np.log(yMax/CMM)
          / (3.0*np.sqrt(2.0)*np.log(sigmag))
        )
      - erf
        (
            np.log(yMin/CMM)
          / (3.0*np.sqrt(2.0)*np.log(sigmag))
        )
    ) \
  / (1.0/3.0*np.log(yMax/yMin))
