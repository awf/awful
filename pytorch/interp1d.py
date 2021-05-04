def interp1d(x, y, xnew):
    """
    PyTorch interp1d, as in scipy.interpolate(x,y,assume_sorted=True)(xnew)
    
    Limited to 1D-1D for now, should not be hard to enhance if needed
    
    github.com/awf
    """
    # For each point in x, we want to find index of the knot above it, and hence that below it
    # so knots_x[ind-1] <= x < knots_x[ind]
    inds = torch.bucketize(xnew, x, right=True)

    # Call those points xlo, xhi
    xlo = x[inds-1]
    xhi = x[inds]
    ylo = y[inds-1]
    yhi = y[inds]
    
    dx = xhi - xlo
    dy = yhi - ylo

    # Then t = (xnew - xlo)/dx
    # ynew = ylo + t * dy
    t = (xnew - xlo) / dx
    return ylo + t * dy

def demo_interp1d():
    # scipy example from https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.interp1d.html
    import scipy.interpolate
    import numpy as np
    x = np.arange(0, 10)
    y = np.exp(-x/3.0)
    f = scipy.interpolate.interp1d(x, y)
    xnew = np.arange(0, 9, 0.1)
    ynew = f(xnew)   # use interpolation function returned by `interp1d`

    # Now do it with torch
    tx = torch.arange(0, 10)
    ty = torch.exp(-tx/3.0)
    txnew = torch.arange(0, 9, 0.1)
    tynew = interp1d(tx,ty,txnew)
    plt.figure()
    plt.plot(x, y, 'ro', xnew, ynew, 'r-',
             tx, ty, 'kx', txnew, tynew, 'k:')
    plt.show()

if __name__ == "__main__":
  demo_interp1d()
