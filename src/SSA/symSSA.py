# symSSA.py: offered as a symbolic toolbox for the structural sensitivity analysis.

import numpy as np
import sympy as sp
from string import ascii_uppercase

class CRN(object):

    """
    An object that defines a Chemical Reaction Network with a given stoichiometric matrix.

        stoi        (2D numpy.ndarray)
                    the stoichiometrix matrix with its entries assuming to be integers.
        varinames   (1D iterable)
                    the names/notations of the chemical species.
        kineticmask (string, or 2D array-like)
                    the condition in which the partial derivatives can be omitted.

    Example:
    >>> nu
    array([[ 1, -1,  1, -1],
        [ 0,  1, -1,  0],
        [ 0,  1, -1,  0]])
    >>> crn = CRN(nu)
    >>> crn.A
    Matrix([
    [         0,          0,          0,  0, -1],
    [r_{2x_{1}},          0,          0, -1,  0],
    [         0, r_{3x_{2}}, r_{3x_{3}}, -1,  0],
    [r_{4x_{1}},          0,          0,  0, -1],
    [         0,         -1,          1,  0,  0]])
    """

    def __init__(self, stoi, varinames = None, kineticmask = ">="):

        assert stoi.dtype == np.dtype('int64')
        self._nu = stoi
        self._M, self._N = stoi.shape

        if varinames is None:
            self._X = [sp.symbols('x_{%d}'%(i+1)) for i in range(stoi.shape[0])]
        else:
            self._X = [sp.symbols(sp.latex('{}'.format(var))) for var in varinames]

        if type(kineticmask) == np.ndarray:
            self._kineticmask = kineticmask.astype(bool)
        else:
            assert isinstance(kineticmask, str)
            self._kineticmask = np.array([
                [eval("{} {} 0".format(stoi[m, n], kineticmask)) for n in range(self._N)]
            for m in range(self._M)])
        
        self._A = _AugMat(stoi, self._X, self._kineticmask)

    @property
    def varinames(self):
        return self._X
    
    @varinames.setter
    def varinames(self, varinames):
        self._X = varinames
        self._A = _AugMat(self._nu, varinames, self._kineticmask)
        
    @property
    def kineticmask(self):
        return self._kineticmask

    @kineticmask.setter
    def kineticmask(self, kineticmask):
        self._kineticmask = kineticmask
        self._A = _AugMat(self._nu, self._X, self._kineticmask)

    @property
    def A(self):
        return self._A
    
    def Jacobian(self):
        return sp.Matrix(self._nu.dot(self._A[:self._N,:self._M]))

def __KerImg(mat):
    """
    Bases of the image and null spaces for a given matrix (composed of integers) by Gaussian elimination.
    """
    M, N = mat.shape
    def _gaussian(ErMat, m, n):
        "Gaussian elimination on the m-th row for except for the first n's columns."
        headrow  = ErMat[N+m,n:]
        minval   = np.min(np.abs(headrow)[headrow!=0])
        k        = (np.arange(N-n)[np.abs(headrow)==minval])[0]
        if k != 0:
            ErMat[:, [n,n+k]] = ErMat[:, [n+k,n]] # swap
        for i in np.arange(0, N):
            if i == n:
                pass
            else:
                if ErMat[N+m,i] != 0:
                    lcm = np.lcm(headrow[0], ErMat[N+m,i])
                    ErMat[:,i] = ((lcm/ErMat[N+m,i]) * ErMat[:,i]).astype(int)
                    ErMat[:,i] = ErMat[:,i] \
                            - ((lcm/headrow[0]) * ErMat[:,n]).astype(int)
        return ErMat.astype(int)
    ermat = np.concatenate([np.eye(N), mat.copy()], axis = 0).astype(int)
    # Computation for a basis of the kernel
    nonzeroidx = []
    for n in range(N):
        mlist = (np.arange(M)[np.any(ermat[N:,n:], axis = 1)])
        if len(mlist) == 0:
            break
        ermat = _gaussian(ermat, mlist[0], n)
        nonzeroidx.append([N+mlist[0], n])
    for idx in nonzeroidx:
        if (ermat[idx[0], idx[1]] < 0):
            ermat[:, idx[1]] = -ermat[:, idx[1]]
    indice = ~np.any(ermat[N:,:], axis = 0)
    return ermat[:N, indice], ermat[N:, ~indice]

def _AugMat(nu, Xname, mask):
    """
    To construct a symbolic augmented matrix A given a stoichiometric matrix.
    
    [Input]
        nu      (2D array-like)
                the stoichiometric matrix, of which entries are required to be integers.
        Xname   (1D iterable)
                the name of the system variables.
        mask    (2D array-like)
                pecifies the condition in which partial derivatives of r can be omitted.
    Output:
        spA
            the symbolic augmented matrix A.
    """
    col_stand = lambda mat: (mat / np.gcd.reduce(mat, axis = 0))
    C   = col_stand(__KerImg(nu)[0])
    D   = col_stand(__KerImg(nu.T)[0])
    dim = nu.shape[0]+C.shape[1]
    spA = sp.Matrix(np.zeros([dim, dim]).astype(int))
    if np.all(C.shape):
        spA[:C.shape[0],-C.shape[1]:] = sp.Matrix(C.astype(int))
    if np.all(D.shape):
        spA[-D.shape[1]:,:D.shape[0]] = sp.Matrix(D.T.astype(int))
    for i in range(nu.shape[0]):
        for j in range(nu.shape[1]):
            if not mask[i, j]:
                spA[j,i] = sp.symbols('r_{%s%s}'%(j+1, Xname[i]))
    return spA


##########################################################################
# Displaying by manipulatiung a browser to read HTML with LaTeX embedded #
##########################################################################

from selenium import webdriver
import os, time

_mathjax_  = """<!DOCTYPE html>
<html>
<head>
<script type="text/javascript" id="MathJax-script" async
src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js"></script>
</head>
<body>
\[
    {}
\]
</body>
</html>
""".format

def _set_browser(browsername):
    "To set the browser for displaying MathML"
    if browsername == 'Firefox':
        browser = webdriver.Firefox()
    elif browsername == 'Chrome':
        browser = webdriver.Chrome()
    elif browsername == 'Safari':
        browser = webdriver.Safari()
    elif browsername == 'ChromiumEdge':
        browser = webdriver.ChromiumEdge()
    return browser

class newsheet(object):

    """
    To open a browswer reading _scratchpaper_.html so as to display results.
    Always remember to .close() to dispose the data.
    """

    def __init__(self, mode = "w", browsername = 'Firefox'):
        self._driver  = _set_browser(browsername)
        self._latex   = ""
        self._mode    = mode

    @property
    def mode(self):
        return self._mode

    @mode.setter
    def mode(self, newmode = "w"):
        self._mode = newmode

    def close(self):
        """
        To close the browser and clean the _scratchpaper_.html.
        """
        self._driver.close()
        os.remove("_scratchpaper_.html")

    def update(self, contents, end = ""):
        """
        The input (contents) is assumed to be a string that reads in LaTeX.
        """
        if self._mode == 'w':
            self._latex = ""
        self._latex = self._latex + contents + end
        with open("_scratchpaper_.html", "w") as scppml:
            scppml.write(_mathjax_(self._latex))
            self._driver.get("file://"+os.path.realpath(scppml.name))
        time.sleep(1)
        self._driver.refresh()

    def display(self, spobj, end = ""):
        """
        The input (spobj) is assumed to be a sympy object
        """
        self.update(sp.latex(spobj), end)

## == END OF THE SCRIPT == ##
