import networkx as nx
import numpy as np
import numpy as np
from scipy.sparse import issparse
from shapely.geometry import Polygon,Point, MultiPolygon
import os
import gzip
import pickle
import subprocess  # Start process for git hash
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

sx = np.array([[0, 1], [1, 0]])
sy = np.array([[0, -1.j], [1.j, 0]])
sz = np.array([[1, 0], [0, -1]])
s0 = np.eye(2)


def draw_sphere(r,ax,N=50,stride=1,**kwargs):
    u = np.linspace(0, 2 * np.pi, N)
    v = np.linspace(0, np.pi, N)
    x = r*np.outer(np.cos(u), np.sin(v))
    y = r*np.outer(np.sin(u), np.sin(v))
    z = r*np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(x, y, z, linewidth=0.0, cstride=stride, rstride=stride,**kwargs)


# ========== Convenience functions to handle states ====================

def generate_state(theta, phi):
    """Generate the state
    |\Psi> = cos(theta/1) |0> + exp(i*phi) sin(theta/2) |1>.

    Args:
        theta (float): polar angle
        phi (float): azimutal angle 

    Returns:
        np.ndarray: State |Psi> on the Bloch sphere
    """
    return np.cos(theta/2.) * np.array([1, 0]) + np.exp(1.j*phi)*np.sin(theta/2)*np.array([0, 1])


def generate_state_2qubit(tp1, tp2):
    """Generate a two-qubit state characterized by two theta-phi tuples

    Args:
        tp1 (tuple): (theta1, phi1)
        tp2 (tuple): (theta2, phi2)

    Returns:
        np.ndarray: 2qubit state (4d vector)
    """
    ket0 = generate_state(*tp1)
    ket1 = generate_state(*tp2)
    return np.kron(ket0, ket1)


# ========== Utility Functions ====================

def draw_coordinate_graph(g, flip_theta_phi = False):
    """Draw a graph given by g

    Args:
        g (nx.Graph): Input graph
    """
    if flip_theta_phi:
        pos = {ind:data["coord"][::-1] for ind,data in g.nodes(data=True)}
    else:
        pos = {ind:data["coord"] for ind,data in g.nodes(data=True)}
    nx.draw_networkx(g, pos)
    labels = nx.get_edge_attributes(g, 'weight')
    nx.draw_networkx_edge_labels(g, pos, edge_labels=labels)
    if flip_theta_phi:
        plt.xlabel(r"$\phi$")
        plt.ylabel(r"$\theta$")
    else:
        plt.ylabel(r"$\phi$")
        plt.xlabel(r"$\theta$")


def draw_graph(g):
    """Draw a graph given by g

    Args:
        g (nx.Graph): Input graph
    """
    pos = nx.random_layout(g)
    nx.draw_networkx(g, pos)
    labels = nx.get_edge_attributes(g, 'weight')
    nx.draw_networkx_edge_labels(g, pos, edge_labels=labels)


def subplot_grid(n: int):
    """Generate a grid of n axis as close as possible to a square.

    Args:
        n (int): number of plots

    Returns:
        (int,int),figure,axisvec: (ncols,nrows), figure, axisvec
    """
    nearest_square = int(np.sqrt(n))
    diff = n-nearest_square**2
    if diff > 0:
        nrows = nearest_square
        ncols = nearest_square+int(np.ceil(diff/nrows))
    else:
        nrows = nearest_square
        ncols = nearest_square
    f, ax = plt.subplots(nrows, ncols)
    if len(ax.shape) == 1:
        # This is 1 d array, but we expect a grid
        ax = ax[..., np.newaxis]
    return (ncols, nrows), f, ax


def base_repr(val, base, padding=0):
    base_str = np.base_repr(val, base)
    base_str = base_str.zfill(padding)
    dest = [int(x) for x in base_str]
    return dest



def print_columns(listvals, padding=4, header=False):
    """Print a multi-dimensional list in a table

    Args:
        listvals (list of lists): Input data
        padding (int, optional): Padding of the columns. Defaults to 4.
        header (bool, optional): Print a header on top of the table. Defaults to False.
    """
    col_width = max([len(str(word))
                     for row in listvals for word in row]) + padding
    for ind, row in enumerate(listvals):
        print("".join(str(word).ljust(col_width) for word in row))
        if header and ind == 0:
            print("")


def sizeof_fmt(num, suffix='B'):
    """Pretty print a size as mutliples of 1024."""
    for unit in ['', 'Ki', 'Mi', 'Gi', 'Ti', 'Pi', 'Ei', 'Zi']:
        if abs(num) < 1024.0:
            return "%3.1f %s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%3.1f %s%s" % (num, 'Yi', suffix)


def get_git_hash():
    """Get the git hash of the current commit in the repository.

    Returns:
        str: git hash
    """
    # This assumes that .git is two levels above of util.py
    packagedir = os.path.join(os.path.dirname(
        os.path.realpath(__file__)), os.path.pardir)
    rootdir = os.path.join(packagedir, os.path.pardir)
    gitdir = os.path.join(rootdir, ".git")
    githash = subprocess.check_output(
        ['git', '--git-dir={}'.format(gitdir), 'rev-parse', 'HEAD'])
    return githash.decode("utf-8").strip()


def select_except(arr, ind: int):
    """Return all elements of a list except the indicated one

    Args:
        arr (list/np.array): list of values
        ind (int): index

    Returns:
        np.array: Array with all elements of arr except for arr[ind]
    """
    # This function works only on the outer-most layer
    if isinstance(arr, list):
        arr = np.asarray(arr)
    mask = np.ones(len(arr), dtype=bool)
    mask[ind] = False
    return arr[mask]


def multiply_except(arr, ind: int):
    """Product of all array values except for arr[ind]

    Args:
        arr (list/np.arr): list of values
        ind (int): index

    Returns:
        float: Multiplication of all array values except for arr[ind]
    """
    if len(arr) > 1:
        others = select_except(arr, ind)
        return np.prod(others)
    else:
        # It does not make sense to execute this function with only one element
        return arr[0]


# =========== Matrix Evaluation Functions ====================

def is_hermitian(mat):
    """Returns true if the matrix is hermitian."""
    if issparse(mat):
        return np.allclose(mat.todense(), mat.H.todense())
    else:
        return np.allclose(np.conjugate(np.transpose(mat)), mat)


def is_diagonal(mat):
    """Returns true if the matrix is diagonal."""
    if issparse(mat):
        return np.allclose((mat-mat.diagonal()).todense(), np.zeros(mat.shape))
    else:
        return np.allclose(mat-np.diag(np.diag(mat)), np.zeros_like(mat))


def is_symmetric(mat):
    """Returns true if the matrix is symmetric. """
    if issparse(mat):
        return np.allclose(mat.todense(), mat.T.todense())
    else:
        return np.allclose(np.transpose(mat), mat)


def is_permutation(mat):
    """Returns true if the matrix is a permutation matrix. """
    n, m = mat.shape
    if issparse(mat):
        raise NotImplementedError(
            "Checking for sparse permutation matrices is not implemented.")
    else:
        square = n == m
        id = np.allclose(np.eye(n), mat@np.transpose(mat))
        sum_rows = np.all(np.sum(mat, axis=0) == 1)
        sum_cols = np.all(np.sum(mat, axis=1) == 1)
        return square and id and sum_rows and sum_cols


def is_antisymmetric(mat):
    """Returns true if the matrix is symmetric. """
    if issparse(mat):
        return np.allclose(mat.todense(), -mat.T.todense())
    else:
        return np.allclose(-np.transpose(mat), mat)


def anti_symmetrize(mat):
    """Force a matrix to be anti-symmetirc."""
    if issparse(mat):
        return 0.5*(mat-mat.T)
    else:
        return 0.5*(mat-np.transpose(mat))


def get_nonzero_fraction(mat):
    """Returns fraction of non-zero elements."""
    return np.count_nonzero(mat)/np.prod(mat.shape)


def herm_conj(mat):
    """Returns the hermitian conjugate of a matrix."""
    return np.conjugate(np.transpose(mat))


def commutator(mat1, mat2):
    """Calculate the commutator of two matrices

    Args:
        mat1 (2d np.ndarray): First argument of commutator
        mat2 (2d np.ndarray): Second argument of commutator

    Returns:
        2d np.ndarray: Commutator
    """
    return mat1@mat2-mat2@mat1


def anticommutator(mat1, mat2):
    """Calculate the anti-commutator of two matrices

    Args:
        mat1 (2d np.ndarray): First argument of anti-commutator
        mat2 (2d np.ndarray): Second argument of anti-commutator

    Returns:
        2d np.ndarray: Anti-commutator
    """
    return mat1@mat2+mat2@mat1


# =========================== Cache Server =================================


class CacheServer:
    """Storage Server for arbitrary data that can be stored in dictionaries"""

    def __init__(self):
        self.store = {}

    def add(self, name, mat):
        self.store[name] = mat

    def get(self, name):
        try:
            return self.store[name]
        except KeyError:
            return None

    def load(self, fname):
        if os.path.isfile(fname):
            with gzip.open(fname, "rb") as infile:
                self.store = pickle.load(infile)

    def save(self, fname):
        # We only save if the file does not exist yet
        if not os.path.isfile(fname):
            with gzip.open(fname, "wb") as outfile:
                pickle.dump(self.store, outfile)

    def list(self):
        print(self.store.keys)

    def __str__(self):
        print("CacheServer: {} Entries".format(len(self.store)))


# ========== Debugging Functions ====================

def show_vector(vec, title=None):
    """Display a matrix and interrupt the program. """
    f, ax = plt.subplots(1, 1)
    ax.plot(vec)
    if title is not None and len(title) > 0:
        plt.title(title)
    plt.show()


def show_matrix(mat, title=None, **kwargs):
    """Display a matrix and interrupt the program. """
    show_matrixvec([mat], title=[title], **kwargs)


def show_matrixvec(matvec, title=None, log=False):
    """Display a matrix and interrupt the program. """
    f, axvec = plt.subplots(1, len(matvec))
    if len(matvec) == 1:
        axvec = [axvec]
    for ind, mat in enumerate(matvec):
        if log:
            minval = np.min(mat)
            if minval == 0:
                # This is a dirty hack to display the 0 in a log plot
                mat += 1e-10
                minval += 1e-10
            matax = axvec[ind].matshow(
                mat, norm=LogNorm(vmin=minval, vmax=np.max(mat)))
        else:
            matax = axvec[ind].matshow(mat)
        f.colorbar(matax, ax=axvec[ind])

    if title is not None:
        if type(title) is list and len(title) == len(matvec):
            for ax, titleval in zip(axvec, title):
                ax.set_title(titleval)
        elif type(title) is str and len(title) > 0:
            plt.title(title)
    plt.show()


def print_mat_stats(mat, title=None):
    """Display general information about matrix."""
    print("Min:\t{}".format(np.min(mat)))
    print("Max:\t{}".format(np.max(mat)))
    print("Avg:\t{}".format(np.mean(mat)))
    print("Norm:\t{}".format(np.linalg.norm(mat)))


def show_eigenvalues(mat):
    """Display the eigenvalues of a matrix"""
    if is_hermitian(mat):
        # Plot the real eigenvalues
        f, ax = plt.subplots(1, 1)
        eigvals = np.linalg.eigvalsh(mat)
        ax.plot(eigvals, 'o')
    else:
        # Plot the real eigenvalues
        f, ax = plt.subplots(1, 2)
        eigvals = np.linalg.eigvals(mat)
        ax[0].set_title("Real part")
        ax[0].plot(np.real(eigvals), 'o')
        ax[0].set_title("Imaginary part")
        ax[1].plot(np.imag(eigvals), 'o')
    plt.show()


def show_polygon(multipolygon:MultiPolygon):
    f,ax = plt.subplots()
    from descartes import PolygonPatch
    for polygon in multipolygon:
        patch = PolygonPatch(polygon, alpha=0.5, zorder=2)
        ax.add_patch(patch)
        ax.set_xlabel(r"$\phi$")
        ax.set_ylabel(r"$\theta$")
        ax.set_xlim(0,2*np.pi)
        ax.set_ylim(0,np.pi)
    plt.show()

# ========== Testing Functions ====================

def compare_array_elementwise(testcase, ref, res, print_vals=True):
    testcase.assertEqual(ref.shape, res.shape)
    if print_vals:
        for i in range(ref.shape[0]):
            for j in range(ref.shape[1]):
                if not np.isclose(ref[i, j], res[i, j]):
                    print("{},{}: ref: {},res:{}".format(
                        i, j, ref[i, j], res[i, j]))
    testcase.assertTrue(np.allclose(ref, res))
