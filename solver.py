import numpy as np
import numpy.typing as npt

from scipy.sparse import csr_matrix, isspmatrix
from scipy.sparse.linalg import spsolve


def solve(
    K: csr_matrix | npt.NDArray[np.float64],
    F: npt.NDArray[np.float64],
    u_fixed_idx: list[int],
) -> npt.NDArray[np.float64] | None:
    """
    Lösen des linearen Systems KU = F 
    Parameter
    ----------
    K : csr_matrix 
        Global stiffness matrix 
    F : ndarray
        Global force vector.
    u_fixed_idx : list

    gibt U zurück oder None wenn das System singulär ist.
    -------
    u : ndarray or None
       Verchiebungsvektor oder None wenn das System singulär ist.
    """

    n = F.shape[0]
    assert K.shape[0] == K.shape[1] == n

    # Umwandeln in sparse matrix falls nötig
    if not isspmatrix(K):
        K = csr_matrix(K)

    fixed = np.array(u_fixed_idx, dtype=int)
    free = np.setdiff1d(np.arange(n), fixed)

    if free.size == 0:
        return np.zeros(n)

    # System für freie Freiheitsgrade lösen
    K_ff = K[free, :].tocsc()[:, free].tocsr()
    F_f = F[free]

    try:
        u_f = spsolve(K_ff, F_f)
    except Exception:
        return None

    # Vollständigen Verschiebungsvektor erstellen
    u = np.zeros(n)
    u[free] = u_f
    u[fixed] = 0.0

    return u



def test_case_horizontal():
    # Horizontal Element in x-Richtung
    e_n = np.array([1.0, 0.0])
    e_n /= np.linalg.norm(e_n)

    k = 1.0
    K_local = k * np.array([[1.0, -1.0], [-1.0, 1.0]])

    O = np.outer(e_n, e_n)
    K_global = np.kron(K_local, O)

    K_global = csr_matrix(K_global)

    u_fixed_idx = [0, 1]  # node i fixiert
    F = np.array([0.0, 0.0, 10.0, 0.0])

    u = solve(K_global, F, u_fixed_idx)
    print("Horizontal test:", u)


def test_case_diagonal():
    # Diagonales element in 45° Richtung
    e_n = np.array([1.0, 1.0])
    e_n /= np.linalg.norm(e_n)

    k = 1.0 / np.sqrt(2.0)
    K_local = k * np.array([[1.0, -1.0], [-1.0, 1.0]])

    O = np.outer(e_n, e_n)
    K_global = np.kron(K_local, O)

    K_global = csr_matrix(K_global)

    u_fixed_idx = [0, 1]
    F = np.array([0.0, 0.0, 1.0, 1.0])

    u = solve(K_global, F, u_fixed_idx)
    print("Diagonal test:", u)


if __name__ == "__main__":
    test_case_horizontal()
    test_case_diagonal()