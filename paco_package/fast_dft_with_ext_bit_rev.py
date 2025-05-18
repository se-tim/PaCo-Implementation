from sage.all import log
import numpy as np
from .ckks_in_sagemath.ckks_package.bit_rev import bit_rev
from .ckks_in_sagemath.ckks_package.fast_dft import (
    get_E,
    get_roots_of_unity,
    Multidiags,
)


def get_matrix_D(n, l, inverse=False):
    """
    Compute the matrix D_{n, 2^l} (or its inverse) representing the ring
    isomorphism from R_{n, 2^{l+1}} to R_{n, 2^l}.

    Args:
        n (int):
            The size of the matrix (must be a power of two).
        l (int):
            The index of the matrix, a value in [0, log(n, 2) - 1].
        inverse (bool, optional):
            If True, computes the inverse of the matrix. Defaults to False.

    Returns:
        Multidiags:
            The matrix D_{n, 2^l} or its inverse.
    """
    return get_E(n, l, inverse)  # Not the matrix E from the paper!


def get_matrix_E(n, l, B, inverse=False):
    """
    Compute the matrix E_{n, 2^l}^{(B)} (or its inverse) representing the ring
    isomorphism from R_{n, 2^{l+1}} to R_{n, 2^l}, but with interposed
    permutation matrices corresponding to extended bit-reversing with respect
    to B.

    Args:
        n (int):
            The size of the matrix (must be a power of two).
        l (int):
            The index of the matrix, a value in [0, log(n, 2) - 1].
        B (int):
            The power of two with respect to which the extended bit-reversing
            is computed.
        inverse (bool, optional):
            If True, computes the inverse of the matrix. Defaults to False.

    Returns:
        Multidiags:
            The matrix E_{n, 2^l}^{(B)} or its inverse.
    """
    half_k = B // (2 * 2**l)

    if half_k < 1:
        return get_E(n, l, inverse)

    k = 2 * half_k

    roots = get_roots_of_unity(4 * n)
    half_roots = roots / 2

    n_over_B = n // B
    log_n_over_B = log(n_over_B, 2)

    diag0 = np.zeros(n, dtype=np.complex128)
    diag1 = np.zeros(n, dtype=np.complex128)
    diag2 = np.zeros(n, dtype=np.complex128)

    if inverse == False:
        for i in range(n):
            i_div, i_mod = divmod(i, k)
            if i_mod < half_k:
                index = (
                    (2**l)
                    * 5
                    ** (
                        (
                            i_mod * n_over_B
                            + bit_rev(i_div // 2**l, log_n_over_B)
                        )
                    )
                ) % (4 * n)
                diag0[i] = roots[0]
                diag1[i] = roots[index]
            else:
                index = (
                    (2**l)
                    * 5
                    ** (
                        (
                            (i_mod - half_k) * n_over_B
                            + bit_rev(i_div // 2**l, log_n_over_B)
                        )
                    )
                ) % (4 * n)
                diag0[i] = -roots[index]
                diag2[i] = roots[0]

    else:
        for i in range(n):
            i_div, i_mod = divmod(i, k)
            if i_mod < half_k:
                diag0[i] = half_roots[0]
                diag1[i] = half_roots[0]
            else:
                index = (
                    (2**l)
                    * 5
                    ** (
                        (i_mod - half_k) * n_over_B
                        + bit_rev(i_div // 2**l, log_n_over_B)
                    )
                ) % (4 * n)
                diag0[i] = -half_roots[-index]
                diag2[i] = half_roots[-index]

    if half_k == n // 2:
        return Multidiags(n, {0: diag0, half_k: diag1 + diag2})
    else:
        return Multidiags(n, {0: diag0, half_k: diag1, -half_k: diag2})
