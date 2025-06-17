from sage.all import log
import numpy as np
from .ckks_in_sagemath.ckks_package.bit_rev import bit_rev


def ext_bit_rev(i, B, log_B=None):
    """
    Perform extended bit-reversing of the integer i with respect to B.

    Args:
        i (int):
            The integer to reverse.
        B (int):
            The power of two with respect to which the extended bit-reversing
            is computed.
        log_B (int, optional):
            The logarithm of B to the base 2.

    Returns:
        int:
            The bit-reversed integer.
    """
    if log_B is None:
        log_B = log(B, 2)
    q, r = divmod(i, B)
    return q * B + bit_rev(r, log_B)


def ext_bit_rev_vector(v, B):
    """
    Rearrange the entries of a vector into extended bit-reversed order.

    Args:
        v (np.ndarray):
            The input vector whose entries should be rearranged.
        B (int):
            The power of two with respect to which the extended bit-reversing
            is computed.

    Returns:
        np.ndarray:
            A new vector with entries rearranged in extended bit-reversed
            order.
    """
    n = len(v)
    log_B = log(B, 2)
    return np.array([v[ext_bit_rev(i, B, log_B)] for i in range(n)])
