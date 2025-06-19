from sage.all import log, randint, ZZ
import numpy as np
import time
from .ckks_in_sagemath.ckks_package.ckks import CKKS
from .ckks_in_sagemath.ckks_package.fast_dft import group_matrices
from .ckks_in_sagemath.ckks_package.poly import Poly
from .ext_bit_rev import ext_bit_rev_vector
from .fast_dft_with_ext_bit_rev import get_matrix_D, get_matrix_E


class CKKS_x(CKKS):
    """
    Class for CKKS scheme with PaCo bootstrapping.
    """

    @classmethod
    def config(cls, N, L_boot, q0, p, delta, print_messages=False):
        """
        Configure basic CKKS parameters.

        Args:
            N (int):
                Ring dimension (must be a power of two).
            L_boot (int):
                Maximal level during bootstrapping.
            q0 (int):
                Smallest modulus.
            p (int):
                Scaling factor outside of bootstrapping.
            delta (int):
                Scaling factor during bootstrapping.
            print_messages (bool, optional):
                Whether to print information messages. Defaults to False.
        """
        super().config(
            N, N // 2, L_boot, q0, p, delta, print_messages=print_messages
        )

    @classmethod
    def key_gen(cls, h=None, print_messages=False):
        """
        Generate the secret key, public key and evaluation key for the scheme.

        Args:
            h (int, optional):
                Hamming weight of the secret key. Must be a power of two.
                Defaults to 2 ** min(6, (log(N, 2) // 2)).
            print_messages (bool, optional):
                Whether to print information messages. Defaults to False.
        """
        if h is not None and h & (h - 1) != 0:
            raise ValueError("h must be a power of two.")
        if h is None:
            h = 2 ** min(6, 2 ** (log(cls.N, 2) // 2))

        # Generate the PaCo secret key (skGen in paper)

        B = cls.N // (4 * h)
        u = [0] + [randint(0, B - 1) for _ in range(4 * h - 1)]
        d = [1] + [0] * (4 * h - 1)
        for v in range(1, h):
            t = randint(0, 3)
            d[t * h + v] = 1
        cls.u = u
        cls.d = d
        sk_coeffs = [0] * cls.N
        for v in range(4 * h):
            sk_coeffs[v + u[v] * 4 * h] = d[v]
        sk = Poly(sk_coeffs, cls.N)

        q = cls.moduli_boot[-2]  # Largest modulus for evaluation key
        P = cls.moduli_boot[-2]  # Extra factor added to evaluation key modulus
        super().key_gen(sk=sk, q=q, P=P, print_messages=print_messages)

    @classmethod
    def get_security(cls, check_primal_hybrid=False):
        """
        Use the LWE estimator to compute the security level.

        Args:
            check_primal_hybrid (bool, optional):
                Whether to check the primal hybrid security level. Defaults to
                False.

        Returns:
            float:
                Logarithm (base 2) of estimated number of operations required
                to break the scheme.
        """
        return super().get_security(
            N=cls.N - cls.N // cls.h,
            q=cls.moduli_boot[-2] ** 2,
            sk_plus=cls.h - 1,
            sk_minus=0,
            check_primal_hybrid=check_primal_hybrid,
        )

    @classmethod
    def config_PaCo(cls, sk, C, g0, g1, print_messages=False):
        """
        Perform precomputations needed for PaCo bootstrapping.

        Args:
            sk (Poly):
                Secret key.
            C (int):
                Number of coefficients to bootstrap. Must be a power of two and
                at least 2.
            g0 (int):
                Grouping parameter for partial CoeffToSlot.
            g1 (int):
                Grouping parameter for conventional SlotToCoeff.
            print_messages (bool, optional):
                Whether to print information messages. Defaults to False.
        """
        if 4 * C * cls.h > cls.N:
            raise ValueError("4 * C * h is bounded by N.")
        if C & (C - 1) != 0:
            raise ValueError("C must be a power of two.")
        if C < 2:
            raise ValueError("C must be at least 2.")

        cls._check_key_gen()

        # Generate the bootstrapping keys (bskGen in paper)

        if print_messages:
            print("Generating bootstrapping keys...")
        cls.C = C
        k = cls.N // (4 * cls.h * C)
        n = 2 * cls.h * C
        bsk_list = []
        for t in range(4):
            sigma_t = []
            for r in range(k):
                z = np.zeros(n).astype(np.complex128)
                for v in range(cls.h):
                    i = cls.u[t * cls.h + v] + r
                    if i % k == 0:
                        z[v * 2 * cls.C + i // k] = cls.d[t * cls.h + v]
                for l in range(log(cls.C, 2) + 1):
                    z = get_matrix_D(n, log(cls.C, 2) - l) * z
                sigma_t.append(z)
            sigma_t = np.concatenate(sigma_t)
            rev_sigma_t = ext_bit_rev_vector(sigma_t, cls.C // 2)
            bsk_list.append(
                cls.enc_poly_with_sk(
                    cls.encode(rev_sigma_t, is_boot=True), sk, is_boot=True
                )
            )
        cls.bsk = bsk_list

        # Generate and encode vectors mu and eta

        if print_messages:
            print("Encoding relevant vectors as polynomials...")
        mu = np.array([(i % (2 * cls.C) < cls.C) for i in range(n)])
        eta = (
            -cls.q0
            * 1j
            / (4 * cls.delta * np.pi)
            * np.array(
                [
                    (i < cls.C // 2) + 1j * (i >= cls.C // 2)
                    for i in range(cls.C)
                ]
            )
        )
        cls.encoded_mu = cls.encode(mu, is_boot=True)
        cls.encoded_eta = cls.encode(eta, is_boot=True)

        # Generate and encode matrices for partial CoeffToSlot and conventional
        # SlotToCoeff

        if print_messages:
            print(
                "Generating and encoding matrices required for partial "
                "CoeffToSlot and conventional SlotToCoeff..."
            )
        g0 = min(g0, log(cls.C, 2) + 1)
        g1 = min(g1, log(cls.C, 2) - 1)
        iE_list = [
            get_matrix_E(n, l, cls.C // 2, True)
            for l in range(log(cls.C, 2) + 1)
        ]
        E_list = [
            get_matrix_E(cls.C // 2, l, cls.C // 2)
            for l in range(log(cls.C, 2) - 1)
        ]
        grouped_iE_list = group_matrices(iE_list, g0, True)
        grouped_E_list = group_matrices(E_list, g1)
        cls.grouped_poly_iE_list = [
            cls.get_poly_matrix(A, is_boot=True) for A in grouped_iE_list
        ]
        cls.grouped_poly_E_list = [
            cls.get_poly_matrix(A, is_boot=True) for A in grouped_E_list
        ]

        # Generate missing switching keys

        if print_messages:
            print("Generating missing switching keys...")
        cls.get_galois_swk(
            -1, sk, cls.moduli_boot[-2], cls.P
        )  # For conjugation
        cls.get_galois_swk(
            -(5**cls.C), sk, cls.moduli_boot[-2], cls.P
        )  # For conjugation with rotation by C slots
        rotation_indices = [
            2**i for i in range(log(cls.N, 2))
        ]  # For rotation by powers of two
        for poly_matrix in cls.grouped_poly_iE_list + cls.grouped_poly_E_list:
            rotation_indices += cls.get_BSGS_rotation_indices(poly_matrix)
        for i in rotation_indices:
            cls.get_galois_swk(5**i, sk, cls.moduli_boot[-2], cls.P)

        if print_messages:
            print("The PaCo bootstrapping configuration is done!")

    def get_coeff_encodings(self):
        """
        Generate the coefficient encodings (getCoeffEnc in paper).

        Returns:
            list:
                List of coefficient encodings from the ciphertext self.

        """
        k = self.N // (4 * self.h * self.C)
        n = 2 * self.h * self.C

        def psi(a):
            return np.exp(2 * np.pi * 1j * ZZ(a) / self.q0)

        self = self % self.q0
        ct0, ct1 = self.b.coeffs, self.a.coeffs

        tilde_b_dict = {}

        for v in range(4 * self.h):
            for r in range(k):
                tilde_b_v_r = np.zeros(2 * self.C).astype(np.complex128)
                for i in range(self.C):
                    if v == 0:
                        index = 4 * self.h * (i * k + r)
                        tilde_b_v_r[i] = psi(ct0[index] + ct1[index])
                    elif i == 0 and r == 0:
                        tilde_b_v_r[i] = psi(-ct1[self.N - v])
                    else:
                        tilde_b_v_r[i] = psi(ct1[4 * self.h * (i * k + r) - v])
                tilde_b_dict[(v, r)] = tilde_b_v_r

        coeff_encoding_list = []
        for t in range(4):
            beta_t = []
            for r in range(k):
                z = [tilde_b_dict[(t * self.h + v, r)] for v in range(self.h)]
                z = np.concatenate(z)
                for l in range(log(self.C, 2) + 1):
                    # These matrices have already been precomputed for the bsk
                    z = get_matrix_D(n, log(self.C, 2) - l) * z
                beta_t.append(z)
            beta_t = np.concatenate(beta_t)
            rev_beta_t = ext_bit_rev_vector(beta_t, self.C // 2)
            coeff_encoding_list.append(self.encode(rev_beta_t, is_boot=True))
        return coeff_encoding_list

    def seq_PaCo(self):
        """
        Apply the PaCo bootstrapping procedure (seqPaCo in paper). The
        ciphertext self is refreshed by increasing its level. Only the
        coefficients indexed by multiples of N / C are bootstrapped.

        Returns:
            CKKS:
                A new ciphertext with increased level.
        """
        n = 2 * self.h * self.C
        coeff_encoding_list = self.get_coeff_encodings()

        ct = sum(coeff_encoding_list[t] @ self.bsk[t] for t in range(4))
        ct = ct.trace(self.N // 2, n)
        for iE in self.grouped_poly_iE_list:
            # Partial CtS
            ct = ct.BSGS_left_mult(iE)
        ct_x = ct.conj_rotate(self.C)
        ct = ct + ct_x
        ct = self.encoded_mu @ ct
        ct = ct.product(n, 2 * self.C)
        ct = ct.trace(2 * self.C, self.C)
        ct = ct - ct.conjugate()
        ct = self.encoded_eta @ ct
        ct = ct.trace(self.C, self.C // 2)
        for E in self.grouped_poly_E_list:
            # Decomposed StC
            ct = ct.BSGS_left_mult(E)
        ct = ct.boot_to_nonboot()
        return ct

    def parallel_PaCo(self, kappa, get_time=False):
        """
        Apply the PaCo bootstrapping procedure, but kappa times “in parallel”
        (parallelPaCo in paper). Only the coefficients indexed by multiples of
        N / (kappa * C) are bootstrapped.

        Args:
            kappa (int):
                Number of processors.
            get_time (bool, optional):
                Whether to return the time for a parallel simulation.
                Defaults to False.

        Returns:
            CKKS:
                A new ciphertext with increased level.
        """
        ct_list = []
        times = []
        for r in range(kappa):
            # Each iteration of this loop can be executed independently in
            # parallel
            time_start = time.time()
            monomial0 = Poly.get_monomial(
                -r * self.N // (self.C * kappa), self.N, self.q0
            )
            monomial1 = Poly.get_monomial(
                r * self.N // (self.C * kappa), self.N, self.moduli_boot[1]
            )
            ct = monomial0 * self
            ct = ct.seq_PaCo()
            ct = monomial1 * ct
            ct_list.append(ct)
            times.append(time.time() - time_start)
        parallel_time = max(times)
        time_start = time.time()
        ct = sum(ct_list)
        parallel_time += time.time() - time_start
        if get_time:
            return ct, parallel_time
        return ct

    def get_precision(self, other, sk, num_coeffs=None):
        """
        Compute the average number of bits of precision preserved in the
        coefficients of the underlying plaintext polynomial of other, relative
        to self.

        Args:
            other (CKKS):
                The ciphertext to compare against self.
            sk (Poly):
                The secret key.
            num_coeffs (int, optional):
                Number of coefficients to consider. Defaults to C.

        Returns:
            float:
                The average number of bits of precision preserved from self to
                other.

        """
        if num_coeffs is None:
            num_coeffs = self.C
        return super().get_precision(other, sk, num_coeffs // 2)
