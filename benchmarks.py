from sage.all import ceil, log, randint
import time
from paco_package.ckks_x import CKKS, CKKS_x, Poly

# Configuration

print(
    "This script benchmarks the PaCo bootstrapping against the original one."
)

parameter_sets = [
    [2**15, 2**29, 2**22, 2**32, 64, 3, 7, 8],  # Set 1
    [2**16, 2**44, 2**32, 2**48, 64, 3, 7, 9],  # Set 2
]  # Under the form [N, q, p, delta, h, g, d, r]

if input("Use a parameter set from the paper (y / n)? ") == "y":
    set_number = int(input("Which parameter set (1 / 2)? "))
    N, q, p, delta, h, g, d, r = parameter_sets[set_number - 1]
else:
    N = 2 ** int(input("Ring degree: N = 2^"))
    q = 2 ** int(input("Base modulus: q = 2^"))
    p = 2 ** int(input("Scaling factor outside of bootstrapping: p = 2^"))
    delta = 2 ** int(input("Scaling factor during bootstrapping: delta = 2^"))
    h = 2 ** int(input("Hamming weight of secret key: h = 2^"))
    g = int(input("Grouping parameter for CoeffToSlot and SlotToCoeff: g = "))
    d = int(input("Degree for Taylor polynomial in EvalMod: d = "))
    r = int(input("Number of successive squarings in EvalMod: r = "))
num_tests = int(input("Number of tests for each configuration: "))

B = N // (4 * h)
log_N = log(N, 2)
log_h = log(h, 2)
log_d = ceil(log(d, 2))
log_B = log(B, 2)

if input("Test sequential or parallel bootstrapping (s / p)? ") == "s":
    print()
    header = (
        f" {'C':^5} || "
        + " | ".join(
            [f"{'PaCo try ' + str(i+1):^12}" for i in range(num_tests)]
        )
        + " || "
        + " | ".join(
            [f"{'Orig try ' + str(i+1):^12}" for i in range(num_tests)]
        )
    )
    print(header)
    print("-" * len(header))

    for log_C in range(1, log_B + 1):
        C = 2**log_C
        print(f" {C:^5} ", end="|", flush=True)

        # PaCo
        L = ceil((log_C + 1) / g) + ceil((log_C - 1) / g) + log_h + 3
        CKKS_x.config(N, L, q, p, delta)
        CKKS_x.key_gen(h)
        CKKS_x.config_PaCo(CKKS_x.sk, C, g, g)
        precision_paco = 0
        for _ in range(num_tests):
            coeffs = [
                randint(-p, p) if i % (N // C) == 0 else 0 for i in range(N)
            ]
            m = Poly(coeffs, N)
            ct = CKKS_x.enc_poly_with_sk(m, CKKS_x.sk) % q
            t_start = time.time()
            ct_boot = ct.seq_PaCo()
            t_end = time.time()
            t = f"{float(round(t_end - t_start, 2))} s"
            print(f"| {t:^12} ", end="", flush=True)
            if log_C == log_B:
                precision_paco += (
                    ct.get_precision(ct_boot, CKKS_x.sk) / num_tests
                )
        print("|", end="", flush=True)

        # Original
        L = 2 * ceil(max(log_C - 1, 1) / g) + log_d + r + 1
        CKKS.config(N, 2 ** (log_C - 1), L, q, p, delta)
        CKKS.key_gen(h)
        CKKS.config_bootstrap(CKKS.sk, d, r, g)
        precision_orig = 0
        for _ in range(num_tests):
            coeffs = [
                randint(-p, p) if i % (N // C) == 0 else 0 for i in range(N)
            ]
            m = Poly(coeffs, N)
            ct = CKKS.enc_poly_with_sk(m, CKKS.sk) % q
            t_start = time.time()
            ct_boot = ct.bootstrap()
            t_end = time.time()
            t = f"{float(round(t_end - t_start, 2))} s"
            print(f"| {t:^12} ", flush=True)
            if log_C == log_B:
                precision_orig += (
                    ct.get_precision(ct_boot, CKKS.sk) / num_tests
                )
    print()
    print(f"PaCo precision: {float(round(precision_paco, 2))} bits.")
    print(f"Original precision: {float(round(precision_orig, 2))} bits.")
else:
    kappa = int(input("Number of simulated processors (for PaCo): kappa = "))
    log_kappa = log(kappa, 2)

    print()
    header = (
        f" {'C':^5} || "
        + " | ".join(
            [f"{'PaCo try ' + str(i+1):^12}" for i in range(num_tests)]
        )
        + " || "
        + " | ".join(
            [f"{'Orig try ' + str(i+1):^12}" for i in range(num_tests)]
        )
    )
    print(header)
    print("-" * len(header))

    for log_D in range(log_kappa + 1, min(log_B + log_kappa, log_N) + 1):
        D = 2**log_D
        C = D // kappa
        log_C = log(C, 2)
        print(f" {D:^5} |", end="", flush=True)

        # PaCo
        L = ceil((log_C + 1) / g) + ceil((log_C - 1) / g) + log_h + 3
        CKKS_x.config(N, L, q, p, delta)
        CKKS_x.key_gen(h)
        CKKS_x.config_PaCo(CKKS_x.sk, C, g, g)
        for _ in range(num_tests):
            coeffs = [
                randint(-p, p) if i % (N // D) == 0 else 0 for i in range(N)
            ]
            m = Poly(coeffs, N)
            ct = CKKS_x.enc_poly_with_sk(m, CKKS_x.sk) % q
            ct_boot, t = ct.parallel_PaCo(kappa, get_time=True)
            t = f"{float(round(t, 2))} s"
            print(f"| {t:^12} ", end="", flush=True)
        print("|", end="", flush=True)

        # Original
        L = 2 * ceil(max(log_D - 1, 1) / g) + log_d + r + 1
        CKKS.config(N, 2 ** (log_D - 1), L, q, p, delta)
        CKKS.key_gen(h)
        CKKS.config_bootstrap(CKKS.sk, d, r, g)
        for _ in range(num_tests):
            coeffs = [
                randint(-p, p) if i % (N // D) == 0 else 0 for i in range(N)
            ]
            m = Poly(coeffs, N)
            ct = CKKS.enc_poly_with_sk(m, CKKS.sk) % q
            t_start = time.time()
            ct_boot = ct.bootstrap()
            t_end = time.time()
            t = f"{float(round(t_end - t_start, 2))} s"
            print(f"| {t:^12} ", end="", flush=True)
        print()
