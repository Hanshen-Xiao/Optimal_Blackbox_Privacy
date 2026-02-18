import numpy as np
from math import sqrt, pi, log, exp
import hashlib
import gc
import time as _time
from math import gcd
from Crypto.Util.number import getPrime

class DeterministicRand:
    def __init__(self, seed: bytes):
        self.seed = seed
        self.counter = 0
        self.buf = b""

    def __call__(self, n: int) -> bytes:
        out = bytearray()
        while len(out) < n:
            if not self.buf:
                self.counter += 1
                self.buf = hashlib.sha256(
                    self.seed + self.counter.to_bytes(8, "big")
                ).digest()
            need = n - len(out)
            take = self.buf[:need]
            out.extend(take)
            self.buf = self.buf[len(take):]
        return bytes(out)


def generate_rsa_key_from_seed(key_bits=512, seed=0, e=65537):
    seed_bytes = f"rsa-seed:{int(seed)}".encode("utf-8")
    drand = DeterministicRand(seed_bytes)

    while True:
        p = getPrime(key_bits // 2, randfunc=drand)
        q = getPrime(key_bits // 2, randfunc=drand)
        if p == q:
            continue

        phi = (p - 1) * (q - 1)
        if gcd(e, phi) != 1:
            continue

        N = p * q
        d = pow(int(e), -1, int(phi))
        return N, int(e), int(d), p, q
    

def RSA_total_time(
    seed,
    v,
    num_encrypt=1,
    key_bits=512,
    e=65537,
    trials=1,
    trim=0.0,
    disable_gc=True,
):
    if disable_gc:
        gc_was_enabled = gc.isenabled()
        gc.disable()
    else:
        gc_was_enabled = None

    times = np.empty(int(trials), dtype=np.int64)

    try:
        for t in range(int(trials)):
            t1 = _time.perf_counter_ns()

            N, e2, d, p, q = generate_rsa_key_from_seed(
                key_bits=int(key_bits),
                seed=seed,   
                e=int(e),
            )

            N = int(N)
            e2 = int(e2)

            v_mod = int(v) % N
            c_last = 0
            for _ in range(int(num_encrypt)):
                c_last = pow(v_mod, e2, N)

            t2 = _time.perf_counter_ns()
            times[t] = (t2 - t1)*1e-5

        times.sort()
        k = int(trim * trials)
        if 2 * k < trials:
            times = times[k: trials - k]

        return int(np.median(times))

    finally:
        if disable_gc and gc_was_enabled:
            gc.enable()

def hw_leakage(x, v):
    return RSA_total_time(x, v)

def gaussian_pdf(o, mean, sigma2):
    return (1.0 / sqrt(2 * pi * sigma2)) * np.exp(-(o - mean)**2 / (2 * sigma2))

def alpha_divergence(p, q):
    return np.sum((p**alpha) * (q**(1-alpha)))

alpha = 10
T = 5
m = 100
num_bits = 8
keyspace = np.arange(2**num_bits)
prior = np.ones_like(keyspace, dtype=float) / len(keyspace)

r_t = [0.1] * T
r_t_prime = [0.1] * T
gamma_t = [0.0001] * T
R_t = [(r_t[i] - (1 - gamma_t[i]) * r_t_prime[i]) / gamma_t[i] for i in range(T)]
T_max = 8

def A_t_minus_1(x, transcript):

    A = prior[x]**(alpha - 1)

    for info in transcript:

        v_s = info["v"]
        o_s = info["o"]
        sigma2_s = info["sigma2"]
        output_grid = info["output_grid"]
        W_density = info["W_density"]

        mean = hw_leakage(x, v_s)

        pM = gaussian_pdf(o_s, mean, sigma2_s)
        pM = max(pM, 1e-300)

        pW = np.interp(o_s, output_grid, W_density)
        pW = max(pW, 1e-300)

        A *= (pM / pW)**alpha

    return A

secret = np.random.choice(keyspace, p=prior)
print("True secret:", secret)

transcript = []

v_t = np.random.randint(0, 256)

for t in range(T):

    print("\nRound", t+1)

    X_sample = np.random.choice(keyspace, size=m, p=prior)

    def compute_W(sigma2):

        output_grid = np.linspace(-2, 12, 2000)
        mixture = []

        for o in output_grid:
            vals = []
            for x in X_sample:
                A_val = A_t_minus_1(x, transcript)
                mean = hw_leakage(x, v_t)
                p = gaussian_pdf(o, mean, sigma2)
                vals.append(A_val * (p**alpha))

            mixture.append((np.mean(vals))**(1/alpha))

        mixture = np.array(mixture)
        C = np.trapz(mixture, output_grid)
        C = max(C, 1e-300)

        W_density = mixture / C
        return output_grid, W_density, C


    sigma2_low = 1e-6
    sigma2_high = 10

    X_val = np.random.choice(keyspace, size=m, p=prior)

    def check_feasible(sigma2):

        output_grid, W_density, C_t = compute_W(sigma2)

        D_vals = []
        A_vals = []

        for x in X_val:
            A_val = A_t_minus_1(x, transcript)
            A_vals.append(A_val)

            mean = hw_leakage(x, v_t)
            p = gaussian_pdf(output_grid, mean, sigma2)
            D = alpha_divergence(p, W_density)
            D_vals.append(D)

        L = np.mean([A_vals[i] * D_vals[i] for i in range(m)])
        K = np.mean(A_vals)

        # Compute Bound_D
        log_Bound_D = (alpha - 1) * np.log(C_t) + \
                      alpha*(alpha-1)/(2*sigma2)*T_max**2

        if log_Bound_D > 50:
            Bound_D = np.inf
        else:
            Bound_D = np.exp(log_Bound_D)

        # Compute Bound_A
        max_prior = np.max(prior)
        Bound_A = max_prior**(alpha - 1)
        
        for info in transcript:

            o_s = info["o"]
            sigma2_s = info["sigma2"]
            output_grid_s = info["output_grid"]
            W_density_s = info["W_density"]

            # compute p_W_s(o_s)
            pW_s = np.interp(o_s, output_grid_s, W_density_s)
            pW_s = max(pW_s, 1e-300)

            # define d(o_s)
            if 0 < o_s < T_max:
                d_os = 0
            elif o_s >= T_max:
                d_os = o_s - T_max
            else:
                d_os = -o_s

            Bound_A *= np.exp(
                alpha * (
                    - d_os**2 / (2 * sigma2_s)
                    - np.log(pW_s)
                )
            )

        beta = Bound_A * Bound_D * sqrt(log(1/gamma_t[t])/(2*m))

        feasible = (L <= r_t_prime[t]*K - beta) and (Bound_D <= R_t[t])

        return feasible


    while True:
        if check_feasible(sigma2_high):
            break
        sigma2_high *= 2
        if sigma2_high > 1e4:
            break

    for _ in range(40):
        sigma2_mid = (sigma2_low + sigma2_high) / 2
        if check_feasible(sigma2_mid):
            sigma2_high = sigma2_mid
        else:
            sigma2_low = sigma2_mid

    sigma2 = sigma2_high
    print("Calibrated sigma^2:", sigma2)

    noise = np.random.normal(0, sqrt(sigma2))
    o_t = hw_leakage(secret, v_t) + noise

    output_grid, W_density, _ = compute_W(sigma2)

    transcript.append({
        "v": v_t,
        "o": o_t,
        "sigma2": sigma2,
        "output_grid": output_grid,
        "W_density": W_density
    })

    print("Observed o_t:", o_t)

print("\nTranscript:")
for r in transcript:
    print((r["v"], r["o"]))
