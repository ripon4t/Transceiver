"""
Stand-alone Python model of the NRZ transceiver link (no Lumerical needed):

 PRBS-9 -> NRZ drive (finite rise time) -> Mach-Zehnder modulator (analytical phase model,
 finite extinction ratio, insertion loss) -> attenuator -> PIN photodiode (shot + thermal noise)
 -> 4th-order Bessel low-pass filter

It writes waveforms in the same format as the INTERCONNECT export and runs the same analysis,
so it can be used to check the analysis code and to see what the Lumerical result should look like.

    python link_reference.py
"""
import os
import sys
import numpy as np
from scipy import signal
from scipy.ndimage import gaussian_filter1d

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "analysis"))
import eye_analysis  # noqa: E402

Q_E = 1.602176634e-19
# ---- link parameters (illustrative) ----
BITRATE = 10e9
SPS = 64                       # samples per bit
N_BITS = 512
P_LASER = 1e-3                 # W
IL_DB, ER_DB = 6.0, 20.0       # modulator insertion loss / extinction ratio
V_PI, V_DRIVE = 4.0, 4.0       # modulator V_pi and NRZ drive swing (V)
RISE_FRACTION = 0.3            # 10-90 rise time as a fraction of a bit period
RESPONSIVITY = 1.0             # A/W
THERMAL_NOISE = 20e-12         # A/sqrt(Hz)
BESSEL_BW = 0.75 * BITRATE
ATTENUATIONS_DB = [0, 3, 6, 9, 12, 15]


def prbs9(n):
    state = [1] * 9
    out = []
    for _ in range(n):
        new = state[8] ^ state[4]
        out.append(state[8])
        state = [new] + state[:-1]
    return np.array(out, dtype=float)


def simulate(att_db, seed):
    rng = np.random.default_rng(seed)
    fs = BITRATE * SPS
    bits = prbs9(N_BITS)
    x = np.repeat(bits, SPS)
    sigma = RISE_FRACTION * SPS / 2.563
    x = gaussian_filter1d(x, sigma, mode="wrap")                  # finite rise/fall time
    # Mach-Zehnder: power transmission from the phase difference between the arms
    delta = 10 ** (-ER_DB / 20)
    dphi = np.pi * (V_DRIVE * x) / V_PI
    trans = ((1 + delta**2) / 2 + (1 - delta**2) / 2 * np.cos(dphi)) * 10 ** (-IL_DB / 10)
    p_opt = P_LASER * trans * 10 ** (-att_db / 10)
    i = RESPONSIVITY * p_opt
    sig_shot = np.sqrt(2 * Q_E * np.clip(i, 0, None) * fs / 2)
    sig_th = THERMAL_NOISE * np.sqrt(fs / 2)
    i_noisy = i + rng.normal(0, 1, i.size) * np.sqrt(sig_shot**2 + sig_th**2)
    sos = signal.bessel(4, BESSEL_BW, btype="low", fs=fs, output="sos", norm="mag")
    v = signal.sosfilt(sos, i_noisy)
    t = np.arange(v.size) / fs
    return t, v


def main():
    data = os.path.join(HERE, "data")
    os.makedirs(data, exist_ok=True)
    for k, a in enumerate(ATTENUATIONS_DB):
        t, v = simulate(a, seed=100 + k)
        np.savez(os.path.join(data, f"wave_att{a:05.1f}dB.npz"), t=t, v=v, att_db=a, bitrate=BITRATE)
    rows = eye_analysis.analyze(data, os.path.join(HERE, "images"), label="(Python reference link)")
    print(f"{'att(dB)':>8} {'Q':>8} {'BER':>10}")
    for r in rows:
        print(f"{r['att_db']:8.1f} {r['Q']:8.2f} {r['BER']:10.2e}")


if __name__ == "__main__":
    main()
