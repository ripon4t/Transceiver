"""
Eye-diagram and signal-quality analysis of receiver waveforms.

Input  : data/wave_att*.npz   (one file per channel attenuation), each containing
           t        time vector (s)
           v        receiver output after the Bessel filter (A)
           att_db   channel attenuation (dB)
           bitrate  bit rate (bit/s)
Output : figures in images/ and a metrics table in data/metrics.json

The same code is used for the waveforms exported from Lumerical INTERCONNECT and for the
Python reference link, so the two can be compared directly.

    python eye_analysis.py [--data ../data] [--out ../images]
"""
import argparse
import glob
import json
import os
import numpy as np
from scipy.special import erfc
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
BLUE, ORANGE, RED = "#2a6fbb", "#e07a1f", "#c0392b"
plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3})
SKIP_BITS = 8          # discard start-up transient


# ----------------------------------------------------------------------------- loading
def load_runs(data_dir):
    runs = []
    for f in sorted(glob.glob(os.path.join(data_dir, "wave_att*.npz"))):
        d = np.load(f)
        t, v = np.asarray(d["t"], float).ravel(), np.asarray(d["v"], float).ravel()
        dt = np.median(np.diff(t))
        tu = np.arange(t[0], t[-1], dt)                      # make the time base uniform
        runs.append(dict(att=float(d["att_db"]), t=tu, v=np.interp(tu, t, v), br=float(d["bitrate"])))
    if not runs:
        raise FileNotFoundError(f"No wave_att*.npz files found in {data_dir}")
    runs.sort(key=lambda r: r["att"])
    return runs


# ----------------------------------------------------------------------------- sampling
def sample(run, phase):
    T = 1.0 / run["br"]
    n_bits = int((run["t"][-1] - run["t"][0]) / T)
    tau = run["t"][0] + (np.arange(SKIP_BITS, n_bits - 2) + phase) * T
    return np.interp(tau, run["t"], run["v"])


def _q_of(samples, thr):
    hi, lo = samples[samples > thr], samples[samples <= thr]
    if hi.size < 5 or lo.size < 5:
        return -np.inf
    return (hi.mean() - lo.mean()) / (hi.std() + lo.std() + 1e-30)


def best_phase(run, n_phase=128):
    """Sampling phase (fraction of a bit) that opens the eye the most."""
    best, best_q = 0.0, -np.inf
    for ph in np.linspace(0, 1, n_phase, endpoint=False):
        s = sample(run, ph)
        thr = 0.5 * (np.percentile(s, 98) + np.percentile(s, 2))
        q = _q_of(s, thr)
        if q > best_q:
            best, best_q = ph, q
    return best


def metrics(run, phase, ref_bits):
    s = sample(run, phase)
    n = min(s.size, ref_bits.size)
    s, b = s[:n], ref_bits[:n]
    hi, lo = s[b], s[~b]
    q = (hi.mean() - lo.mean()) / (hi.std() + lo.std())
    return dict(att_db=run["att"], mean_1=float(hi.mean()), mean_0=float(lo.mean()),
                sigma_1=float(hi.std()), sigma_0=float(lo.std()), Q=float(q),
                BER=float(0.5 * erfc(q / np.sqrt(2))),
                eye_height=float((hi.mean() - 3 * hi.std()) - (lo.mean() + 3 * lo.std())))


# ----------------------------------------------------------------------------- driver
def analyze(data_dir, out_dir, label=""):
    os.makedirs(out_dir, exist_ok=True)
    runs = load_runs(data_dir)
    clean = runs[0]                                   # lowest attenuation = reference for the bit pattern
    phase = best_phase(clean)
    s_clean = sample(clean, phase)
    thr = 0.5 * (np.percentile(s_clean, 98) + np.percentile(s_clean, 2))
    ref_bits = s_clean > thr
    rows = [metrics(r, phase, ref_bits) for r in runs]
    T = 1.0 / clean["br"]

    # ---- Fig 1: eye diagrams ----
    pick = runs if len(runs) <= 6 else [runs[i] for i in np.linspace(0, len(runs) - 1, 6).astype(int)]
    ncol = 3 if len(pick) > 3 else len(pick)
    nrow = int(np.ceil(len(pick) / ncol))
    fig, axs = plt.subplots(nrow, ncol, figsize=(4.6 * ncol, 3.6 * nrow), squeeze=False)
    for ax, r in zip(axs.ravel(), pick):
        row = next(rw for rw in rows if rw["att_db"] == r["att"])
        lo_lvl, span = row["mean_0"], row["mean_1"] - row["mean_0"]
        m = r["t"] >= r["t"][0] + SKIP_BITS * T
        tau = r["t"][m] - r["t"][0] - phase * T
        x = ((tau + T) % (2 * T)) / T - 1.0
        ax.set_facecolor("black")
        ax.hist2d(x, (r["v"][m] - lo_lvl) / span, bins=(120, 140), range=[[-1, 1], [-0.5, 1.5]],
                  cmap="inferno", cmin=1)
        q = next(row["Q"] for row in rows if row["att_db"] == r["att"])
        ax.set_title(f"{r['att']:.0f} dB attenuation   Q = {q:.1f}", fontsize=10)
        ax.set_xlabel("Time (bit periods)"); ax.set_ylabel("Signal (norm. to own eye amplitude)"); ax.grid(False)
    for ax in axs.ravel()[len(pick):]:
        ax.axis("off")
    fig.suptitle(f"Eye diagrams vs channel attenuation {label}".strip())
    fig.tight_layout(); fig.savefig(os.path.join(out_dir, "01_eye_diagrams.png"), dpi=200); plt.close(fig)

    # ---- Fig 2: Q and BER ----
    att = np.array([r["att_db"] for r in rows]); q = np.array([r["Q"] for r in rows])
    ber = np.array([max(r["BER"], 1e-30) for r in rows])
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.2))
    axs[0].plot(att, q, "o-", color=BLUE); axs[0].axhline(6, ls="--", color="gray", lw=1)
    axs[0].text(att.min(), 6.3, "Q = 6 (BER ≈ 1e-9)", fontsize=9, color="gray")
    axs[0].set_xlabel("Channel attenuation (dB)"); axs[0].set_ylabel("Q factor"); axs[0].set_title("Signal quality")
    axs[1].semilogy(att, ber, "o-", color=RED); axs[1].axhline(1e-9, ls="--", color="gray", lw=1)
    axs[1].set_xlabel("Channel attenuation (dB)"); axs[1].set_ylabel("Estimated BER (Gaussian)"); axs[1].set_title("Bit error rate")
    fig.tight_layout(); fig.savefig(os.path.join(out_dir, "02_q_and_ber_vs_attenuation.png"), dpi=200); plt.close(fig)

    # ---- Fig 3: time waveforms ----
    fig, axs = plt.subplots(2, 1, figsize=(9, 5.6), sharex=True)
    for ax, r, col in zip(axs, (runs[0], runs[-1]), (BLUE, ORANGE)):
        t0 = r["t"][0] + SKIP_BITS * T
        m = (r["t"] >= t0) & (r["t"] < t0 + 24 * T)
        ax.plot((r["t"][m] - t0) * 1e9, r["v"][m] * 1e6, color=col, lw=1)
        ax.set_ylabel("Receiver output (µA)"); ax.set_title(f"{r['att']:.0f} dB attenuation", fontsize=10)
    axs[1].set_xlabel("Time (ns)")
    fig.tight_layout(); fig.savefig(os.path.join(out_dir, "03_time_waveforms.png"), dpi=200); plt.close(fig)

    with open(os.path.join(data_dir, "metrics.json"), "w") as f:
        json.dump(dict(sampling_phase=phase, bitrate=clean["br"], rows=rows), f, indent=2)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default=os.path.join(HERE, "..", "data"))
    ap.add_argument("--out", default=os.path.join(HERE, "..", "images"))
    args = ap.parse_args()
    rows = analyze(args.data, args.out)
    print(f"{'att(dB)':>8} {'Q':>8} {'BER':>10} {'eye height (uA)':>16}")
    for r in rows:
        print(f"{r['att_db']:8.1f} {r['Q']:8.2f} {r['BER']:10.2e} {r['eye_height']*1e6:16.2f}")


if __name__ == "__main__":
    main()
