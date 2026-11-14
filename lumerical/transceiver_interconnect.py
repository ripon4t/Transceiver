"""
NRZ optical transceiver link in Ansys Lumerical INTERCONNECT (time-domain data-flow simulation),
built and run from Python with lumapi.

 CW laser -> Mach-Zehnder modulator <- NRZ pulse generator <- PRBS generator
 -> optical attenuator (channel) -> PIN photodetector -> Bessel low-pass filter
 -> oscilloscope + eye-diagram analyzer

The channel attenuation is swept; for every value the receiver waveform is exported to
../data/wave_att<dB>.npz so analysis/eye_analysis.py can compute eyes, Q and BER.

    python transceiver_interconnect.py [--hide] [--bitrate 10e9] [--atten 0 3 6 9 12 15]

Element and property names follow the INTERCONNECT element library. A property that does not
exist in your version is reported as "[warning] ..." and skipped; check it in the element's
Properties tab and adjust the CONFIG block below.
"""
import argparse
import os
import numpy as np
from lum_utils import import_lumapi, set_props

lumapi = import_lumapi()
C = 299792458.0
HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")

# --------------------------------------------------------------------------- CONFIG
WAVELENGTH = 1550e-9
LASER_POWER = 1e-3            # W
MZM_IL_DB, MZM_ER_DB = 6.0, 20.0
V_PI = 4.0                    # modulator V_pi (V)
NRZ_AMPLITUDE = 4.0           # drive swing (V) -> 0..V_pi gives a full on/off swing
RESPONSIVITY = 1.0            # A/W
THERMAL_NOISE = 20e-12        # A/sqrt(Hz)
BESSEL_CUTOFF_FACTOR = 0.75   # cutoff = factor * bitrate
SEQUENCE_BITS = 512
SAMPLES_PER_BIT = 64

# Element types (Element Library names) and the port names used to wire them
ELEMENTS = [
    # name,      type,                      x,   y,   properties
    ("LASER",   "CW Laser",                 100, 200, {"frequency": C / WAVELENGTH, "power": LASER_POWER}),
    ("PRBS",    "PRBS Generator",           100, 400, {}),
    ("NRZ",     "NRZ Pulse Generator",      300, 400, {"amplitude": NRZ_AMPLITUDE, "bias": 0.0}),
    ("MZM",     "Mach-Zehnder Modulator",   500, 250, {"modulator type": "unbalanced single drive",
                                                      "bias voltage 1": 0.0, "bias voltage 2": 0.0,
                                                      "pi dc voltage": V_PI, "pi rf voltage": V_PI,
                                                      "extinction ratio": MZM_ER_DB,
                                                      "insertion loss": MZM_IL_DB}),
    ("ATT",     "Optical Attenuator",       700, 250, {"attenuation": 0.0}),
    ("PIN",     "PIN Photodetector",        900, 250, {"responsivity": RESPONSIVITY,
                                                      "thermal noise": THERMAL_NOISE,
                                                      "enable thermal noise": True,
                                                      "enable shot noise": True}),
    ("FILT",    "Low Pass Bessel Filter",  1100, 250, {}),     # cutoff set from the bit rate below
    ("OSC_RX",  "Oscilloscope",            1300, 150, {}),
    ("EYE",     "Eye Diagram Analyzer",    1300, 350, {}),
    ("OSC_TX",  "Oscilloscope",             500, 520, {}),
]
CONNECTIONS = [
    # (from element, from port, to element, to port)
    ("LASER", "output", "MZM", "input"),
    ("PRBS", "output", "NRZ", "input"),
    ("NRZ", "output", "MZM", "modulation 1"),
    ("MZM", "output", "ATT", "input"),
    ("ATT", "output", "PIN", "input"),
    ("PIN", "output", "FILT", "input"),
    ("FILT", "output", "OSC_RX", "input"),
    ("FILT", "output", "EYE", "input"),
    ("NRZ", "output", "EYE", "reference"),      # clock recovery reference, as in the Ansys example
    ("NRZ", "output", "OSC_TX", "input"),
]


def build(ic, bitrate):
    ic.switchtodesign()
    # Global time-domain settings live on the root element
    for key, val in {"bitrate": bitrate, "simulation input": "sequence length",
                     "sequence length": SEQUENCE_BITS, "samples per bit": SAMPLES_PER_BIT}.items():
        try:
            ic.setnamed("::Root Element", key, val)
        except Exception as exc:
            print(f"[warning] root element: could not set '{key}': {exc}")

    for name, etype, x, y, props in ELEMENTS:
        ic.addelement(etype)
        props = dict(props)
        props.update({"name": name, "x position": x, "y position": y})
        set_props(ic, props, name)
    try:
        ic.setnamed("FILT", "cutoff frequency", BESSEL_CUTOFF_FACTOR * bitrate)
    except Exception as exc:
        print(f"[warning] FILT: could not set 'cutoff frequency': {exc}")

    for a, pa, b, pb in CONNECTIONS:
        try:
            ic.connect(a, pa, b, pb)
        except Exception as exc:
            print(f"[warning] could not connect {a}.{pa} -> {b}.{pb}: {exc}")


def dataset_to_waveform(res):
    """Turn a lumapi dataset (dict) into (time, signal) 1-D arrays."""
    if not isinstance(res, dict):
        raise TypeError(f"unexpected result type {type(res)}")
    keys = [k for k in res.keys() if not str(k).startswith("Lumerical")]
    tkey = next((k for k in keys if "time" in str(k).lower()), None)
    if tkey is None:
        raise KeyError(f"no time axis in result; available keys: {keys}")
    t = np.asarray(res[tkey], float).ravel()
    for k in keys:
        if k == tkey:
            continue
        v = np.asarray(res[k]).ravel()
        if v.size == t.size:
            return t, np.real(v).astype(float)
    raise KeyError(f"no signal array matching the time axis; available keys: {keys}")


def read_scope(ic, element, result_names=("signal", "waveform", "input")):
    last = None
    for rn in result_names:
        try:
            return dataset_to_waveform(ic.getresult(element, rn))
        except Exception as exc:
            last = exc
    raise RuntimeError(f"could not read a waveform from {element}: {last}. "
                       "Open the analyzer's Results tab and pass its result name via --result-name.")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hide", action="store_true")
    ap.add_argument("--bitrate", type=float, default=10e9)
    ap.add_argument("--atten", type=float, nargs="+", default=[0, 3, 6, 9, 12, 15],
                    help="channel attenuation values (dB)")
    ap.add_argument("--result-name", default=None, help="oscilloscope result name if not 'signal'")
    args = ap.parse_args()
    os.makedirs(DATA, exist_ok=True)

    names = (args.result_name,) if args.result_name else ("signal", "waveform", "input")
    with lumapi.INTERCONNECT(hide=args.hide) as ic:
        build(ic, args.bitrate)
        ic.save(os.path.join(DATA, "transceiver.icp"))
        for k, att in enumerate(args.atten):
            ic.switchtodesign()
            ic.setnamed("ATT", "attenuation", float(att))
            ic.run()
            t, v = read_scope(ic, "OSC_RX", names)
            np.savez(os.path.join(DATA, f"wave_att{att:05.1f}dB.npz"), t=t, v=v, att_db=att, bitrate=args.bitrate)
            print(f"attenuation {att:5.1f} dB: {t.size} samples saved")
            if k == 0:
                try:
                    tt, vt = read_scope(ic, "OSC_TX", names)
                    np.savez(os.path.join(DATA, "tx_drive.npz"), t=tt, v=vt)
                except Exception as exc:
                    print(f"[warning] could not export the transmitter drive waveform: {exc}")


if __name__ == "__main__":
    main()
