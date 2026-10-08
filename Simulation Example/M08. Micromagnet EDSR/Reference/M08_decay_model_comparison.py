__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

# [VALIDATION] M08 2차 라운드 Step 2 — Decay-law 모델 비교
#
# 분류: B. Parameter sensitivity / numerical validation study.
#   EDSR_noise.py, MOS_EDSR.py는 수정하지 않는다.
#
# [BASELINE] EDSR_noise.py와 동일한 noise 모델(Lorentzian, S0=0.01, num_runs=100)로
#   20/50/100/200 Rabi cycle에 대해 raw ensemble-averaged sigma_z(t)를 다시 계산하고
#   (이번엔 raw 배열을 저장), 세 모델(Gaussian/Exponential/Stretched-exponential)을
#   각각 fit해서 비교한다.
#
# [WARNING] 특정 모델이 "더 그럴듯해 보인다"고 임의로 선택하지 않는다 — RMSE/AIC/BIC로
#   정량 비교하고, fit-window dependence로 단일 decay constant가 잘 정의되는지 확인한다.

import json
import pathlib
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import qutip
from qtcad.device import constants as ct
from qtcad.qubit import spectra, dynamics, noise

script_dir = pathlib.Path(__file__).parent.resolve()
output_dir = script_dir / "output"
export_dir = output_dir / "analysis_exports"
export_dir.mkdir(parents=True, exist_ok=True)

# --- [BASELINE] h0/u 로드 (Step 1 파이프라인 재사용) ---
h0_u_file = output_dir / "mos_edsr_h0_u.npz"
_loaded = np.load(h0_u_file)
h0 = _loaded["h0"]
u = _loaded["u"]
omega_rabi_loaded = float(_loaded["omega_rabi"])

H0 = h0 / ct.hbar
delta_V = u / ct.hbar
omega = np.absolute(H0[0, 0] - H0[1, 1])
omega_Rabi = np.absolute(delta_V[0, 1])
_rel_err = abs(omega_Rabi - omega_rabi_loaded) / abs(omega_rabi_loaded)
if _rel_err > 1e-3:
    raise RuntimeError(f"h0/u consistency check FAILED (rel.err={_rel_err:.2e})")
print(f"[consistency check] PASSED (rel.err={_rel_err:.2e})")
T_Rabi = 2 * np.pi / omega_Rabi
print(f"T_Rabi = {T_Rabi*1e9:.2f} ns, f_Rabi = {omega_Rabi/(2*np.pi)*1e-6:.4f} MHz\n")

# --- [BASELINE] Noise 모델, EDSR_noise.py와 동일 ---
S0 = 0.01
wc = omega_Rabi
w0 = 0
POINTS_PER_CYCLE = 50
NUM_RUNS = 100

dyn = dynamics.Dynamics()
N = noise.Noise()

cycles_list = [20, 50, 100, 200]
raw_data = {}

for n_cycles in cycles_list:
    t0 = time.time()
    T = n_cycles * T_Rabi
    npts = int(POINTS_PER_CYCLE * n_cycles)
    times = np.linspace(0, T, npts)
    T0 = max(50 * T_Rabi, 2 * T)
    omega_max = 5 * wc
    wk = np.arange(0, omega_max, 2 * np.pi / T0) + w0
    spectrum = spectra.lorentz(wk, S0, wc, w0=w0)
    psi0 = qutip.basis(2, 0)
    sigz = N.dynamics(H0, delta_V, omega, spectrum, T0, psi0, times, qutip.sigmaz(),
                       vec_omega=wk, num_runs=NUM_RUNS)
    elapsed = time.time() - t0
    raw_data[n_cycles] = (times, sigz)
    np.savez(export_dir / f"m08_raw_ensemble_{n_cycles}cycles.npz", times=times, sigz=sigz)
    print(f"[n_cycles={n_cycles}] ensemble computed in {elapsed:.1f}s, npts={npts}")

# --- [VALIDATION] 3개 decay 모델 정의 ---
def model_gaussian(t, A, T, phi, C):
    return A * np.cos(omega_Rabi * t + phi) * np.exp(-(t / T) ** 2) + C

def model_exponential(t, A, T, phi, C):
    return A * np.cos(omega_Rabi * t + phi) * np.exp(-t / T) + C

def model_stretched(t, A, T, phi, C, p):
    return A * np.cos(omega_Rabi * t + phi) * np.exp(-(t / T) ** p) + C

MODELS = {
    "Gaussian (p=2 fixed)": (model_gaussian, [1.0, None, 0.0, 0.0], 4),
    "Exponential (p=1 fixed)": (model_exponential, [1.0, None, 0.0, 0.0], 4),
    "Stretched exponential": (model_stretched, [1.0, None, 0.0, 0.0, 1.0], 5),
}


def fit_and_score(t, y, model_func, p0_template, k_params, T_guess):
    p0 = [x if x is not None else T_guess for x in p0_template]
    try:
        popt, pcov = curve_fit(model_func, t, y, p0=p0, maxfev=50000)
        perr = np.sqrt(np.diag(pcov))
        y_fit = model_func(t, *popt)
        resid = y - y_fit
        rmse = np.sqrt(np.mean(resid ** 2))
        n = len(t)
        rss = np.sum(resid ** 2)
        # Gaussian-noise likelihood assumption for AIC/BIC
        aic = n * np.log(rss / n + 1e-300) + 2 * k_params
        bic = n * np.log(rss / n + 1e-300) + k_params * np.log(n)
        T_idx = 1  # T is always popt[1] by construction above
        return {
            "T_us": popt[T_idx] * 1e6, "T_err_us": perr[T_idx] * 1e6,
            "p_fit": popt[4] if k_params == 5 else np.nan,
            "p_err": perr[4] if k_params == 5 else np.nan,
            "RMSE": rmse, "AIC": aic, "BIC": bic, "fit_valid": True,
            "popt": popt, "resid": resid,
        }
    except Exception as exc:
        return {"T_us": np.nan, "T_err_us": np.nan, "p_fit": np.nan, "p_err": np.nan,
                "RMSE": np.nan, "AIC": np.nan, "BIC": np.nan, "fit_valid": False,
                "popt": None, "resid": None, "error": str(exc)}


results_rows = []
fit_cache = {}
for n_cycles in [50, 100, 200]:  # 20 cycles는 diagnostic만 (지침 6절)
    t_arr, y_arr = raw_data[n_cycles]
    T_guess = n_cycles * T_Rabi * 2
    for model_name, (func, p0_t, k) in MODELS.items():
        res = fit_and_score(t_arr, y_arr, func, p0_t, k, T_guess)
        res["cycles"] = n_cycles
        res["model"] = model_name
        fit_cache[(n_cycles, model_name)] = res
        results_rows.append({
            "cycles": n_cycles, "model": model_name, "T_us": res["T_us"],
            "T_err_us": res["T_err_us"], "p_fit": res["p_fit"], "p_err": res["p_err"],
            "RMSE": res["RMSE"], "AIC": res["AIC"], "BIC": res["BIC"],
            "fit_valid": res["fit_valid"],
        })
        print(f"[{n_cycles} cycles, {model_name}] T={res['T_us']:.2f}+-{res['T_err_us']:.2f}us, "
              f"p={res['p_fit']}, RMSE={res['RMSE']:.4f}, AIC={res['AIC']:.1f}, BIC={res['BIC']:.1f}")

results_df = pd.DataFrame(results_rows)
print("\n" + "=" * 100)
print("MODEL COMPARISON TABLE")
print("=" * 100)
print(results_df.to_string(index=False))
results_df.to_csv(export_dir / "m08_decay_fit_results.csv", index=False)

# --- Plot: fit overlay + residuals, for 200 cycles (가장 데이터가 많은 case) ---
t200, y200 = raw_data[200]
fig, axs = plt.subplots(2, 1, figsize=(10, 7), sharex=True, gridspec_kw={"height_ratios": [2, 1]})
axs[0].plot(t200 * 1e6, y200, ".", ms=2, alpha=0.4, color="0.5", label="data (200 cycles)")
colors = {"Gaussian (p=2 fixed)": "tab:blue", "Exponential (p=1 fixed)": "tab:green",
          "Stretched exponential": "tab:red"}
for model_name, (func, _, k) in MODELS.items():
    res = fit_cache[(200, model_name)]
    if res["fit_valid"]:
        y_fit = func(t200, *res["popt"])
        axs[0].plot(t200 * 1e6, y_fit, "-", lw=1.2, color=colors[model_name],
                    label=f"{model_name}: T={res['T_us']:.1f}us")
        axs[1].plot(t200 * 1e6, res["resid"], ".", ms=2, alpha=0.5, color=colors[model_name],
                    label=model_name)
axs[0].set_ylabel(r"$\langle\sigma_z(t)\rangle$"); axs[0].legend(fontsize=8)
axs[1].set_xlabel("t (us)"); axs[1].set_ylabel("residual"); axs[1].legend(fontsize=7)
axs[1].axhline(0, color="k", lw=0.5)
fig.suptitle("M08 decay-model comparison (200 Rabi cycles, 100-run ensemble)")
fig.tight_layout()
fig.savefig(export_dir / "m08_decay_model_fits_200cycles.png", dpi=150)
plt.close(fig)

# --- T_decay vs simulation window (각 model) ---
fig, axs = plt.subplots(1, 2, figsize=(12, 5))
for model_name in MODELS:
    cyc, Ts, perrs = [], [], []
    for n_cycles in [50, 100, 200]:
        res = fit_cache[(n_cycles, model_name)]
        cyc.append(n_cycles); Ts.append(res["T_us"]); perrs.append(res["T_err_us"])
    axs[0].errorbar(cyc, Ts, yerr=perrs, marker="o", label=model_name, color=colors[model_name])
axs[0].set_xlabel("simulation window (Rabi cycles)"); axs[0].set_ylabel("fitted T (us)")
axs[0].set_title("T_decay vs. simulation window"); axs[0].legend(fontsize=8); axs[0].grid(True)

cyc, ps, perrs = [], [], []
for n_cycles in [50, 100, 200]:
    res = fit_cache[(n_cycles, "Stretched exponential")]
    cyc.append(n_cycles); ps.append(res["p_fit"]); perrs.append(res["p_err"])
axs[1].errorbar(cyc, ps, yerr=perrs, marker="o", color="tab:red")
axs[1].axhline(1, color="tab:green", ls="--", lw=1, label="p=1 (exponential)")
axs[1].axhline(2, color="tab:blue", ls="--", lw=1, label="p=2 (Gaussian)")
axs[1].set_xlabel("simulation window (Rabi cycles)"); axs[1].set_ylabel("stretched exponent p")
axs[1].set_title("Stretched exponent vs. simulation window"); axs[1].legend(fontsize=8); axs[1].grid(True)
fig.tight_layout()
fig.savefig(export_dir / "m08_T_and_p_vs_window.png", dpi=150)
plt.close(fig)

# --- Fit-window dependence test (200-cycle data, 0-50%/0-75%/0-100%) ---
print("\n" + "=" * 100)
print("FIT-WINDOW DEPENDENCE TEST (200 cycles)")
print("=" * 100)
window_rows = []
for frac in [0.5, 0.75, 1.0]:
    n_end = int(len(t200) * frac)
    t_sub, y_sub = t200[:n_end], y200[:n_end]
    for model_name, (func, p0_t, k) in MODELS.items():
        T_guess = t_sub[-1] * 2
        res = fit_and_score(t_sub, y_sub, func, p0_t, k, T_guess)
        window_rows.append({"fit_window_frac": frac, "model": model_name,
                            "T_us": res["T_us"], "T_err_us": res["T_err_us"],
                            "p_fit": res["p_fit"], "RMSE": res["RMSE"]})
        print(f"[window={frac*100:.0f}%, {model_name}] T={res['T_us']:.2f}us")

window_df = pd.DataFrame(window_rows)
window_df.to_csv(export_dir / "m08_fit_window_dependence.csv", index=False)

# 판정: T가 fit window에 얼마나 민감한가
for model_name in MODELS:
    sub = window_df[window_df["model"] == model_name]
    T_vals = sub["T_us"].values
    if np.all(~np.isnan(T_vals)) and T_vals[0] != 0:
        rel_range = (np.max(T_vals) - np.min(T_vals)) / np.min(T_vals) * 100
        print(f"\n[{model_name}] fit-window에 따른 T 상대 변화 범위: {rel_range:.1f}%")

print(f"\nExports: {export_dir}")
