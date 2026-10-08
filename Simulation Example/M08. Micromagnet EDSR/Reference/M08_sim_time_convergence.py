__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

# [VALIDATION] M08 Step 2 — 82us "T2*" fit 수렴성 검증 (Physics-oriented validation study)
#
# 분류: B. Parameter sensitivity / numerical convergence study.
#   Official tutorial(EDSR_noise.py)도, device prediction도 아니다.
#   목적: 현재 분석 notebook(M08_Micromagnet_EDSR_Analysis.ipynb)이 보고하는
#   damped-cosine fit 값(~82 us)이 시뮬레이션 창(20 Rabi cycles ≈ 13.2 us)에 비해
#   6배 이상 긴 외삽인지, 시뮬레이션 시간을 늘리면 fit이 수렴하는지 확인한다.
#
# [BASELINE] 원본 EDSR_noise.py/MOS_EDSR.py, 분석 notebook은 전혀 수정하지 않는다.
#   이 스크립트는 EDSR_noise.py와 동일한 noise 모델(Lorentzian, S0=0.01)과 h0/u
#   파이프라인(Step 1에서 수정된 mos_edsr_h0_u.npz 로드)을 재사용하되, 20/50/100/200
#   Rabi cycle로 시뮬레이션 창만 늘려 비교한다.
#
# [NUMERICAL] fit model은 notebook의 damped_rabi와 동일: cos(omega_Rabi*t+phase) *
#   exp(-(t/Tdecay)^2) (진폭/오프셋은 1/0으로 고정, notebook과 동일 조건으로 비교
#   하기 위함).
#
# Fit-validity gate (project-level practical criterion, QTCAD 공식 기준 아님):
#   시뮬레이션 창 끝에서 fitted envelope가 0.7 미만으로 떨어져야 "decay resolved"로
#   인정한다. 그렇지 않으면 "decay unresolved within simulation window"로 보고하고
#   숫자를 공식 결과로 제시하지 않는다.

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

# --- [VALIDATION] h0/u 로드 (Step 1에서 만든 파이프라인, EDSR_noise.py와 동일 로직) ---
h0_u_file = output_dir / "mos_edsr_h0_u.npz"
if not h0_u_file.exists():
    raise FileNotFoundError(
        f"{h0_u_file} 가 없습니다. 먼저 MOS_EDSR.py를 실행해야 합니다."
    )
_loaded = np.load(h0_u_file)
h0 = _loaded["h0"]
u = _loaded["u"]
omega0_loaded = float(_loaded["omega0"])
omega_rabi_loaded = float(_loaded["omega_rabi"])

H0 = h0 / ct.hbar  # [NUMERICAL] rad/s 단위, plain numpy array (Qobj 아님 — EDSR_noise.py와 동일 관례)
delta_V = u / ct.hbar

omega = np.absolute(H0[0, 0] - H0[1, 1])
omega_Rabi = np.absolute(delta_V[0, 1])
_rel_err = abs(omega_Rabi - omega_rabi_loaded) / abs(omega_rabi_loaded)
if _rel_err > 1e-3:
    raise RuntimeError(f"h0/u consistency check FAILED (rel.err={_rel_err:.2e})")
print(f"[consistency check] PASSED (rel.err={_rel_err:.2e})\n")

T_Rabi = 2 * np.pi / omega_Rabi
print(f"T_Rabi = {T_Rabi*1e9:.2f} ns, omega_Rabi = {omega_Rabi:.6e} rad/s "
      f"(f_Rabi = {omega_Rabi/(2*np.pi)*1e-6:.4f} MHz)\n")

# --- [BASELINE] Noise 모델, EDSR_noise.py와 동일 ---
S0 = 0.01  # Total power
wc = omega_Rabi  # Cutoff frequency
w0 = 0  # Central frequency
POINTS_PER_CYCLE = 50  # baseline: 1000 points / 20 cycles = 50 points/cycle, 동일 밀도 유지
NUM_RUNS = 100  # baseline과 동일

dyn = dynamics.Dynamics()
N = noise.Noise()


def damped_rabi(t, t2star, phase):
    """notebook의 fit 모델과 동일 (진폭=1, 오프셋=0 고정)."""
    return np.cos(omega_Rabi * t + phase) * np.exp(-(t / t2star) ** 2)


# [NUMERICAL] Computational Cost Guard (지침 42절) — 각 case 실행 전 예상 비용을 출력.
cycles_list = [20, 50, 100, 200]
print("=" * 80)
print("COMPUTATIONAL COST ESTIMATE")
print("=" * 80)
for n_cycles in cycles_list:
    npts = int(POINTS_PER_CYCLE * n_cycles)
    print(f"  n_cycles={n_cycles:4d}: T_sim={n_cycles*T_Rabi*1e6:7.2f} us, "
          f"{npts:6d} time points x {NUM_RUNS} noise realizations")
print("=" * 80 + "\n")

results = []
fig, ax = plt.subplots(figsize=(9, 6))

for n_cycles in cycles_list:
    t_start = time.time()
    T = n_cycles * T_Rabi
    npts = int(POINTS_PER_CYCLE * n_cycles)
    times = np.linspace(0, T, npts)

    # [NUMERICAL] T0 must satisfy 2*pi/T0 << wc to resolve the spectrum, and
    # T0 should be >> largest timescale of the problem (here, T itself).
    T0 = max(50 * T_Rabi, 2 * T)
    omega_max = 5 * wc
    wk = np.arange(0, omega_max, 2 * np.pi / T0) + w0
    spectrum = spectra.lorentz(wk, S0, wc, w0=w0)

    psi0 = qutip.basis(2, 0)
    sigz = N.dynamics(
        H0, delta_V, omega, spectrum, T0, psi0, times, qutip.sigmaz(),
        vec_omega=wk, num_runs=NUM_RUNS,
    )

    # Fit (project convention: notebook's damped_rabi model, amp/offset fixed)
    try:
        popt, pcov = curve_fit(damped_rabi, times, sigz, p0=[T_Rabi * 5, 0.0],
                                maxfev=20000)
        t2star_fit = abs(popt[0])
        perr = np.sqrt(np.diag(pcov))[0]
        fit_ok = True
    except Exception as exc:
        t2star_fit, perr, fit_ok = np.nan, np.nan, False
        print(f"  [n_cycles={n_cycles}] fit FAILED: {exc}")

    envelope_at_end = np.exp(-(T / t2star_fit) ** 2) if fit_ok else np.nan
    # [VALIDATION] Fit-validity gate (project practical criterion, 지침 39절):
    # envelope가 0.7 미만으로 떨어져야 "resolved"
    resolved = fit_ok and envelope_at_end < 0.7
    status = "decay resolved" if resolved else "decay unresolved within simulation window"

    elapsed = time.time() - t_start
    results.append({
        "n_cycles": n_cycles,
        "T_sim_us": T * 1e6,
        "n_points": npts,
        "num_runs": NUM_RUNS,
        "t2star_fit_us": t2star_fit * 1e6 if fit_ok else np.nan,
        "t2star_fit_rel_uncertainty": (perr / t2star_fit) if fit_ok else np.nan,
        "ratio_fit_to_window": (t2star_fit / T) if fit_ok else np.nan,
        "envelope_at_window_end": envelope_at_end,
        "status": status,
        "wall_time_s": elapsed,
    })
    print(f"[n_cycles={n_cycles:4d}] T_sim={T*1e6:8.2f} us | "
          f"t2star_fit={t2star_fit*1e6:8.2f} us | "
          f"envelope@end={envelope_at_end:.4f} | {status} | ({elapsed:.1f}s)")

    ax.plot(times * 1e6, sigz, ".", ms=2, alpha=0.5, label=f"n_cycles={n_cycles} data")
    if fit_ok:
        ax.plot(times * 1e6, damped_rabi(times, *popt), "-", lw=1,
                label=f"n_cycles={n_cycles} fit (T={t2star_fit*1e6:.1f}us)")

ax.set_xlabel("t (us)")
ax.set_ylabel(r"$\langle \sigma_z(t) \rangle$")
ax.set_title("M08 decay-time convergence study (Parameter sensitivity / numerical validation)")
ax.legend(fontsize=7, ncol=2)
fig.tight_layout()
fig.savefig(export_dir / "m08_sim_time_convergence.png", dpi=150)
plt.close(fig)

df = pd.DataFrame(results)
print("\n" + "=" * 80)
print("CONVERGENCE TABLE")
print("=" * 80)
print(df.to_string(index=False))
df.to_csv(export_dir / "m08_sim_time_convergence_results.csv", index=False)

# [VALIDATION] 수렴 판정: 마지막 두 case의 변화율 확인 (project practical criterion)
if len(df) >= 2 and df["t2star_fit_us"].notna().iloc[-2:].all():
    last_two = df["t2star_fit_us"].iloc[-2:].values
    rel_change = abs(last_two[1] - last_two[0]) / abs(last_two[0])
    print(f"\n마지막 두 case(n_cycles={df['n_cycles'].iloc[-2]} -> "
          f"{df['n_cycles'].iloc[-1]}) 간 t2star_fit 변화율: {rel_change*100:.1f}%")
    if rel_change < 0.05 and df["envelope_at_window_end"].iloc[-1] < 0.7:
        print("판정: 수렴 가능성 있음 (상대변화<5% AND envelope<0.7) — 단, project-level "
              "practical criterion이며 QTCAD 공식 기준 아님.")
    else:
        print("판정: 아직 수렴하지 않음 — fit이 여전히 외삽(extrapolation)일 가능성이 높음.")

metadata = {
    "S0": S0, "wc_rad_s": float(wc), "w0": w0, "num_runs": NUM_RUNS,
    "points_per_cycle": POINTS_PER_CYCLE, "cycles_tested": cycles_list,
    "fit_model": "cos(omega_Rabi*t+phase)*exp(-(t/t2star)^2), amp=1, offset=0 (notebook과 동일)",
    "fit_validity_gate": "envelope_at_window_end < 0.7 (project practical criterion)",
    "classification": "B. Parameter sensitivity / numerical convergence study",
}
with open(export_dir / "m08_sim_time_convergence_metadata.json", "w", encoding="utf-8") as f:
    json.dump(metadata, f, indent=2)

print(f"\nExports: {export_dir}")
