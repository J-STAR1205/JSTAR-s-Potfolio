__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

# [VALIDATION] M08 2차 라운드 Step 3 — Noise ensemble 수렴성 + realization 독립성 검증
#
# 분류: B. Parameter sensitivity / numerical validation study.
#   EDSR_noise.py는 수정하지 않는다.
#
# [NUMERICAL] 비용을 낮추기 위해 20 Rabi cycle 창(가장 빠른 case, 이전 실험에서
#   100-run에 6.9s)을 사용한다 — 여기서는 decay fit이 아니라 "ensemble 크기가
#   평균/분산에 미치는 영향"만 보는 것이 목적이므로 긴 창이 필요 없다.
#
# [WARNING] 각 N마다 단 1회의 ensemble만 보면 "N=100이 충분한지"를 판단할 근거가
#   부족하다 (그 자체가 하나의 확률표본일 뿐). R=5회 반복해서 ensemble-mean 자체의
#   run-to-run 변동(=실질적 standard error)을 직접 측정한다.

import json
import pathlib
import time

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import qutip
from qtcad.device import constants as ct
from qtcad.qubit import spectra, dynamics, noise

script_dir = pathlib.Path(__file__).parent.resolve()
output_dir = script_dir / "output"
export_dir = output_dir / "analysis_exports"
export_dir.mkdir(parents=True, exist_ok=True)

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

S0 = 0.01
wc = omega_Rabi
w0 = 0
N_CYCLES = 20  # [NUMERICAL] 빠른 창으로 ensemble-size 효과만 분리해서 봄
T = N_CYCLES * T_Rabi
npts = int(50 * N_CYCLES)
times = np.linspace(0, T, npts)
T0 = max(50 * T_Rabi, 2 * T)
omega_max = 5 * wc
wk = np.arange(0, omega_max, 2 * np.pi / T0) + w0
spectrum = spectra.lorentz(wk, S0, wc, w0=w0)
psi0 = qutip.basis(2, 0)

N_list = [1, 10, 50, 100]
R_REPEATS = 5  # 각 N마다 5회 독립 반복 -> ensemble-mean의 run-to-run 변동 측정

print("=" * 80)
print("COMPUTATIONAL COST ESTIMATE")
print("=" * 80)
total_calls = sum(N_list) * R_REPEATS
print(f"Total noise realizations to compute: {total_calls} "
      f"(sum(N)={sum(N_list)} x R={R_REPEATS})")
print("=" * 80 + "\n")

all_results = {}  # N -> list of R waveforms
t0_total = time.time()
for Nval in N_list:
    waveforms = []
    for r in range(R_REPEATS):
        t0 = time.time()
        sigz = noise.Noise().dynamics(H0, delta_V, omega, spectrum, T0, psi0, times,
                                       qutip.sigmaz(), vec_omega=wk, num_runs=Nval)
        waveforms.append(sigz)
        elapsed = time.time() - t0
        print(f"[N={Nval}, repeat {r+1}/{R_REPEATS}] computed in {elapsed:.1f}s")
    all_results[Nval] = np.array(waveforms)  # shape (R_REPEATS, npts)

print(f"\nTotal wall time: {time.time()-t0_total:.1f}s")

# --- N=100 (가장 큰 ensemble)의 R회 반복 평균을 '기준'으로 사용 ---
reference_mean = all_results[100].mean(axis=0)

rows = []
for Nval in N_list:
    waveforms = all_results[Nval]  # (R, npts)
    mean_waveform = waveforms.mean(axis=0)
    std_across_repeats = waveforms.std(axis=0)  # time-pointwise std across R repeats
    sem_estimate = std_across_repeats  # 이것 자체가 "ensemble mean의 R-반복 표준오차"
    max_diff_vs_ref = np.max(np.abs(mean_waveform - reference_mean))
    rms_diff_vs_ref = np.sqrt(np.mean((mean_waveform - reference_mean) ** 2))
    rows.append({
        "N_realizations": Nval,
        "mean_of_mean_sigz_at_t0": mean_waveform[0],
        "max_pointwise_std_across_repeats": np.max(std_across_repeats),
        "mean_pointwise_std_across_repeats": np.mean(std_across_repeats),
        "max_abs_diff_vs_N100": max_diff_vs_ref,
        "rms_diff_vs_N100": rms_diff_vs_ref,
    })
    print(f"[N={Nval}] max|mean-mean_N100|={max_diff_vs_ref:.4f}, "
          f"RMS diff={rms_diff_vs_ref:.4f}, "
          f"mean pointwise std(across {R_REPEATS} repeats)={np.mean(std_across_repeats):.4f}")

conv_df = pd.DataFrame(rows)
print("\n" + "=" * 80)
print(conv_df.to_string(index=False))
conv_df.to_csv(export_dir / "m08_ensemble_convergence.csv", index=False)

# --- realization 독립성 검증: num_runs=1을 여러 번 호출해서 서로 다른지 확인 ---
print("\n" + "=" * 80)
print("REALIZATION INDEPENDENCE CHECK (num_runs=1, called 3 times)")
print("=" * 80)
single_runs = []
for i in range(3):
    sigz1 = noise.Noise().dynamics(H0, delta_V, omega, spectrum, T0, psi0, times,
                                    qutip.sigmaz(), vec_omega=wk, num_runs=1)
    single_runs.append(sigz1)
for i in range(len(single_runs)):
    for j in range(i + 1, len(single_runs)):
        diff = np.max(np.abs(single_runs[i] - single_runs[j]))
        print(f"max|run_{i} - run_{j}| = {diff:.6f} "
              f"({'서로 다름 (정상)' if diff > 1e-6 else '동일 (의심 — seed 재사용?)'})")

# --- Plot: N별 waveform 비교 ---
fig, axs = plt.subplots(1, 2, figsize=(13, 5))
for Nval, c in zip(N_list, ["tab:gray", "tab:orange", "tab:green", "tab:blue"]):
    mean_waveform = all_results[Nval].mean(axis=0)
    axs[0].plot(times * 1e9, mean_waveform, label=f"N={Nval}", color=c, lw=1.2)
axs[0].set_xlabel("t (ns)"); axs[0].set_ylabel(r"$\langle\sigma_z(t)\rangle$")
axs[0].set_title(f"Ensemble-mean waveform vs. N ({N_CYCLES} cycles)"); axs[0].legend(fontsize=8)

axs[1].plot(conv_df["N_realizations"], conv_df["rms_diff_vs_N100"], "o-", label="RMS diff vs N=100")
axs[1].plot(conv_df["N_realizations"], conv_df["mean_pointwise_std_across_repeats"], "s-",
            label=f"mean pointwise std ({R_REPEATS} repeats)")
axs[1].set_xlabel("N (ensemble size)"); axs[1].set_ylabel("magnitude")
axs[1].set_xscale("log"); axs[1].set_title("Ensemble convergence metrics vs. N")
axs[1].legend(fontsize=8); axs[1].grid(True)
fig.tight_layout()
fig.savefig(export_dir / "m08_ensemble_convergence.png", dpi=150)
plt.close(fig)

print(f"\nExports: {export_dir}")
