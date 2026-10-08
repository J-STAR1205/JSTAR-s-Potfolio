__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import pathlib
import numpy as np
from matplotlib import pyplot as plt
from qtcad.device import constants as ct
import qutip
from qtcad.qubit import spectra, dynamics, noise

# [VALIDATION] [수정 — M08 Step 1, stale-data pipeline fix]
# 과거에는 이 자리에 h0, u를 직접 하드코딩했는데(손으로 복사한 값), MOS_EDSR.py의
# 파라미터가 바뀐 뒤에도 갱신되지 않아 Rabi 주파수가 1.6899MHz(현재 MOS_EDSR.py 실행)
# vs 1.5177MHz(과거 하드코딩) 로 약 11% 불일치하는 stale-data 문제가 있었다
# (M03_M06_M08_Physics_Review_Report.md 7절 — RWA 컨벤션 차이가 아니라는 것을
# 재현 테스트로 확인했음). 이제 h0, u는 MOS_EDSR.py가 저장한 파일에서 직접 로드한다
# — 사람이 숫자를 복사해 붙여넣는 과정이 없으므로 이 종류의 불일치는 구조적으로
# 재발하지 않는다.
script_dir = pathlib.Path(__file__).parent.resolve()
h0_u_file = script_dir / "output" / "mos_edsr_h0_u.npz"
if not h0_u_file.exists():
    raise FileNotFoundError(
        f"{h0_u_file} 가 없습니다. 먼저 MOS_EDSR.py를 실행해 h0/u를 생성해야 합니다 "
        "(EDSR_noise.py는 더 이상 h0/u를 하드코딩하지 않습니다)."
    )
_loaded = np.load(h0_u_file)
h0 = _loaded["h0"]
H0 = h0 / ct.hbar  # Rewrite in units of hbar = 1 ([NUMERICAL] h0 단위: J)
u = _loaded["u"]
delta_V = u / ct.hbar  # Rewrite in units of hbar = 1 ([NUMERICAL] u 단위: J)
omega0_loaded = float(_loaded["omega0"])  # rad/s, MOS_EDSR.py의 transition_2_levels() 반환값
omega_rabi_loaded = float(_loaded["omega_rabi"])  # rad/s, 동일

omega = np.absolute(H0[0, 0] - H0[1, 1])  # qubit frequency (rad/s), h0로부터 재계산
omega_Rabi = np.absolute(delta_V[0, 1])  # Rabi frequency on resonance (rad/s), u로부터 재계산

# [VALIDATION] Consistency check (M06.md 36절 acceptance criterion):
# h0/u로부터 재계산한 omega/omega_Rabi가 MOS_EDSR.py가 저장한 omega0/omega_rabi와
# 1e-3 상대오차 이내로 일치해야 한다. 불일치하면(예: npz가 손상되었거나 다른 실행의
# 파일과 혼동된 경우) noise simulation을 진행하지 않고 즉시 중단한다.
_rel_err_omega = abs(omega - omega0_loaded) / abs(omega0_loaded)
_rel_err_rabi = abs(omega_Rabi - omega_rabi_loaded) / abs(omega_rabi_loaded)
print(f"[consistency check] qubit omega: recomputed={omega:.6e} rad/s, "
      f"loaded={omega0_loaded:.6e} rad/s, rel.err={_rel_err_omega:.2e}")
print(f"[consistency check] Rabi omega : recomputed={omega_Rabi:.6e} rad/s, "
      f"loaded={omega_rabi_loaded:.6e} rad/s, rel.err={_rel_err_rabi:.2e}")
if _rel_err_omega > 1e-3 or _rel_err_rabi > 1e-3:
    raise RuntimeError(
        f"h0/u consistency check FAILED (qubit rel.err={_rel_err_omega:.2e}, "
        f"Rabi rel.err={_rel_err_rabi:.2e}, threshold=1e-3). "
        "mos_edsr_h0_u.npz가 손상되었거나 MOS_EDSR.py의 다른 실행 결과와 섞였을 "
        "가능성이 있습니다. Noise simulation을 중단합니다."
    )
print("[consistency check] PASSED — proceeding with noise simulation.\n")

T_Rabi = 2 * np.pi / omega_Rabi  # Rabi period

T = 20 * T_Rabi  # Total time of the simulations
times = np.linspace(0, T, 1000)  # Times at which dynamics are simulated.

# Initialize objects that will handle dynamics and noise.
dyn = dynamics.Dynamics()
N = noise.Noise()

# The Hamiltonian H = H0 + delta_V cos(omega*t) in a frame rotating
# at frequency omega.
H = dyn.H_RF(H0, delta_V, omega)
# Solve for the dynamics without any noise.
psi0 = qutip.basis(2, 0)  # initial state
result0 = qutip.mesolve(H, psi0, times)

# Plot expectation values of sigma_z and sigma_y as a function of time.
fig, axes = plt.subplots(2)
fig.set_size_inches((8, 8))
# <sigma_z (t)>
axes[0].plot(result0.times, qutip.expect(qutip.sigmaz(), result0.states))
# <sigma_y (t)>
axes[1].plot(result0.times, qutip.expect(qutip.sigmay(), result0.states))
# Labels
axes[0].set_ylabel(f"$\left< \sigma_z (t) \\right>$", fontsize=20)
axes[1].set_ylabel(f"$\left< \sigma_y (t) \\right>$", fontsize=20)
axes[0].set_xlabel(f"$t \, (s)$", fontsize=20)
axes[1].set_xlabel(f"$t \, (s)$", fontsize=20)
plt.show()

# Include charge noise which is assumed to affect only the gate which
# generates the drive delta_V.
# Lorentzian spectrum
S0 = 0.01  # Total power
wc = omega_Rabi  # Cutoff frequency
w0 = 0  # Central frequency
# T0 should be much larger than the largest time scale of the problem under consideration.
# Used to properly define the Fourier series of the time series related to the noise.
T0 = 50 * T_Rabi  # Must take T0 such that 2*pi/T0 << wc to resolve the spectrum
omega_max = 5 * wc
wk = np.arange(0, omega_max, 2 * np.pi / T0) + w0
spectrum = spectra.lorentz(wk, S0, wc, w0=w0)

# Plot spectrum
fig, axes = plt.subplots(1)
fig.set_size_inches((10, 6))
axes.plot(wk, spectrum)
axes.set_ylabel(f"$S_\\beta(\omega) \, (s)$", fontsize=20)
axes.set_xlabel("$\omega \, (s^{-1})$", fontsize=20)
axes.set_title("Lorentzian spectral density", fontsize=30)
plt.show()

# Plot time series for the charge noise affecting the gate
# responsible for delta_V.
beta = N.gen_process(times, spectrum, wk, T0)
fig, axes = plt.subplots(1)
fig.set_size_inches((8, 5))
axes.plot(times, beta)
axes.set_ylabel(f"$\\beta(t)$", fontsize=20)
axes.set_xlabel("$t \, (s)$", fontsize=20)
plt.show()

# Simulate dynamic including noise.
num_runs1 = 1  # number of runs to average over
sigz_1run = N.dynamics(
    H0,
    delta_V,
    omega,
    spectrum,
    T0,
    psi0,
    times,
    qutip.sigmaz(),
    vec_omega=wk,
    num_runs=num_runs1,
)
num_runs2 = 100  # number of runs to average over
sigz_multiruns = N.dynamics(
    H0,
    delta_V,
    omega,
    spectrum,
    T0,
    psi0,
    times,
    qutip.sigmaz(),
    vec_omega=wk,
    num_runs=num_runs2,
)

# Plot dynamics
fig, axes = plt.subplots(2)
axes[0].plot(
    result0.times, qutip.expect(qutip.sigmaz(), result0.states), label="Without noise"
)
axes[0].plot(result0.times, sigz_1run, label="With noise")
axes[1].plot(
    result0.times, qutip.expect(qutip.sigmaz(), result0.states), label="Without noise"
)
axes[1].plot(result0.times, sigz_multiruns, label="With noise")
# Format Plots
axes[0].set_ylabel(f"$\left< \sigma_z (t) \\right>$", fontsize=20)
axes[1].set_ylabel(f"$\left< \sigma_z (t) \\right>$", fontsize=20)
axes[0].set_xlabel(f"$t \, (s)$", fontsize=20)
axes[1].set_xlabel(f"$t \, (s)$", fontsize=20)
axes[0].set_title(f"num_runs = {num_runs1}")
axes[1].set_title(f"num_runs = {num_runs2}")
axes[0].legend(loc="center left", fontsize=16)
axes[1].legend(loc="center left", fontsize=16)
fig.tight_layout()

plt.show()
