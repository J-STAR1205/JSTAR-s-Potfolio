"""M07 Study 1단계 — 1D MVEMT로 Si/SiGe 우물의 valley splitting

무엇을 바꾸나 (Reference 대비)
------------------------------
Reference(Device 12)는 Si/SiO2 계면 하나(장벽 3 eV, 계단형)였다.
여기서는 Philips 2022 적층에 가깝게
  - 장벽 높이: 3 eV → Si/Si0.7Ge0.3 전도대 오프셋 ΔEc (기본 0.15 eV)
  - 계면 수: 1개 → 2개 (Si 우물 두께 WELL_WIDTH_NM)
  - 전기장: 5~50 MV/m (MOS) → 1~10 MV/m (Si/SiGe에서 흔한 범위)
  - 계면 폭: 계단 → tanh 프로파일, 폭 w를 스윕
로 바꾸고, valley splitting(VS)이 계면 폭과 전기장에 어떻게 반응하는지 본다.

물리 요약
------------------------------
VS는 대략 "계면에서 envelope가 느끼는 퍼텐셜의 2k0 푸리에 성분"에 비례한다
(k0 ≈ 0.84·2π/a, 주기 약 0.32 nm). 계면이 계단이면 이 성분이 크고,
계면이 수 원자층에 걸쳐 완만해지면 2k0 성분이 급격히 줄어 VS가 작아진다.
MVEMT는 평균 프로파일만 보므로 합금 무질서에 의한 VS 산포는 나오지 않는다
(그건 3단계 TB에서 본다).

좌표·부호 규약 (Reference와 동일)
------------------------------
1D 메시의 축은 x이고, 이것이 성장(구속) 방향이다 → ±x valley가 갈라진다.
퍼텐셜 에너지 V(x) = 장벽 - e·E·x 이므로 E > 0이면 전자는 +x 쪽,
즉 상부 계면(x = 0, 스페이서 쪽)으로 눌린다.

  x:  -(W+B) ── buffer ── -W ── Si well ── 0 ── spacer ── S

근사
------------------------------
Reference에서 ff 근사는 full 대비 평균 오차 1.9 %, 중앙값 기준 약 2.8배 빠름 →
스윕은 ff로 한다. 장벽 영역도 Si 유효질량·Bloch 함수를 쓰는 근사다
(SiGe 장벽 안의 파동함수 꼬리가 작으므로 VS에 대한 영향은 작다고 가정).

실행: (qtcad) python 01_mvemt_sige_well.py
출력: output/01_mvemt_vs.csv, output/01_vs_vs_E.png, output/01_vs_vs_width.png
예상 시간: 1점당 수십 초~수 분 × (검증 1점 + 폭 4 × 전기장 4 = 17점)
"""

import pathlib
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from qtcad.device.mesh1d import Mesh
from qtcad.device import Device
from qtcad.device import materials as mt
from qtcad.device import constants as ct
from qtcad.device.valleycoupling.bloch import silicon as bloch_amplitudes
from qtcad.device.multivalley_EMT import Solver, SolverParams

from atoms_builder import well_window

# =============================================================================
# 바꿔볼 파라미터 (모두 여기서만 수정)
# =============================================================================
WELL_WIDTH_NM = 8.0          # Si 우물 두께 (Philips 2022: 8 nm)
SPACER_NM = 6.0              # 우물 위 SiGe (1D 계산 영역용, 실제 30 nm 중 일부)
BUFFER_NM = 12.0             # 우물 아래 SiGe
DEC_EV = 0.15                # 전도대 오프셋 ΔEc (eV) — QTCAD 재료값으로 확인할 것
E_FIELDS_MVM = [1.0, 2.0, 5.0, 10.0]        # 수직 전기장 (MV/m)
INTERFACE_WIDTHS_NM = [0.0, 0.25, 0.5, 1.0]  # tanh 특성 길이 w (10-90 % 폭 ≈ 2.2 w)
APPROX = "ff"                # "trivial" / "ff" / None(full)
H_NM = 0.01                  # 메시 간격 (Reference와 동일). 계단 계면엔 0.01 이하 권장
NUM_STATES = 6               # VS에는 2개면 충분, 여유 있게 6
RUN_VALIDATION = True        # Reference 재현 1점(ΔEc=3 eV, 계단, E=50 MV/m, ff) 먼저 실행

# =============================================================================
# 경로·상수
# =============================================================================
script_dir = pathlib.Path(__file__).parent.resolve()
mesh_dir = script_dir / "meshes"
out_dir = script_dir / "output"
mesh_dir.mkdir(exist_ok=True)
out_dir.mkdir(exist_ok=True)
path_mesh = mesh_dir / "qw_1d_SiGe.msh"
path_csv = out_dir / "01_mvemt_vs.csv"

scale = 1e-9
a = 0.543 * scale
valleys = mt.Si.valleys[0:6]
me_inv = mt.Si.Me_inv_valley[0:6]
bloch_path = pathlib.Path(bloch_amplitudes.__path__[0])
bloch_paths = [bloch_path / f for f in
               ["BA_px_Si.data", "BA_mx_Si.data", "BA_py_Si.data",
                "BA_my_Si.data", "BA_pz_Si.data", "BA_mz_Si.data"]]

REF_VS_FF_50 = 1.0372  # meV, Reference ff 근사 E = 50 MV/m 결과


# =============================================================================
# 1. 메시: buffer / well / spacer 세 구간
# =============================================================================
def build_mesh(path):
    W, S, B, h = WELL_WIDTH_NM, SPACER_NM, BUFFER_NM, H_NM
    xs = [-(W + B), -W, 0.0, S]
    names = ["Buffer", "Well", "Spacer"]
    try:
        import gmsh
        gmsh.initialize()
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.model.add("qw_1d_SiGe")
        pts = [gmsh.model.geo.addPoint(xv, 0, 0, h) for xv in xs]
        lines = [gmsh.model.geo.addLine(pts[i], pts[i + 1]) for i in range(3)]
        gmsh.model.geo.synchronize()
        for tag, (name, ln) in enumerate(zip(names, lines), start=1):
            gmsh.model.addPhysicalGroup(1, [ln], tag)
            gmsh.model.setPhysicalName(1, tag, name)
        gmsh.model.mesh.generate(1)
        gmsh.write(str(path))
        gmsh.finalize()
    except ImportError:
        # gmsh 파이썬 모듈이 없으면 .geo를 쓰고 gmsh 실행 파일로 메시 생성
        import subprocess
        geo = path.with_suffix(".geo")
        body = [f"h = {h};"]
        body += [f"Point({i+1}) = {{{xv}, 0, 0, h}};" for i, xv in enumerate(xs)]
        body += [f"Line({i+1}) = {{{i+1}, {i+2}}};" for i in range(3)]
        body += [f'Physical Line("{n}", {i+1}) = {{{i+1}}};' for i, n in enumerate(names)]
        body += ["Mesh 1;", f'Save "{path.name}";']
        geo.write_text("\n".join(body) + "\n")
        subprocess.run(["gmsh", geo.name, "-"], cwd=geo.parent, check=True)


# =============================================================================
# 2. 소자와 퍼텐셜
# =============================================================================
def create_device(E_Vm, w_nm, dEc_eV):
    mesh = Mesh(scale, str(path_mesh))
    d = Device(mesh, conf_carriers="e")
    for name in ["Buffer", "Well", "Spacer"]:
        d.new_region(name, mt.Si)   # 장벽은 아래 V(x)로 넣는다

    W = WELL_WIDTH_NM * scale
    w = w_nm * scale

    def V(x):
        barrier = dEc_eV * ct.e * (1.0 - well_window(x, 0.0, -W, w))
        return barrier - ct.e * E_Vm * x

    d.set_V(V)
    return d


def solve(d, approx):
    p = SolverParams()
    p.num_states = NUM_STATES
    p.approx = approx
    p.val = valleys
    p.bloch_files = bloch_paths
    p.m_inv_tensors = me_inv
    p.lattice_const = a
    p.method = "fast"
    p.tol = 1e-7
    p.maxiter = 5000
    t0 = time.perf_counter()
    Solver(d, solver_params=p).solve()
    return time.perf_counter() - t0


def ground_state_position_nm(d):
    """바닥 상태 확률밀도(6 valley 합)의 평균 위치 — 부호 규약 확인용."""
    x = d.mesh.glob_nodes[:, 0] / scale
    rho = np.sum(np.abs(d.eigenfunctions[:, 0, :]) ** 2, axis=1)
    return float(np.sum(x * rho) / np.sum(rho))


def run_point(E_mvm, w_nm, dEc_eV, approx):
    d = create_device(E_mvm * 1e6, w_nm, dEc_eV)
    dt = solve(d, approx)
    vs = (d.energies[1] - d.energies[0]) / ct.e * 1e3
    return {"E_MVm": E_mvm, "w_nm": w_nm, "dEc_eV": dEc_eV,
            "approx": str(approx), "VS_meV": vs,
            "E0_meV": d.energies[0] / ct.e * 1e3,
            "x_mean_nm": ground_state_position_nm(d), "time_s": dt}


# =============================================================================
# 3. 실행
# =============================================================================
build_mesh(path_mesh)
rows = []

if RUN_VALIDATION:
    r = run_point(50.0, 0.0, 3.0, "ff")
    print(f"[검증] VS = {r['VS_meV']:.4f} meV (Reference ff: {REF_VS_FF_50:.4f} meV), "
          f"<x> = {r['x_mean_nm']:.2f} nm, {r['time_s']:.1f} s")
    print("       → 수 % 안이면 메시·퍼텐셜 설정이 Reference와 일치. "
          "<x>가 0에 가까운 음수(상부 계면 바로 아래)여야 부호가 맞다.")
    r["case"] = "validation"
    rows.append(r)

for w_nm in INTERFACE_WIDTHS_NM:
    for E in E_FIELDS_MVM:
        r = run_point(E, w_nm, DEC_EV, APPROX)
        r["case"] = "sweep"
        rows.append(r)
        print(f"w = {w_nm:.2f} nm, E = {E:5.1f} MV/m → VS = {r['VS_meV']:.4f} meV, "
              f"<x> = {r['x_mean_nm']:.2f} nm, {r['time_s']:.1f} s")
        pd.DataFrame(rows).to_csv(path_csv, index=False)  # 중간 저장

df = pd.DataFrame(rows)
df.to_csv(path_csv, index=False)
sw = df[df["case"] == "sweep"]

# =============================================================================
# 4. 그림
# =============================================================================
fig, ax = plt.subplots(figsize=(7, 4.5))
for w_nm, g in sw.groupby("w_nm"):
    ax.plot(g["E_MVm"], g["VS_meV"] * 1e3, "-o", label=f"w = {w_nm:g} nm")
ax.set_xlabel("Electric field (MV/m)")
ax.set_ylabel("Valley splitting (µeV)")
ax.set_title(f"1D MVEMT, Si/SiGe {WELL_WIDTH_NM:g} nm well, ΔEc = {DEC_EV} eV ({APPROX})")
ax.legend(); ax.grid(True)
fig.tight_layout(); fig.savefig(out_dir / "01_vs_vs_E.png", dpi=150)

fig, ax = plt.subplots(figsize=(7, 4.5))
for E, g in sw.groupby("E_MVm"):
    ax.semilogy(g["w_nm"], g["VS_meV"] * 1e3, "-o", label=f"E = {E:g} MV/m")
ax.set_xlabel("Interface width parameter w (nm)  [10-90 % ≈ 2.2 w]")
ax.set_ylabel("Valley splitting (µeV)")
ax.set_title("Valley splitting vs. interface width")
ax.legend(); ax.grid(True, which="both")
fig.tight_layout(); fig.savefig(out_dir / "01_vs_vs_width.png", dpi=150)

print("\n요약 (µeV):")
print((sw.pivot(index="w_nm", columns="E_MVm", values="VS_meV") * 1e3).round(2))
print(f"\n중앙값 계산 시간: {sw['time_s'].median():.1f} s/점")
plt.show()
