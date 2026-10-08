__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import pathlib
import time
import numpy as np
import matplotlib.pyplot as plt
from qtcad.device.mesh1d import Mesh
from qtcad.device import Device
from qtcad.device import materials as mt
from qtcad.device import constants as ct

from qtcad.device.valleycoupling.bloch import silicon as bloch_amplitudes
from qtcad.device.multivalley_EMT import Solver, SolverParams

# Setup -----------------------------------------------------------------------
# File paths
script_dir = pathlib.Path(__file__).parent.resolve()
path_mesh = script_dir / "meshes" / "Si_well.msh"
path_out = script_dir / "output"

# Outputs
## Plots
path_valley = path_out / f"valley_wavefunctions.png"
path_vs_E = path_out / f"valley_splitting_vs_E.png"
## Data
path_triv = path_out / f"valley_splitting_qw_trivial.txt"
path_ff = path_out / f"valley_splitting_qw_ff.txt"
path_None = path_out / f"valley_splitting_qw_None.txt"
## [수정] 점별 계산 시간 기록 (근사 수준, E, 풀이 시간)
path_timing = path_out / "valley_splitting_timing.txt"

# [수정] 결과 파일은 아래 루프에서 덧쓰기("a")로 저장되므로,
# 재실행 시 행이 누적되지 않도록 시작할 때 기존 파일을 지운다.
path_out.mkdir(exist_ok=True)
for _p in (path_triv, path_ff, path_None, path_timing):
    _p.unlink(missing_ok=True)

# Inputs
## Bloch amplitudes
bloch_path = pathlib.Path(bloch_amplitudes.__path__[0])
bloch_paths = [
    bloch_path / "BA_px_Si.data",
    bloch_path / "BA_mx_Si.data",
    bloch_path / "BA_py_Si.data",
    bloch_path / "BA_my_Si.data",
    bloch_path / "BA_pz_Si.data",
    bloch_path / "BA_mz_Si.data",
]

## Constants
scale = 1e-9  # nm
a = 0.543 * scale  # Lattice constant of Si
k0 = 0.84 * 2 * np.pi / a  # Wavevector
valleys = mt.Si.valleys[0:6]  # Valley wavevectors
me_inv = mt.Si.Me_inv_valley[0:6]  # Inverse effective mass tensors


# Create the device ------------------------------------------------------------
def create_well(E) -> Device:
    """Creates a quantum-well device with a given electric field.

    Args:
        E (float): Electric field in V/m.

    Returns:
        Device: Device object.
    """

    # Mesh
    mesh = Mesh(scale, path_mesh)

    # Device
    d = Device(mesh, conf_carriers="e")
    d.new_region("Domain", mt.Si)
    d.new_region("Barrier", mt.SiO2)

    # Set potential
    def V(z):
        v = -ct.e * E * z

        return v

    d.set_V(V)
    d.add_to_V(3 * ct.e, "Barrier")  # Add barrier offset

    return d


# Solving MVEMT Schrödinger equation -------------------------------------------
def solve_MVEMT(d: Device, approx: str = None) -> None:
    """Solve the multivalley effective mass theory for the device.

    Args:
        d (Device): Device for which we wish to solve the MVEMT.
        approx (None or str): Level of approximation to use in the MVEMT.
            - "trivial" for trivial approximation.
            - "ff" for form-factor approximation.
            - None for full calculation.
    """

    # Compute valley coupling (using form-factor approximation)
    params = SolverParams()
    params.num_states = 20
    params.approx = approx
    params.val = valleys
    params.bloch_files = bloch_paths
    params.m_inv_tensors = me_inv
    params.lattice_const = a
    params.method = "fast"
    params.tol = 1e-7
    params.maxiter = 5000

    mvemt = Solver(d, solver_params=params)
    mvemt.solve()
    d.print_energies()


# Calculation -----------------------------------------------------------------

# Levels of approximation
approximations = ["trivial", "ff", None]
for approx in approximations:
    # Electric field values
    E_fields = np.linspace(5, 50, 10)
    for E in E_fields:
        d = create_well(E * 1e6)
        t0 = time.perf_counter()  # [수정] 계산 시간 측정
        solve_MVEMT(d, approx)
        dt = time.perf_counter() - t0
        with open(path_timing, "a") as file:  # [수정]
            file.write(f"{approx}\t{E:.3f}\t{dt:.2f}\n")

        save_file = path_out / f"valley_splitting_qw_{approx}.txt"
        with open(save_file, "a") as file:
            out = np.zeros((1, 2))
            out[0, 0] = E
            out[0, 1] = (d.energies[1] - d.energies[0]) / ct.e
            np.savetxt(file, out)

# Visualizing the results ------------------------------------------------------


def plot_env(wf: np.ndarray, mesh: Mesh, path_fig: pathlib.Path) -> None:
    """Plot the envelope contributions to the wavefunctions from each valley.

    Args:
        wf (np.ndarray): Wavefunctions computed from a MVEMT solver.
        mesh (Mesh): Mesh object over which the wavefunctions are defined.
        path_fig (pathlib.Path): Path to save the plot.
    """
    # Sort the wavefunctions over the global nodes of the mesh
    ix = np.argsort(mesh.glob_nodes, axis=0)[:, 0]  # index
    x = mesh.glob_nodes[ix, 0] / scale  # sorted nodes
    env = wf[ix, :, :]  # sorted wf

    # Labels for the valley states
    labels = ["+x", "-x", "+y", "-y", "+z", "-z"]

    # Create the figure
    fig, axs = plt.subplots(3, 2, figsize=(10, 10))

    # Plotting the valley contributions
    for i in range(6):
        ax = axs[i // 2, i % 2]
        ax2 = ax.twinx()
        for j in range(env.shape[1]):
            ax.plot(
                x, np.real(env[:, j, i]), label=f"Re[$F_{'{' + labels[i] + '}'}^{j}$]"
            )
            ax.plot(
                x, np.imag(env[:, j, i]), label=f"Im[$F_{'{' + labels[i] + '}'}^{j}$]"
            )
            ax2.plot(
                x,
                np.abs(env[:, j, i]) ** 2,
                "--",
                label=f"$|F_{'{' + labels[i] + '}'}^{j}|^2$",
            )
        ax.set_title(f"{labels[i]} valley")
        # [참고] 1D 메시의 좌표축이 x이므로 구속 방향도 x → ±x valley가 갈라진다
        ax.set_xlabel("$x$ [nm] (confinement axis)")
        ax.set_ylabel(f"$F_{'{' + labels[i] + '}'}$ [1/m$^{{3/2}}$]")
        ax2.set_ylabel(f"$|F_{'{' + labels[i] + '}'}|^2$ [1/m$^{{3}}$]")
        ax.legend(loc="upper left")
        ax2.legend(loc="lower left")
        ax.grid()

    fig.tight_layout()
    plt.savefig(path_fig)


plot_env(d.eigenfunctions[:, 0:2, :], d.mesh, path_valley)


def plot_vs(
    data_None: np.ndarray,
    data_ff: np.ndarray,
    data_triv: np.ndarray,
    save_file: pathlib.Path,
) -> None:
    """Plot the valley splitting for the different approximations.
    Args:
        data_None (np.ndarray): Valley splitting data for full calculation.
        data_ff (np.ndarray): Valley splitting data for form-factor
            approximation.
        data_triv (np.ndarray): Valley splitting data for trivial
            approximation.
        save_file (pathlib.Path): Path to save the plot.
    """

    # Create Figure
    fig, ax = plt.subplots(nrows=1, ncols=1, layout="tight")
    fig.set_size_inches(8, 4)

    # Plot data in [meV]
    ax.plot(
        data_ff[:, 0], data_ff[:, 1] * 1000, "-r", linewidth=2, label="form factors"
    )
    ax.plot(data_triv[:, 0], data_triv[:, 1] * 1000, "-k", linewidth=2, label="trivial")
    ax.plot(data_None[:, 0], data_None[:, 1] * 1000, "-g", linewidth=2, label="full")

    # Settings
    ax.set_title(f"Fig. (S1) [Gamble et al. APL 109 (2016)] using QTCAD", fontsize=18)
    ax.set_xlabel("Electric field [MV/m]", fontsize=16)
    ax.set_ylabel("valley splitting [meV]", fontsize=16)

    ax.tick_params(axis="x", labelsize=14)
    ax.tick_params(axis="y", labelsize=14)

    ax.legend(fontsize=16)
    ax.grid()

    # Save
    plt.savefig(save_file, dpi=300)


# Reload Data
data_None = np.loadtxt(path_None)
data_ff = np.loadtxt(path_ff)
data_triv = np.loadtxt(path_triv)

# Plot valley splitting vs. E
plot_vs(data_None, data_ff, data_triv, path_vs_E)
