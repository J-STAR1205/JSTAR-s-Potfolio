__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

import pathlib
import numpy as np
from numpy.typing import NDArray
from qtcad.device import constants as ct
from qtcad.device import materials as mt
from qtcad.device import Device, SubDevice
from qtcad.device.schrodinger_poisson import Solver
from qtcad.device.schrodinger_poisson import SolverParams
from qtcad.device.schrodinger import SolverParams as SchrodSolverParams
from qtcad.device.poisson import SolverParams as PoissonSolverParams
from qtcad.device import io
from double_dot_fdsoi import save_slice, state_probability_density

# Script directory
script_dir = pathlib.Path(__file__).parent.resolve()

# Degeneracy factor for silicon
# Start from the default confined-electron degeneracy and remove the
# valley factor, since valley splitting is treated explicitly.
g = int(mt.Si.gqc / 2)

# Solver parameters
# (Poisson)
p_params = PoissonSolverParams()
p_params.tol = 1e-8
# (Schrödinger)
schrod_params = SchrodSolverParams()
schrod_params.tol = 1e-8
schrod_params.num_states = 10
# (Schrödinger-Poisson)
sp_params = SolverParams()
sp_params.tol = 1e-5
sp_params.maxiter = 2000
sp_params.bound_state_charges_only = True
sp_params.sc_method = "underrelax"
sp_params.initialization = False
sp_params.mixing_algo = "adaptive_linear"
sp_params.mixing_param = 0.1
sp_params.poisson_solver_params = p_params
sp_params.schrod_solver_params = schrod_params


def compute_chemical_potential(
    V_SET: NDArray,
    Nmin: int,
    Nmax: int,
    d: Device,
    set_device: SubDevice,
    out_path: pathlib.Path = script_dir,
) -> NDArray:
    """Compute the chemical potential for different number of particles and
    different voltage configurations. This function also saves the potential
    to HDF5 files, eigenstate probability densities to VTU files, and the
    energies and population factors to text files. All output files are saved
    in the specified output directory.

    Args:
        V_SET: Array of SET gate voltages.
        Nmin: Minimum number of electrons.
        Nmax: Maximum number of electrons.
        d: Device object representing the entire FD-SOI device.
        set_device: SubDevice representing the SET region of the FD-SOI device.
        out_path: Output directory path.

    Returns:
        2D array of chemical potentials for each voltage configuration.
        The first index runs over the applied voltages to the SET gate, and the
        second index runs over the number of particles.

    Note:
        The chemical potential is calculated as the difference in many-body
        ground-state energies when adding an electron to the system. The
        chemical potentials are computed for the electron numbers in the range
        [Nmin+1, Nmax].
    """

    # output files
    # (results)
    phi_out_file = str(out_path / "phi_vset{vset:.2f}_N{n}.hdf5")
    mu_file = out_path / "mu.txt"
    pop_file = out_path / "population.txt"
    E_file = out_path / "energies.txt"

    # (visualization)
    wf_vtu_file = str(out_path / "wf_vset{vset:.2f}_N{n}.vtu")
    rho_slice_file = str(out_path / "rho_vset{vset:.2f}_N{n}_slice.png")
    wf_slice_file = str(out_path / "wf_vset{vset:.2f}_N{n}_slice.png")

    # Range of number of particles
    N = range(Nmin, Nmax + 1)
    # Many-body energies
    mb_E = np.zeros((len(V_SET), len(N)))
    # Chemical potential
    Mu = np.zeros((len(V_SET), len(N) - 1))

    for j, vset in enumerate(V_SET):
        # Print info
        print("\n")
        print("==================================================")
        print(f"Setting V_SET = {vset} V")
        print("==================================================")
        print("\n")

        # Set the potential in the SET region
        d.set_applied_potential("plunger_gate_1_bnd", vset)

        for i, n in enumerate(N):
            # Print info
            print("--------------------------------------------------")
            print(f"Solving for N = {n} electrons...")
            print("--------------------------------------------------")
            print("\n")

            # Solve Schrödinger-Poisson in the SET region
            sp_solver = Solver(d, subdevice=set_device, solver_params=sp_params)
            sp_solver.solve(N=n, g=g)
            print("Energies and population factors of SET SubDevice:")
            set_device.print_energies()

            # Save phi, eigenstate densities, and charge density
            io.save(phi_out_file.format(vset=vset, n=n), {"phi": d.phi})
            wf_dict = {
                f"WF_{i}": state_probability_density(set_device, i)
                for i in range(set_device.eigenfunctions.shape[1])
            }
            io.save(wf_vtu_file.format(vset=vset, n=n), wf_dict, set_device.mesh)

            save_slice(
                set_device,
                state_probability_density(set_device, 0),
                wf_slice_file.format(vset=vset, n=n),
                label="Ground-state probability density [1 / m$^3$]",
            )
            save_slice(
                d,
                d.rho / ct.e,
                rho_slice_file.format(vset=vset, n=n),
                label="$\\rho$ [$e$ / m$^3$]",
            )

            # Save energies and population factors
            with open(pop_file, "a") as f:
                out = [vset, n] + list(set_device.population_factors)
                np.savetxt(f, np.array(out)[np.newaxis, :], fmt="%.6f")
            with open(E_file, "a") as f:
                out = [vset, n] + list(set_device.energies / ct.e)
                np.savetxt(f, np.array(out)[np.newaxis, :], fmt="%.8f")

            # Compute many-body ground-state energy
            mb_E[j, i] = np.sum(g * set_device.population_factors * set_device.energies)
            print(f"  Many-body ground-state energy (eV): {mb_E[j, i] / ct.e:.10f}")

            if i > 0:
                # Compute chemical potential
                Mu[j, i - 1] = mb_E[j, i] - mb_E[j, i - 1]
                print(f"  Chemical potentials (eV): {Mu[j, i - 1] / ct.e}")
                with open(mu_file, "a") as f:
                    np.savetxt(
                        f, np.array([[vset, n, Mu[j, i - 1] / ct.e]]), fmt="%.10f"
                    )

    return Mu
