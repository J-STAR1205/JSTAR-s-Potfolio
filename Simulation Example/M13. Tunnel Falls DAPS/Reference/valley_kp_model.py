__copyright__ = "Copyright 2022-2026, Nanoacademic Technologies Inc."

from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import RegularGridInterpolator
from qtcad.device import ElectronKPModel
from qtcad.device import ElectronKPParameter
from qtcad.device import constants as ct


def gaussian_covariance(
    positions: np.ndarray,
    *,
    sigma_delta_complex: float,
    correlation_length_x: float,
    correlation_length_y: float,
) -> np.ndarray:
    """Return the covariance for one component of the complex field Delta.

    Args:
        positions: Array of ``(x, y)`` positions where the field is sampled.
        sigma_delta_complex: Square root of the variance of the complex
            valley-coupling field, defined so that
            ``<|Delta|**2> = sigma_delta_complex**2``.
        correlation_length_x: Correlation length along x.
        correlation_length_y: Correlation length along y.

    Returns:
        Gaussian covariance matrix for either the real or imaginary part of
        the complex valley-coupling field.
    """
    if correlation_length_x <= 0.0:
        raise ValueError("correlation_length_x must be positive.")
    if correlation_length_y <= 0.0:
        raise ValueError("correlation_length_y must be positive.")

    dx = positions[:, None, 0] - positions[None, :, 0]
    dy = positions[:, None, 1] - positions[None, :, 1]
    r2 = dx * dx / correlation_length_x**2 + dy * dy / correlation_length_y**2
    return 0.5 * sigma_delta_complex**2 * np.exp(-r2)


def sample_valley_coupling_map(
    positions: np.ndarray,
    *,
    mean_valley_splitting: float,
    correlation_length_x: float,
    correlation_length_y: float,
    seed: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample a correlated valley-coupling field and its valley splitting.

    Args:
        positions: Array of ``(x, y)`` positions where the field is sampled.
        mean_valley_splitting: Target mean valley splitting ``<E_V>``.
        correlation_length_x: Correlation length along x.
        correlation_length_y: Correlation length along y.
        seed: Random seed used for reproducibility.

    Returns:
        Tuple containing the derived valley splitting ``E_V = 2*abs(Delta)``
        and the sampled complex intervalley coupling ``Delta``.
    """
    if mean_valley_splitting <= 0.0:
        raise ValueError("mean_valley_splitting must be positive.")

    # See Marcks et al., Nat. Commun. 16, 11381 (2025).
    sigma_delta_complex = mean_valley_splitting / np.sqrt(np.pi)
    # Get covariance matrix assuming Gaussian distribution.
    covariance = gaussian_covariance(
        positions,
        sigma_delta_complex=sigma_delta_complex,
        correlation_length_x=correlation_length_x,
        correlation_length_y=correlation_length_y,
    )
    # Sample multivariate Gaussian to account for correlation lengths.
    rng = np.random.default_rng(seed)
    mean = np.zeros(covariance.shape[0])
    delta_real = rng.multivariate_normal(mean=mean, cov=covariance)
    delta_imag = rng.multivariate_normal(mean=mean, cov=covariance)
    delta = delta_real + 1j * delta_imag
    # Compute valley splitting from valley coupling.
    valley_splitting = 2.0 * np.abs(delta)

    return valley_splitting, delta


def interpolate_map_to_nodes(
    x_coords: np.ndarray,
    y_coords: np.ndarray,
    map_values: np.ndarray,
    xy_nodes: np.ndarray,
) -> np.ndarray:
    """Interpolate a regular x-y map onto submesh nodes.

    Args:
        x_coords: x coordinates of the regular map.
        y_coords: y coordinates of the regular map.
        map_values: Values sampled on the regular ``(x, y)`` grid.
        xy_nodes: Mesh-node coordinates where the map should be evaluated.

    Returns:
        Interpolated map values on the mesh nodes.
    """
    interpolator = RegularGridInterpolator(
        (x_coords, y_coords),
        map_values,
        method="linear",
        bounds_error=False,
    )
    clipped_points = np.column_stack(
        (
            np.clip(xy_nodes[:, 0], x_coords[0], x_coords[-1]),
            np.clip(xy_nodes[:, 1], y_coords[0], y_coords[-1]),
        )
    )
    return interpolator(clipped_points)


def save_valley_map(
    file_name: str | Path,
    x_coords: np.ndarray,
    y_coords: np.ndarray,
    valley_splitting_grid: np.ndarray,
) -> None:
    """Save an image of the sampled valley-splitting map.

    Args:
        file_name: Output image path.
        x_coords: x coordinates of the regular map.
        y_coords: y coordinates of the regular map.
        valley_splitting_grid: Valley splitting sampled on the regular
            ``(x, y)`` grid.
    """
    figure, axis = plt.subplots(figsize=(5.5, 3.5))
    color_mesh = axis.pcolormesh(
        x_coords / 1e-9,
        y_coords / 1e-9,
        valley_splitting_grid.T / ct.e * 1e6,
        shading="auto",
    )
    axis.set_aspect("equal")
    axis.set_xlabel("x [nm]")
    axis.set_ylabel("y [nm]")

    divider = make_axes_locatable(axis)
    colorbar_axis = divider.append_axes("right", size="5%", pad=0.08)
    figure.colorbar(
        color_mesh,
        cax=colorbar_axis,
        label=r"Valley splitting [$\mu$eV]",
    )

    figure.tight_layout()
    figure.savefig(file_name, dpi=200)
    plt.close(figure)


def prefactor_x(me_inv: np.ndarray) -> np.ndarray:
    """Return the x kinetic-energy prefactor.

    Args:
        me_inv: Inverse electron effective-mass tensor on the mesh nodes.

    Returns:
        Prefactor multiplying the ``k_x^2`` operator.
    """
    return ct.hbar**2 / 2.0 * np.asarray(me_inv)[0, 0]


def prefactor_y(me_inv: np.ndarray) -> np.ndarray:
    """Return the y kinetic-energy prefactor.

    Args:
        me_inv: Inverse electron effective-mass tensor on the mesh nodes.

    Returns:
        Prefactor multiplying the ``k_y^2`` operator.
    """
    return ct.hbar**2 / 2.0 * np.asarray(me_inv)[1, 1]


def prefactor_z(me_inv: np.ndarray) -> np.ndarray:
    """Return the z kinetic-energy prefactor.

    Args:
        me_inv: Inverse electron effective-mass tensor on the mesh nodes.

    Returns:
        Prefactor multiplying the ``k_z^2`` operator.
    """
    return ct.hbar**2 / 2.0 * np.asarray(me_inv)[2, 2]


def build_two_valley_model(off_diagonal: np.ndarray) -> ElectronKPModel:
    """Build an explicit two-valley electron k·p model.

    Args:
        off_diagonal: Complex intervalley coupling ``Delta`` evaluated on
            the mesh nodes.

    Returns:
        Electron k·p model with anisotropic kinetic terms and local valley
        coupling.
    """
    delta_nodes = np.asarray(off_diagonal, dtype=np.complex128).reshape(-1)
    valley_matrix = np.zeros((delta_nodes.size, 2, 2), dtype=np.complex128)
    valley_matrix[:, 0, 1] = delta_nodes
    valley_matrix[:, 1, 0] = np.conjugate(delta_nodes)

    return ElectronKPModel(
        bands=2,
        constant=valley_matrix,
        quadratic={
            ("x", "x"): ElectronKPParameter("Me_inv", transform=prefactor_x),
            ("y", "y"): ElectronKPParameter("Me_inv", transform=prefactor_y),
            ("z", "z"): ElectronKPParameter("Me_inv", transform=prefactor_z),
        },
    )


def generate_correlated_two_valley_model(
    xy_nodes: np.ndarray,
    *,
    mean_valley_splitting: float,
    correlation_length_x: float,
    correlation_length_y: float,
    map_nx: int,
    map_ny: int,
    random_seed: int,
    valley_map_filename: str | Path | None = None,
) -> ElectronKPModel:
    """Build the two-valley k·p model.

    Args:
        xy_nodes: Mesh-node coordinates where the k·p model is evaluated.
        mean_valley_splitting: Mean value used to set the random-field scale.
        correlation_length_x: Correlation length along x.
        correlation_length_y: Correlation length along y.
        map_nx: Number of x-grid points in the sampled map.
        map_ny: Number of y-grid points in the sampled map.
        random_seed: Random seed used for reproducibility.
        valley_map_filename: Optional output path for a diagnostic plot of the
            sampled valley-splitting map.

    Returns:
        Electron k·p model with a correlated intervalley coupling evaluated on
        the mesh nodes.
    """
    x_coords = np.linspace(
        np.min(xy_nodes[:, 0]),
        np.max(xy_nodes[:, 0]),
        map_nx,
    )
    y_coords = np.linspace(
        np.min(xy_nodes[:, 1]),
        np.max(xy_nodes[:, 1]),
        map_ny,
    )
    x_grid, y_grid = np.meshgrid(x_coords, y_coords, indexing="ij")
    positions = np.stack((x_grid, y_grid), axis=-1).reshape((-1, 2))

    # Sample the valley coupling over the device.
    vs, delta = sample_valley_coupling_map(
        positions,
        mean_valley_splitting=mean_valley_splitting,
        correlation_length_x=correlation_length_x,
        correlation_length_y=correlation_length_y,
        seed=random_seed,
    )
    if valley_map_filename is not None:
        save_valley_map(
            file_name=valley_map_filename,
            x_coords=x_coords,
            y_coords=y_coords,
            valley_splitting_grid=vs.reshape(map_nx, map_ny),
        )

    # Interpolate onto mesh nodes.
    delta_mesh = interpolate_map_to_nodes(
        x_coords,
        y_coords,
        delta.reshape(map_nx, map_ny),
        xy_nodes,
    )

    # Create the k·p model.
    model = build_two_valley_model(delta_mesh)

    return model
