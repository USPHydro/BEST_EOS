import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator


PATHS = {
    "Critical EOS": "/storage/home/vfrancao/BEST_EOS/Files_PAR_143_350_3_93_143_286",
    "LAT-only HRG": "/storage/home/vfrancao/BEST_EOS/LATonly_smooth_HRG",
}

MU_VALUES = [0.20, 0.25, 0.30, 0.35, 0.40, 0.45]  # GeV
T_MIN, T_MAX = 0.05, 0.40  # GeV
MU_TOL = 0.001  # GeV; inverted points shown within 1 MeV of target mu_B


def read_original(base_path):
    data = np.loadtxt(f"{base_path}/CorrLength.dat")
    mu = data[:, 0] / 1000.0  # MeV -> GeV
    temperature = data[:, 1] / 1000.0  # MeV -> GeV
    xi = data[:, 2]  # fm
    return mu, temperature, xi


def read_inverted(base_path):
    tables = [
        np.loadtxt(f"{base_path}/visuEOS_{table}.dat")
        for table in range(1, 8)
    ]
    data = np.vstack(tables)

    # visuEOS columns: T, muB, e, nB, s, p, cs2, xi, valid.
    temperature = data[:, 0]  # GeV
    mu = data[:, 1]  # GeV
    xi = data[:, 7]  # fm
    valid = data[:, 8] > 0.5
    return mu, temperature, xi, valid


def original_interpolator(mu, temperature, xi):
    mu_grid = np.unique(mu)
    temperature_grid = np.unique(temperature)

    # CorrLength.dat varies mu_B first for each fixed T.
    xi_grid = xi.reshape(temperature_grid.size, mu_grid.size).T
    return RegularGridInterpolator(
        (mu_grid, temperature_grid),
        xi_grid,
        method="linear",
        bounds_error=False,
        fill_value=np.nan,
    )


fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharex=True)

for ax, (title, base_path) in zip(axes, PATHS.items()):
    mu_original, T_original, xi_original = read_original(base_path)
    mu_inverted, T_inverted, xi_inverted, valid = read_inverted(base_path)
    interpolate_original = original_interpolator(
        mu_original, T_original, xi_original
    )

    comparison_mask = (
        valid
        & (T_inverted >= T_MIN)
        & (T_inverted <= T_MAX)
    )
    original_at_inverted = interpolate_original(
        np.column_stack((mu_inverted[comparison_mask], T_inverted[comparison_mask]))
    )
    absolute_error = np.abs(xi_inverted[comparison_mask] - original_at_inverted)
    finite_error = absolute_error[np.isfinite(absolute_error)]

    print(title)
    print(f"  compared valid points: {finite_error.size}")
    print(f"  max |xi_inverted - xi_original|: {finite_error.max():.6g} fm")
    print(f"  mean absolute error: {finite_error.mean():.6g} fm")

    for index, mu_value in enumerate(MU_VALUES):
        color = f"C{index}"

        original_mask = np.isclose(mu_original, mu_value, atol=1.0e-12)
        order_original = np.argsort(T_original[original_mask])
        ax.plot(
            T_original[original_mask][order_original],
            xi_original[original_mask][order_original],
            color=color,
            linewidth=2,
            label=rf"$\mu_B={mu_value:.2f}$ GeV",
        )

        inverted_mask = (
            valid
            & (np.abs(mu_inverted - mu_value) <= MU_TOL)
            & (T_inverted >= T_MIN)
            & (T_inverted <= T_MAX)
        )
        order_inverted = np.argsort(T_inverted[inverted_mask])
        ax.scatter(
            T_inverted[inverted_mask][order_inverted],
            xi_inverted[inverted_mask][order_inverted],
            color=color,
            marker="o",
            s=10,
            alpha=0.45,
        )

    ax.set_title(title, fontsize=18)
    ax.set_xlabel(r"$T$ (GeV)", fontsize=20)
    ax.set_xlim(T_MIN, T_MAX)
    ax.tick_params(axis="both", which="major", labelsize=14)
    ax.grid(alpha=0.15)

axes[0].set_ylabel(r"$\xi$ (fm)", fontsize=20)
axes[0].legend(fontsize=11, frameon=False, ncol=2)

# Solid curves: original regular (mu_B,T) table.
# Points: inverted (e,n_B) table, selected within MU_TOL of each target mu_B.
fig.text(
    0.5,
    0.01,
    r"Linhas: CorrLength.dat; pontos: visuEOS ($|\Delta\mu_B|\leq 1$ MeV)",
    ha="center",
    fontsize=13,
)
fig.tight_layout(rect=(0, 0.05, 1, 1))
fig.savefig("xi_original_vs_inverted.pdf", bbox_inches="tight")
fig.savefig("xi_original_vs_inverted.png", dpi=300, bbox_inches="tight")
plt.show()
