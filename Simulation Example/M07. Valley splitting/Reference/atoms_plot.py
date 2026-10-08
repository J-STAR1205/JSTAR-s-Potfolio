"""M07 Study 공통 모듈: 원자 배치도(Si/Ge 위치를 점으로 직접 그린 그림).

VESTA·ParaView 없이 파이썬만으로 원자 배열을 본다.
  1) 측면 단면도 (x–z): 얇은 y 슬라이스 안의 원자 — 우물, 계면 거칠기, 장벽의 Ge 배치
  2) 평면도 (x–y): 원자층 하나 — 계면 바로 위 층과 우물 중앙 층 비교
  3) 3D 산점도: 계면 부근 Ge 원자만
  4) 격자 확대도: 계면을 포함한 수 nm 창, 원자 + 결합선 (다이아몬드 격자 구조)

입력 좌표는 nm 단위 (N, 3) 배열, ge는 Ge 여부 bool 배열.
"""

import numpy as np
import matplotlib.pyplot as plt

SI_COLOR = "#7aa6d8"   # Si: 연한 파랑
GE_COLOR = "#d1453b"   # Ge: 빨강
LAYER_NM = 0.5431 / 4  # [001] 원자층 간격


def _dot_size(n_points, base=12.0):
    """점이 많을수록 작게."""
    return float(np.clip(base * 2000.0 / max(n_points, 1), 0.3, base))


def plot_side_view(pos, ge, well_top, well_bot, path, y0=0.0, slice_nm=0.3,
                   z_window=None, title=""):
    """x–z 단면 배치도. |y - y0| < slice_nm/2 인 원자만 그림."""
    sel = np.abs(pos[:, 1] - y0) < slice_nm / 2
    if z_window is not None:
        sel &= (pos[:, 2] >= z_window[0]) & (pos[:, 2] <= z_window[1])
    p, g = pos[sel], ge[sel]
    s = _dot_size(len(p))

    xr = np.ptp(p[:, 0]) if len(p) else 1.0
    zr = np.ptp(p[:, 2]) if len(p) else 1.0
    fig_w = 10.0
    fig_h = float(np.clip(fig_w * zr / max(xr, 1e-9), 3.0, 14.0))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    ax.scatter(p[~g, 0], p[~g, 2], s=s, c=SI_COLOR, lw=0, label=f"Si ({(~g).sum():,})")
    ax.scatter(p[g, 0], p[g, 2], s=s, c=GE_COLOR, lw=0, label=f"Ge ({g.sum():,})")
    for zz, name in ((well_top, "well top (mean)"), (well_bot, "well bottom (mean)")):
        ax.axhline(zz, color="k", ls="--", lw=0.8)
        ax.text(ax.get_xlim()[1] if len(p) else 0, zz, f" {name}", va="center", fontsize=8)
    ax.set_aspect("equal")
    ax.set_xlabel("x (nm)"); ax.set_ylabel("z (nm)")
    ax.set_title(f"Side view (x–z), slice |y − {y0:g}| < {slice_nm/2:g} nm  {title}")
    ax.legend(loc="upper right", markerscale=3, fontsize=8)
    fig.tight_layout(); fig.savefig(path, dpi=200)
    return fig


def _layer_mask(pos, z_target):
    """z_target에 가장 가까운 원자층 하나 (이완 후 흔들림 고려해 ±층간격/2)."""
    zc = pos[:, 2]
    # 실제 존재하는 층 높이 중 가장 가까운 것을 찾는다
    near = np.abs(zc - z_target) < 2 * LAYER_NM
    if not near.any():
        return np.zeros(len(zc), bool), z_target
    z_layer = zc[near][np.argmin(np.abs(zc[near] - z_target))]
    return np.abs(zc - z_layer) < LAYER_NM / 2, z_layer


def plot_top_views(pos, ge, layers, path, title=""):
    """평면 배치도. layers = [(라벨, z_nm), ...] 각 층을 한 칸씩."""
    n = len(layers)
    fig, axs = plt.subplots(1, n, figsize=(5.2 * n, 5.2), squeeze=False)
    for ax, (label, zt) in zip(axs[0], layers):
        m, z_layer = _layer_mask(pos, zt)
        p, g = pos[m], ge[m]
        s = _dot_size(len(p), base=18.0)
        ax.scatter(p[~g, 0], p[~g, 1], s=s, c=SI_COLOR, lw=0)
        ax.scatter(p[g, 0], p[g, 1], s=s, c=GE_COLOR, lw=0)
        frac = g.mean() * 100 if len(g) else float("nan")
        ax.set_aspect("equal")
        ax.set_xlabel("x (nm)"); ax.set_ylabel("y (nm)")
        ax.set_title(f"{label}\nz = {z_layer:.2f} nm, {len(p):,} atoms, Ge {frac:.1f} %",
                     fontsize=10)
    fig.suptitle(f"Top view of single atomic layers (Si blue, Ge red)  {title}")
    fig.tight_layout(rect=[0, 0, 1, 0.94]); fig.savefig(path, dpi=200)
    return fig


def plot_ge_3d(pos, ge, well_top, well_bot, path, z_window, max_points=20000,
               seed=0, title=""):
    """계면 부근 Ge 원자만 3D 산점도 (많으면 무작위로 솎아냄)."""
    m = ge & (pos[:, 2] >= z_window[0]) & (pos[:, 2] <= z_window[1])
    p = pos[m]
    if len(p) > max_points:
        p = p[np.random.default_rng(seed).choice(len(p), max_points, replace=False)]
    fig = plt.figure(figsize=(8, 7))
    ax = fig.add_subplot(projection="3d")
    ax.scatter(p[:, 0], p[:, 1], p[:, 2], s=2, c=p[:, 2], cmap="autumn", depthshade=False)
    xs = np.array([pos[:, 0].min(), pos[:, 0].max()])
    ys = np.array([pos[:, 1].min(), pos[:, 1].max()])
    X, Y = np.meshgrid(xs, ys)
    for zz in (well_top, well_bot):
        ax.plot_surface(X, Y, np.full_like(X, zz), alpha=0.15, color="tab:blue")
    ax.set_xlabel("x (nm)"); ax.set_ylabel("y (nm)"); ax.set_zlabel("z (nm)")
    ax.set_title(f"Ge atoms only, z ∈ [{z_window[0]:g}, {z_window[1]:g}] nm "
                 f"({len(p):,} shown)  {title}")
    fig.tight_layout(); fig.savefig(path, dpi=200)
    return fig


def plot_lattice_closeup(pos, ge, z_center, path, x0=None, width_nm=2.5,
                         height_nm=2.5, y0=0.0, depth_nm=0.5431, bond_nm=0.26,
                         title=""):
    """격자 확대도: 작은 창 안의 원자와 최근접 결합(< bond_nm)을 [010] 방향으로 투영.
    depth_nm = a(격자상수 한 칸)이면 다이아몬드 격자의 지그재그 결합이 보인다.
    Si–Si 0.235 nm, Ge–Ge 0.245 nm이므로 bond_nm = 0.26 nm면 최근접만 잡힌다."""
    if x0 is None:
        x0 = float(np.median(pos[:, 0])) - width_nm / 2
    m = ((pos[:, 0] >= x0) & (pos[:, 0] <= x0 + width_nm)
         & (np.abs(pos[:, 2] - z_center) <= height_nm / 2)
         & (pos[:, 1] >= y0) & (pos[:, 1] < y0 + depth_nm))
    p, g = pos[m], ge[m]
    fig, ax = plt.subplots(figsize=(7, 7 * height_nm / width_nm))
    if 1 < len(p) <= 5000:
        d = np.linalg.norm(p[:, None, :] - p[None, :, :], axis=-1)
        i, j = np.where(np.triu((d > 0.1) & (d < bond_nm)))
        for a_, b_ in zip(i, j):
            ax.plot(p[[a_, b_], 0], p[[a_, b_], 2], "-", color="0.55", lw=1.0, zorder=1)
    # 깊이(y)에 따라 점 크기를 달리해 앞뒤 층 구분
    depth = (p[:, 1] - y0) / depth_nm if len(p) else p[:, 1]
    size = 60 + 90 * (1 - depth)
    ax.scatter(p[~g, 0], p[~g, 2], s=size[~g], c=SI_COLOR, edgecolors="k",
               lw=0.4, zorder=2, label="Si")
    ax.scatter(p[g, 0], p[g, 2], s=size[g], c=GE_COLOR, edgecolors="k",
               lw=0.4, zorder=3, label="Ge")
    ax.axhline(z_center, color="k", ls="--", lw=0.8)
    ax.set_aspect("equal")
    ax.set_xlabel("x (nm)  [100]"); ax.set_ylabel("z (nm)  [001]")
    ax.set_title(f"Lattice close-up, projected along [010] (depth {depth_nm:g} nm)\n"
                 f"big = front, small = back; bonds < {bond_nm:g} nm  {title}", fontsize=10)
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout(); fig.savefig(path, dpi=200)
    return fig


def plot_all(pos, ge, well_top, well_bot, prefix, z_margin=3.0, slice_nm=0.3, title=""):
    """세 가지 배치도를 한 번에 저장. prefix는 경로 앞부분 (확장자 제외)."""
    zw = (well_bot - z_margin, well_top + z_margin)
    mid = 0.5 * (well_top + well_bot)
    plot_side_view(pos, ge, well_top, well_bot, f"{prefix}_side_full.png",
                   slice_nm=slice_nm, title=title)
    plot_side_view(pos, ge, well_top, well_bot, f"{prefix}_side_zoom.png",
                   slice_nm=slice_nm, z_window=zw, title=title)
    plot_top_views(pos, ge, [
        ("1 nm above well top (spacer)", well_top + 1.0),
        ("well top interface", well_top),
        ("well center", mid),
    ], f"{prefix}_top_layers.png", title=title)
    plot_ge_3d(pos, ge, well_top, well_bot, f"{prefix}_ge_3d.png", zw, title=title)
    plot_lattice_closeup(pos, ge, well_top, f"{prefix}_lattice_top_interface.png",
                         title=title)
