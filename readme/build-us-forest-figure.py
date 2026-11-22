"""Draw measured 3DEP point clouds and archived GEDI waveforms on common scales."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LightSource
from matplotlib.cm import ScalarMappable
from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "readme/data/us_forests"
FIG = ROOT / "readme"


def read(path):
    with path.open(newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def numeric(rows, field):
    return np.array([float(r[field]) for r in rows])


meta = read(DATA / "sites.csv")
cloud = read(DATA / "points.csv")
terrain = read(DATA / "terrain.csv")
bins = read(DATA / "bins.csv")
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11, "axes.labelcolor": "#24373C",
                     "text.color": "#24373C", "axes.spines.top": False, "axes.spines.right": False})
fig = plt.figure(figsize=(14.4, 20.5), facecolor="white")
outer = fig.add_gridspec(4, 1, height_ratios=[2.15, 1.2, 1.2, 1.2], hspace=.30,
                        left=.035, right=.98, top=.99, bottom=.105)
mapax = fig.add_subplot(outer[0])
mapax.imshow(Image.open(FIG / "fig-us-forest-sites.png"))
mapax.set_aspect("equal")
mapax.axis("off")
norm = Normalize(0, 60)
for i, site in enumerate(meta):
    row, side = divmod(i, 2)
    inner = outer[row + 1].subgridspec(1, 4, width_ratios=[2.65, 1.0, 2.65, 1.0], wspace=.12)
    ax = fig.add_subplot(inner[side * 2], projection="3d", computed_zorder=False)
    waveax = fig.add_subplot(inner[side * 2 + 1])
    points = [r for r in cloud if r["shot"] == site["shot"]]
    ground = [r for r in terrain if r["shot"] == site["shot"]]
    signal = [r for r in bins if r["shot"] == site["shot"]]
    x, y, z, height = [numeric(points, f) for f in
                        ["relative_easting_m", "relative_northing_m", "relative_elevation_m", "height_m"]]
    tx, ty, tz = [numeric(ground, f) for f in
                   ["relative_easting_m", "relative_northing_m", "relative_elevation_m"]]
    ux, uy = np.unique(tx), np.unique(ty)
    surface = np.full((len(uy), len(ux)), np.nan)
    surface[np.searchsorted(uy, ty), np.searchsorted(ux, tx)] = tz
    gx, gy = np.meshgrid(ux, uy)
    # Illumination is calculated from the measured DTM, without vertical exaggeration.
    hillshade = LightSource(azdeg=315, altdeg=40).hillshade(
        surface, vert_exag=1, dx=float(np.median(np.diff(ux))),
        dy=float(np.median(np.diff(uy))))
    facecolors = plt.cm.Greys_r(.22 + .65 * hillshade)
    ax.plot_surface(gx, gy, np.where(np.hypot(gx,gy)<=12.5,surface,np.nan), facecolors=facecolors, alpha=.50, zorder=1,
                    linewidth=0, rstride=1, cstride=1, shade=False)
    use = (height >= 0) & (np.hypot(x, y) <= 12.5)
    ax.scatter(x[use], y[use], z[use], c=height[use], cmap="viridis", norm=norm,
               s=3.0, alpha=.94, linewidths=0, rasterized=True, depthshade=False, zorder=5)
    ax.set(xlim=(-15, 15), ylim=(-15, 15), zlim=(-10, 60))
    ax.set_box_aspect((30, 30, 70))
    ax.set_proj_type("ortho")
    ax.view_init(elev=35.264, azim=135)
    ax.set_xticks([-10, 0, 10])
    ax.set_yticks([-10, 0, 10])
    ax.set_zticks([0, 20, 40, 60])
    ax.tick_params(labelsize=8.5, pad=-1)
    ax.set_xlabel("x (m)", fontsize=9, labelpad=-1)
    ax.set_ylabel("y (m)", fontsize=9, labelpad=-1)
    ax.set_zlabel("Relative elevation (m)", fontsize=9, labelpad=1)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.fill = False
        axis.pane.set_edgecolor("#E3E7E7")
        axis._axinfo["grid"]["color"] = "#E9EEEE"
        axis._axinfo["grid"]["linewidth"] = .4
    ax.set_title(f'{chr(66+i)}   {site["vegetation"]}', loc="left", fontweight="bold", fontsize=11.5, pad=14)
    ax.text2D(.03, 1.015, f'{site["landscape"]} · 3DEP {int(float(site["als_year"]))}',
              transform=ax.transAxes, fontsize=10)

    amplitude = numeric(signal, "amplitude_pct")
    elevation = numeric(signal, "height_m")
    waveax.fill_betweenx(elevation, 0, amplitude, color="#225C78", alpha=.17)
    waveax.plot(amplitude, elevation, color="#225C78", linewidth=.9)
    rh98 = float(site["rh98"])
    waveax.axhline(0, color="#5F625E", linewidth=1, linestyle=":")
    waveax.text(98, -3.2, "Ground 0 m", fontsize=9, color="#5F625E", ha="right")
    waveax.axhline(rh98, color="#A12B74", linestyle="--", linewidth=.8)
    waveax.text(98, rh98 + 1.2, f"RH98 {rh98:.1f} m", fontsize=9, color="#A12B74", ha="right")
    waveax.set(xlim=(0, 100), ylim=(-10, 60), xticks=[0, 50, 100], yticks=[0, 20, 40, 60])
    waveax.set_xlabel("Peak return (%)", fontsize=10)
    waveax.set_ylabel("Height above GEDI ground (m)", fontsize=9, labelpad=1)
    # Cumulative energy is a display integration, separate from official L2A RH.
    order = np.argsort(elevation)
    zz, aa = elevation[order], amplitude[order]
    selected = (zz >= -10) & (zz <= 60)
    zz, aa = zz[selected], aa[selected]
    energy = np.r_[0, np.cumsum(.5*(aa[1:]+aa[:-1])*np.diff(zz))]
    energy = energy / energy[-1]
    cumulative = waveax.twiny()
    cumulative.plot(energy, zz, color="#303C3B", linewidth=1)
    cumulative.set(xlim=(0,1), xticks=[0,.5,1])
    cumulative.set_xlabel("Cumulative energy", fontsize=9, labelpad=3)
    cumulative.tick_params(labelsize=8.5, pad=1)
    cumulative.spines['top'].set_visible(True)
    cumulative.spines['right'].set_visible(False)
    waveax.set_title(f'GEDI {site["date"]}', fontsize=9.5, pad=33)
    waveax.grid(axis="y", linewidth=.4, color="#E1E7E7")
    waveax.tick_params(labelsize=9)
# Separate 3D views retain the exact plotted measurements and camera.
    fig.canvas.draw()
    bounds = ax.get_tightbbox(fig.canvas.get_renderer()).transformed(fig.dpi_scale_trans.inverted())
    fig.savefig(DATA / f"cloud_{site['id']}.png", bbox_inches=bounds.expanded(1.04,1.08), dpi=300, facecolor="white")
cax = fig.add_axes([.30, .036, .35, .009])
bar = fig.colorbar(ScalarMappable(norm=norm, cmap="viridis"), cax=cax, orientation="horizontal")
bar.set_label("ALS return height above local ground (m)\nAll ALS clips have a 25 m diameter", fontsize=11)
bar.set_ticks([0, 10, 20, 30, 40, 50, 60])
bar.ax.tick_params(labelsize=9)
for suffix in ("png", "pdf"):
    fig.savefig(FIG / f"fig-us-forest-waveforms.{suffix}", dpi=400, facecolor="white")
plt.close(fig)
print("Wrote six ecosystem examples on shared horizontal and vertical scales.")




