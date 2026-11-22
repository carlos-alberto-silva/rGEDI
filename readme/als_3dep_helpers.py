"""Read measured 3DEP point-cloud nodes around each displayed GEDI centre."""
import csv
import io
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import laspy
import numpy as np
import requests
from pyproj import Transformer

out = Path("04_data/fire_revision/forest_waveforms")
catalog = "https://s3-us-west-2.amazonaws.com/usgs-lidar-public"
projects = {
    "North Kaibab": "AZ_NorthKaibabNF_B1_2019",
    "Southern Appalachians": "USGS_LPC_TN_Eastern_2_16_B16_Del2_2016_LAS_2017",
}


def fetch(url):
    response = requests.get(url, timeout=(15, 45))
    response.raise_for_status()
    return response


def around_shot(root, lon, lat, cache=None):
    info = fetch(f"{root}/ept.json").json()
    bounds = info["bounds"]
    native = int(info["srs"]["horizontal"])
    x0, y0 = Transformer.from_crs(4326, native, always_xy=True).transform(lon, lat)
    utm = 32600 + int((lon + 180) // 6) + 1
    convert = Transformer.from_crs(native, utm, always_xy=True)
    east, north = convert.transform(x0, y0)
    halo = 45  # Native Pseudo-Mercator units, safely wider than the 24 m UTM search.

    def intersects(node):
        depth, x, y, _ = map(int, node.split("-"))
        unit = (bounds[3] - bounds[0]) / 2**depth
        nx, ny = bounds[0] + x * unit, bounds[1] + y * unit
        return nx <= x0 + halo and nx + unit >= x0 - halo and ny <= y0 + halo and ny + unit >= y0 - halo

    pages = ["0-0-0-0"]
    nodes = set()
    for page in pages:
        hierarchy = fetch(f"{root}/ept-hierarchy/{page}.json").json()
        for node, count in hierarchy.items():
            if not intersects(node):
                continue
            if count == -1:
                pages.append(node)
            elif count > 0:
                nodes.add(node)

    def read(node):
        path = f"{root}/ept-data/{node}.laz"
        if cache is None:
            cloud = laspy.read(io.BytesIO(fetch(path).content))
        else:
            cache.mkdir(parents=True, exist_ok=True)
            local = cache / f"{node}.laz"
            if not local.exists():
                content = fetch(path).content
                partial = local.with_suffix(".part")
                partial.write_bytes(content)
                partial.replace(local)
            cloud = laspy.read(local)
        x, y = np.asarray(cloud.x), np.asarray(cloud.y)
        near = (abs(x - x0) < halo) & (abs(y - y0) < halo)
        xx, yy = convert.transform(x[near], y[near])
        pts = np.column_stack((xx - east, yy - north, np.asarray(cloud.z)[near],
                               np.asarray(cloud.classification)[near]))
        return pts[(np.abs(pts[:, 0]) <= 24) & (np.abs(pts[:, 1]) <= 24)]

    if not nodes:
        raise RuntimeError(f"No indexed nodes at {lon}, {lat}")
    print(root, len(nodes), "intersecting index nodes", flush=True)
    with ThreadPoolExecutor(max_workers=8) as executor:
        return np.concatenate(list(executor.map(read, sorted(nodes))))


def ground_height(locations, ground):
    output = []
    for batch in np.array_split(locations, max(1, len(locations) // 150)):
        distance = ((batch[:, None, :2] - ground[None, :, :2]) ** 2).sum(axis=2)
        ix = np.argpartition(distance, 3, axis=1)[:, :4]
        weights = 1 / (np.take_along_axis(distance, ix, axis=1) + .01)
        output.append((weights * ground[ix, 2]).sum(axis=1) / weights.sum(axis=1))
    return np.concatenate(output)


def profile(points):
    ground = points[points[:, 3] == 2]
    sample = points[(np.hypot(points[:, 0], points[:, 1]) <= 12.5) &
                    np.isin(points[:, 3], [1, 2, 3, 4, 5])]
    if len(sample) < 150 or len(ground) < 20 or not (sample[:, 3] == 2).any():
        raise RuntimeError("Insufficient point and ground coverage at the GEDI centre")
    terrain_at_returns = ground_height(sample[:, :2], ground)
    origin = ground_height(np.zeros((1, 2)), ground)[0]
    xx, yy = np.meshgrid(np.arange(-14, 14.1, 1), np.arange(-14, 14.1, 1))
    grid = np.column_stack((xx.ravel(), yy.ravel()))
    terrain = np.column_stack((grid, ground_height(grid, ground) - origin))
    heights = sample[:, 2] - terrain_at_returns
    original = np.column_stack((sample[:, :2], sample[:, 2] - origin, heights, sample[:, 3]))
    sample[:, 2] = heights
    return sample, original, terrain


