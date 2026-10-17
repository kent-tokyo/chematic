from pathlib import Path
import io
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import trimesh
import chematic as ch

out = Path("validation/results/crystal-obj-mesh-2026-10-17")
na = [(0, 0, 0), (0, 0.5, 0.5), (0.5, 0, 0.5), (0.5, 0.5, 0)]
cl = [(0.5, 0, 0), (0, 0.5, 0), (0, 0, 0.5), (0.5, 0.5, 0.5)]
source = ch.PeriodicStructure(
    ch.Lattice.cubic(5.6402),
    [ch.Site([("Na", 1.0)], x) for x in na] + [ch.Site([("Cl", 1.0)], x) for x in cl],
)
cif = source.to_cif()
(out / "nacl.cif").write_text(cif)
source = ch.PeriodicStructure.from_cif(cif)
result = ch.write_crystal_obj_checked(source, [0.65] * 4 + [0.9] * 4, subdivisions=2)
(out / "nacl.obj").write_text(result.text)
scene = trimesh.load(
    io.StringIO(result.text),
    file_type="obj",
    force="scene",
    process=False,
    group_material=False,
    skip_materials=True,
)
fig = plt.figure(figsize=(7, 6), facecolor="white")
ax = fig.add_subplot(projection="3d")
light = np.array([1.0, -1.0, 2.0])
light /= np.linalg.norm(light)
for name, mesh in scene.geometry.items():
    base = (
        np.array([0.18, 0.45, 0.85])
        if name.endswith("_Na")
        else np.array([0.25, 0.7, 0.35])
    )
    shade = 0.5 + 0.5 * np.maximum(0, mesh.face_normals @ light)
    faces = Poly3DCollection(
        mesh.triangles,
        facecolors=shade[:, None] * base,
        edgecolors="none",
        zsort="average",
    )
    ax.add_collection3d(faces)
ax.set_xlim(-1, 4)
ax.set_ylim(-1, 4)
ax.set_zlim(-1, 4)
ax.set_box_aspect((1, 1, 1))
ax.view_init(elev=23, azim=35)
ax.set_xlabel("x (Å)")
ax.set_ylabel("y (Å)")
ax.set_zlabel("z (Å)")
ax.set_title("NaCl atom mesh · 8 spheres / 1,024 triangles", pad=20)
fig.text(
    0.5,
    0.035,
    "OBJ read back with trimesh; preview colors added by the renderer.",
    ha="center",
    fontsize=9,
    color="#444444",
)
fig.subplots_adjust(left=0, right=0.95, bottom=0.10, top=0.9)
fig.savefig(out / "nacl-preview.png", dpi=150)
