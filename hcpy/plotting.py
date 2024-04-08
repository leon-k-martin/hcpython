import os
import sys

import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
from matplotlib.colors import ListedColormap
from nilearn import plotting

from hcpy import data


# plotting.py
def plot_bvecs(bvecs, bvals=None, fout="bvecs.html"):
    import os

    import numpy as np
    import plotly.graph_objects as go

    if os.path.exists(bvecs):
        bvecs = np.loadtxt(bvecs)
    if os.path.exists(bvals):
        bvals = np.loadtxt(bvals)
    # Load your b-vectors
    if not isinstance(bvals, type(None)):
        color = bvals
    else:
        color = "red"

    # Define initial figure
    fig = go.Figure(
        data=[
            go.Scatter3d(
                x=bvecs[0, :],
                y=bvecs[1, :],
                z=bvecs[2, :],
                mode="markers",
                marker=dict(size=5, color=color),
            )
        ]
    )
    # Compute limits for each axis
    x_min, x_max = bvecs[0, :].min(), bvecs[0, :].max()
    y_min, y_max = bvecs[1, :].min(), bvecs[1, :].max()
    z_min, z_max = bvecs[2, :].min(), bvecs[2, :].max()

    # Update figure layout to zoom in to the limits of the points and remove all axes
    fig.update_layout(
        title=dict(text="B-vectors", x=0.5, xanchor="center", font=dict(size=20)),
        scene=dict(
            xaxis=dict(visible=False, range=[x_min, x_max]),
            yaxis=dict(visible=False, range=[y_min, y_max]),
            zaxis=dict(visible=False, range=[z_min, z_max]),
            xaxis_showbackground=False,
            yaxis_showbackground=False,
            zaxis_showbackground=False,
        ),
    )
    fig.write_html(fout)


def plot_brain_data():
    import hcp_utils as hcp
    import numpy as np
    import plotly.graph_objects as go

    vertices, triangles = hcp.mesh.inflated_left

    # Extracting x, y, z coordinates for plotting
    x, y, z = vertices.T

    # Defining the mesh
    mesh = go.Mesh3d(
        x=x,
        y=y,
        z=z,
        i=triangles[:, 0],
        j=triangles[:, 1],
        k=triangles[:, 2],
        opacity=1,
        color="red",
    )

    # Creating the figure and adding the mesh
    fig = go.Figure(data=[mesh])
    fig.update_layout(
        scene=dict(
            xaxis=dict(showticklabels=False, showbackground=False),
            yaxis=dict(showticklabels=False, showbackground=False),
            zaxis=dict(showticklabels=False, showbackground=False),
            camera=dict(eye=dict(x=-2.5, y=0, z=0)),  # Adjusting the camera's position
        ),
        width=700,
        margin=dict(r=20, l=10, b=10, t=10),
    )

    # Showing the figure
    fig.show()


def plot_parc2subject(
    subid, parc_name, t1_image_path, parcellation_image_path, output_png_path, roi_info
):

    # Load T1 image and parcellation image
    t1_img = nib.load(t1_image_path)
    parcellation_img = nib.load(parcellation_image_path)

    # Create a custom colormap using the listed colors
    cmap = ListedColormap(
        roi_info.loc[:, ["R", "G", "B", "A"]].values.astype(float) / 255
    )

    disp = plotting.plot_roi(
        parcellation_img,
        t1_img,
        title=f"{subid} {parc_name}",
        display_mode="mosaic",
        cmap=cmap,
        alpha=0.8,
    )

    # Save the overlay image as a PNG file
    disp.savefig(output_png_path)
    plt.close()


def plot_mmp2subject(subid, t1_image_path, parcellation_image_path, output_png_path):

    plot_parc2subject(
        subid,
        "MMP",
        t1_image_path,
        parcellation_image_path,
        output_png_path,
        data.mmp_atlas_info,
    )


def plot_dk2subject(subid, t1_image_path, parcellation_image_path, output_png_path):

    plot_parc2subject(
        subid,
        "DK",
        t1_image_path,
        parcellation_image_path,
        output_png_path,
        data.dk_atlas_info,
    )


def plot_sc(subid, weights, lengths, output_fig):
    import pandas as pd
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    import numpy as np

    fig, axs = plt.subplots(ncols=2, figsize=(16, 9))
    fig.suptitle = f"SC {subid}"
    weights = pd.read_csv(weights, index_col=0)
    lengths = pd.read_csv(lengths, index_col=0)
    im1 = axs[0].imshow(
        weights,
        cmap="viridis",
        aspect="equal",
        interpolation="none",
        # norm=LogNorm()
    )
    axs[0].set_title("Streamline count")
    fig.colorbar(im1, ax=axs[0], shrink=0.5)
    im2 = axs[1].imshow(
        lengths,
        cmap="viridis",
        aspect="equal",
        interpolation="none",
        # norm=LogNorm()
    )
    axs[1].set_title("Streamline length")
    fig.colorbar(im2, ax=axs[1], shrink=0.5)
    fig.savefig(output_fig)


import nibabel as nib

fod = nib.load(
    "/Volumes/bronkodata_work/hcpython/testsub/HCD0001305_V1_MR/T1w/Diffusion/MRTrix/wmfod_norm.nii.gz"
)


def plot_fod_histogram(fod, fod_norm, output_fig=None, ax=None):

    if isinstance(ax, type(None)):
        fig = plt.figure(layout="constrained", figsize=(6, 6))
        ax = fig.add_subplot(111, aspect="auto")
        return_fig = True
    else:
        return_fig = False

    print(ax, type(ax))
    if isinstance(fod, str):
        fod = nib.load(fod)
    if isinstance(fod_norm, str):
        fod_norm = nib.load(fod_norm)

    fod_norm_data = np.abs(fod_norm.get_fdata()).flatten()
    fod_data = np.abs(fod.get_fdata()).flatten()
    # Plot a histogram of the absolute values to visualize the distribution
    ax.hist(fod_data, bins=1000, log=True, alpha=1, label="FOD", color="black")
    ax.hist(
        fod_norm_data, bins=1000, log=True, alpha=0.6, label="Normalized FOD", color="r"
    )

    x1, x2, y1, y2 = -1.5, -0.9, -2.5, -1.9  # subregion of the original image
    axins = ax.inset_axes([0.5, 0.5, 0.47, 0.47], xlim=(0, 0.02))
    axins.hist(fod_data, bins=1000, log=True, alpha=0.8, label="FOD", color="black")
    axins.hist(
        fod_norm_data, bins=1000, log=True, alpha=0.7, label="Normalized FOD", color="r"
    )

    axins.legend()
    ax.indicate_inset_zoom(axins, edgecolor="k")

    plt.title("Histogram of FOD Absolute Values")
    plt.xlabel("Absolute Value")
    plt.ylabel("Log Frequency")

    if output_fig:
        plt.savefig(output_fig)

    if return_fig:
        plt.close()
        return fig


def plot_sift_weights(sift_weights, ax=None, output_fig=None):
    if ax is None:
        fig, ax = plt.subplots()
        return_fig = True
    else:
        return_fig = False

    ax.hist(
        sift_weights,
        bins=1000,
        log=True,
        alpha=0.9,
        label="SIFT Weights",
        color="black",
    )
    ax.set_title("Histogram of SIFT Weights")

    if output_fig:
        plt.savefig(output_fig)
    if return_fig:
        plt.close()
        return fig


def plot_tck_metrics(sift_weights, fod, fod_norm, output_fig=None):
    fig, axs = plt.subplots(1, 2, layout="constrained", figsize=(10, 5))
    plotting.plot_sift_weights(sift_weights, ax=axs[0])
    plotting.plot_fod_histogram(fod, fod_norm, ax=axs[1])

    if output_fig:
        fig.savefig(output_fig)


import numpy as np
from dipy.io.streamline import load_tck
from dipy.io.image import load_nifti
from dipy.viz import window, actor


def scene_plot_tractogram(tckfile, dwi_file, t1_file, slice_opacity=1, zoom=1.9):
    tckfile = load_tck(tckfile, dwi_file)
    data, affine = load_nifti(t1_file)

    shape = data.shape

    # Create a visualization scene
    scene = window.Scene(background="white")

    # Add tractography streamlines to the scene
    stream_actor = actor.line(
        tckfile.streamlines,
        linewidth=0.1,
        opacity=0.8,
        fake_tube=True,
    )
    scene.add(stream_actor)
    # Scale the streamlines actor, if applicable

    # Slicer opacity
    slicer_opacity = slice_opacity

    # Create and configure the sagittal slice (lateral view)
    image_actor_x = actor.slicer(data, affine)
    image_actor_x.opacity(slicer_opacity)
    x_midpoint = int(np.round(shape[0] / 2))
    image_actor_x.display(x_midpoint)

    # Create and configure the coronal slice (frontal view)
    image_actor_y = actor.slicer(data, affine)
    image_actor_y.opacity(slicer_opacity)
    y_midpoint = int(np.round(shape[1] / 2))
    image_actor_y.display(y=y_midpoint)

    # Create and configure the axial slice (superior view)
    image_actor_z = actor.slicer(data, affine)
    image_actor_z.opacity(slicer_opacity)
    z_midpoint = int(np.round(shape[2] / 2))
    image_actor_z.display_extent(
        0, shape[0] - 1, 0, shape[1] - 1, z_midpoint, z_midpoint
    )

    # Add slices to the scene
    scene.add(image_actor_z)
    scene.add(image_actor_x)
    scene.add(image_actor_y)

    scene.set_camera(
        position=(637.38, -130.06, 18.40),
        focal_point=(-0.40, -17.60, 18.40),
        view_up=(0, 0, 1),
    )
    scene.zoom(zoom)

    return scene


def scene2mplarray(scene):
    im = window.snapshot(scene, fname=None, size=(800, 800))
    return im


def plot_tractogram(tckfile, dwi_file, t1_file, out_path, zoom=1.9):
    scene = scene_plot_tractogram(
        tckfile, dwi_file, t1_file, slice_opacity=1, zoom=zoom
    )

    window.record(
        scene,
        out_path=out_path,
        size=(600, 600),
        reset_camera=False,
        n_frames=4,
        az_ang=91,
        path_numbering=True,
        magnification=2,
    )


import os
import imageio
from os.path import dirname, join


def rotating_tractogram(
    tckfile, dwi_file, t1_file, output_gif, zoom=1.9, opacity=0, duration=0.1
):
    scene = scene_plot_tractogram(
        tckfile, dwi_file, t1_file, slice_opacity=opacity, zoom=zoom
    )

    n_frames = 36  # Number of frames for one complete turn
    az_ang = 360 / n_frames  # Azimuthal angle for each frame to complete a 360° turn
    png_dir = join(dirname(output_gif), "gif")
    os.makedirs(png_dir, exist_ok=True)
    window.record(
        scene,
        out_path=f"{png_dir}/gif_tractogram",
        size=(600, 600),
        reset_camera=False,
        n_frames=n_frames,
        az_ang=az_ang,
        path_numbering=True,
        magnification=2,
    )
    del scene

    png_files = sorted(
        [f for f in os.listdir(png_dir) if f.endswith(".png") and f.startswith("gif")]
    )
    # Read each PNG image and append it to the list of images
    images = []
    for filename in png_files:
        image_path = os.path.join(png_dir, filename)
        images.append(imageio.imread(image_path))

    if not output_gif.endswith(".gif"):
        output_gif += ".gif"
    # Save the images as a GIF
    imageio.mimsave(
        output_gif, images, duration=duration, loop=0
    )  # Adjust 'duration' as needed
    os.system(f"rm -rf {png_dir}")
