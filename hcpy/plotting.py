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
    import plotly.graph_objects as go
    import numpy as np
    import hcp_utils as hcp

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
