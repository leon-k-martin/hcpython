def create_graph(f_sc, f_info, f_centers=None, thresh_percentile=0, output_file=None):
    import json
    from matplotlib.colors import to_hex
    import numpy as np
    import pandas as pd
    import numpy as np

    sc = pd.read_csv(
        f_sc,
        index_col=None,
        header=None,
    )

    info = pd.read_csv(
        f_info,
        sep="\s+",
        comment="#",
        header=None,
        names=["id", "label", "r", "g", "b", "a"],
    )
    if f_centers:
        centers = np.loadtxt(f_centers)
    else:
        centers = np.zeros((sc.shape[0], 3))

    thresh = np.percentile(sc, thresh_percentile)
    # graph = dict(directed=False, multigraph=False, graph={}, nodes=set(), links=set())
    nodes = list()
    links = list()
    for i in range(sc.shape[0]):
        i_info = {
            "id": info.loc[i, "id"].astype(float),
            "label": info.loc[i, "label"],
            "color": to_hex(
                np.array((info.loc[i, "r"], info.loc[i, "g"], info.loc[i, "b"])) / 255
            ),
            "x": centers[i, 0],
            "y": centers[i, 1],
            "z": centers[i, 2],
        }

        for j in range(i + 1, sc.shape[1]):
            val = sc.iloc[i, j].astype(float)
            j_info = {
                "id": info.loc[j, "id"].astype(float),
                "label": info.loc[j, "label"],
                "color": to_hex(
                    np.array((info.loc[j, "r"], info.loc[j, "g"], info.loc[j, "b"]))
                    / 255
                ),
                "x": centers[j, 0],
                "y": centers[j, 1],
                "z": centers[j, 2],
            }

            link = {
                "source": info.loc[i, "id"].astype(float),
                "target": info.loc[j, "id"].astype(float),
                "value": val,
            }
            if val > thresh:
                # if i_info not in graph["nodes"]:
                nodes.append(i_info)
                # if j_info not in graph["nodes"]:
                nodes.append(j_info)

                # if link not in graph["links"]:
                links.append(link)
            print(i, j)
    nodes = dict((node["id"], node) for node in nodes).values()  #
    links = dict(((link["source"], link["target"]), link) for link in links).values()
    graph = dict(
        directed=False,
        multigraph=False,
        graph={},
        nodes=list(nodes),
        links=list(links),
    )

    if output_file:
        with open(output_file, "w") as f:
            json.dump(graph, f, indent=4)
    return graph
