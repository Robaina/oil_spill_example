from __future__ import annotations
import networkx as nx
import json
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt


def generate_bipartite_graph(
    exchanges_file_path: str,
    hide_taxa: list[str] = None,
    hide_metabolites: list[str] = None,
    flux_cutoff: float = None,
) -> nx.DiGraph:
    """
    Generates a bipartite graph from a file containing exchange data.

    Args:
        exchanges_file_path (str): The path to the file containing exchange data.
        hide_taxa (list[str], optional): A list of taxa to hide from the graph. Defaults to None.
        hide_metabolites (list[str], optional): A list of metabolites to hide from the graph. Defaults to None.
        flux_cutoff (float, optional): The minimum flux value for an edge to be included in the graph. Defaults to None.

    Returns:
        nx.DiGraph: The generated bipartite graph.
    """
    G = nx.DiGraph()

    with open(exchanges_file_path, "r") as f:
        for line in f:
            if "taxon" in line:
                continue
            cols = line.strip().split("\t")
            taxon = cols[1].replace("_sp", "")
            metabolite = cols[7].replace("_e", "")
            direction = cols[8]
            flux = abs(float(cols[5]))
            if hide_taxa is None or taxon not in hide_taxa:
                G.add_node(taxon, bipartite=0)
            if hide_metabolites is None or metabolite not in hide_metabolites:
                G.add_node(metabolite, bipartite=1)
            if flux_cutoff is None or flux >= flux_cutoff:
                if direction == "export":
                    G.add_edge(taxon, metabolite)
                elif direction == "import":
                    G.add_edge(metabolite, taxon)

    if hide_taxa is not None:
        G.remove_nodes_from(hide_taxa)
    if hide_metabolites is not None:
        G.remove_nodes_from(hide_metabolites)

    G.remove_nodes_from(list(nx.isolates(G)))

    bipartite_graph = G
    data = nx.readwrite.json_graph.node_link_data(G)

    with open("results/micom/graph.json", "w") as f:
        json.dump(data, f)

    return bipartite_graph


def plot_interaction_graph(
    final_df: pd.DataFrame,
    bipartite_graph: nx.Graph,
    target_taxon: str = "Acinetobacter",
    target_compound: str = "tol",
    color_carbon_sources: str = "violet",
    color_acinetobacter: str = "skyblue",
    color_other_edges: str = "lightgrey",
    color_other_nodes: str = "silver",
    large_node_size: int = 6000,
    small_node_size: int = 1000,
    edge_width_acinetobacter: float = 3.5,
    edge_width_other: float = 0.8,
    arrow_size_acinetobacter: int = 10,
    arrow_size_other: int = 5,
) -> None:
    """
    This function plots an interaction graph based on the provided bipartite graph and dataframe.

    Parameters:
    - final_df: A pandas DataFrame containing the final data.
    - bipartite_graph: A NetworkX bipartite graph.
    - target_taxon: The target taxon for the graph.
    - color_carbon_sources: The color for carbon sources in the graph.
    - color_acinetobacter: The color for Acinetobacter in the graph.
    - color_other_edges: The color for other edges in the graph.
    - color_other_nodes: The color for other nodes in the graph.
    - large_node_size: The size for large nodes in the graph.
    - small_node_size: The size for small nodes in the graph.
    - edge_width_acinetobacter: The width for Acinetobacter edges in the graph.
    - edge_width_other: The width for other edges in the graph.
    - arrow_size_acinetobacter: The size for Acinetobacter arrows in the graph.
    - arrow_size_other: The size for other arrows in the graph.
    - hidden_metabolites: A list of metabolites to be hidden in the graph.
    """
    medium_donors = [
        m.replace("_m", "")
        for m in final_df[final_df["donor"] == "medium"]["compound"].unique()
    ]
    direct_nodes = list(bipartite_graph.neighbors(target_compound)) + [target_compound]
    indirect_nodes = [
        neighbor
        for node in direct_nodes
        for neighbor in bipartite_graph.neighbors(node)
    ]
    indirect_nodes_extended = [
        neighbor
        for node in indirect_nodes
        for neighbor in bipartite_graph.neighbors(node)
    ]
    acinetobacter_compounds = [
        n
        for n in bipartite_graph.to_undirected().neighbors(target_taxon)
        if bipartite_graph.nodes[n]["bipartite"] == 1
    ]
    acinetobacter_consumed_compounds = [
        u for u, v in bipartite_graph.edges() if v == target_taxon
    ]
    species_connected_to_acinetobacter_compounds = [
        n
        for compound in acinetobacter_compounds
        for n in bipartite_graph.neighbors(compound)
        if bipartite_graph.nodes[n]["bipartite"] == 0
    ]
    subgraph_nodes = (
        direct_nodes
        + indirect_nodes
        + indirect_nodes_extended
        + acinetobacter_compounds
        + species_connected_to_acinetobacter_compounds
    )
    subgraph = bipartite_graph.subgraph(subgraph_nodes)
    medium_donor_edges = [
        (u, v)
        for u, v in bipartite_graph.edges()
        if u in medium_donors and bipartite_graph.nodes[v]["bipartite"] == 0
    ]
    extended_subgraph = nx.DiGraph(subgraph)
    extended_subgraph.add_edges_from(medium_donor_edges)

    for node in extended_subgraph.nodes():
        if "bipartite" not in extended_subgraph.nodes[node]:
            extended_subgraph.nodes[node]["bipartite"] = 1

    shell_layout_extended_subgraph = nx.shell_layout(
        extended_subgraph,
        [
            set(
                n for n, d in extended_subgraph.nodes(data=True) if d["bipartite"] == 0
            ),
            set(
                n for n, d in extended_subgraph.nodes(data=True) if d["bipartite"] == 1
            ),
        ],
    )

    plt.figure(figsize=(10, 10))

    edge_colors, edge_widths, arrow_sizes = [], [], []
    for u, v in extended_subgraph.edges():
        if v == target_taxon:
            if u in medium_donors:
                edge_colors.append(color_carbon_sources)
                edge_widths.append(edge_width_acinetobacter)
                arrow_sizes.append(arrow_size_acinetobacter)
            else:
                edge_colors.append(color_acinetobacter)
                edge_widths.append(edge_width_acinetobacter)
                arrow_sizes.append(arrow_size_acinetobacter)
        elif u in medium_donors:
            edge_colors.append(color_carbon_sources)
            edge_widths.append(edge_width_other)
            arrow_sizes.append(arrow_size_other)
        else:
            edge_colors.append(color_other_edges)
            edge_widths.append(edge_width_other)
            arrow_sizes.append(arrow_size_other)

    nx.draw(
        extended_subgraph,
        shell_layout_extended_subgraph,
        with_labels=True,
        node_size=[
            large_node_size if bipartite == 0 else small_node_size
            for bipartite in nx.get_node_attributes(
                extended_subgraph, "bipartite"
            ).values()
        ],
        node_color=[
            color_acinetobacter
            if (
                node == target_taxon
                or (
                    node in acinetobacter_consumed_compounds
                    and node not in medium_donors
                )
            )
            else (
                color_other_nodes
                if bipartite == 0
                else color_carbon_sources
                if node in medium_donors
                else color_other_nodes
            )
            for node, bipartite in nx.get_node_attributes(
                extended_subgraph, "bipartite"
            ).items()
        ],
        edge_color=edge_colors,
        arrowsize=arrow_sizes,
        width=edge_widths,
    )
    nx.draw_networkx_edge_labels(
        extended_subgraph,
        shell_layout_extended_subgraph,
        edge_labels=nx.get_edge_attributes(extended_subgraph, "weight"),
    )
    plt.show()
