from __future__ import annotations
from pathlib import Path
import json
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from pandas import DataFrame


def generate_bipartite_graph(
    exchanges_file_path: str,
    hide_taxa: list[str] = None,
    hide_metabolites: list[str] = None,
    flux_cutoff: float = None,
    target_taxon: str = None,
    environmental_carbon_sources: list[str] = None,
    output_graph: Path = None,
    relabel_nodes: dict = None,
    keep_metabolites: list[str] = None,  # New parameter
) -> nx.DiGraph:
    """
    Generates a bipartite graph from a file containing exchange data.

    Args:
        exchanges_file_path (str): The path to the file containing exchange data.
        hide_taxa (list[str], optional): A list of taxa to hide from the graph. Defaults to None.
        hide_metabolites (list[str], optional): A list of metabolites to hide from the graph. Defaults to None.
        flux_cutoff (float, optional): The minimum flux value for an edge to be included in the graph. Defaults to None.
        target_taxon (str, optional): The target taxon to focus on in the graph. Defaults to None.
        environmental_carbon_sources (list[str], optional): A list of environmental carbon sources. Defaults to None.
        keep_metabolites (list[str], optional): A list of metabolites to keep in the graph, even if they are not connected to the target taxon. Defaults to None.

    Returns:
        nx.DiGraph: The generated bipartite graph.
    """
    G = nx.DiGraph()
    if relabel_nodes is not None:
        G = nx.relabel_nodes(G, relabel_nodes)

    with open(exchanges_file_path, "r") as f:
        for line in f:
            if "taxon" in line:
                continue
            cols = line.strip().split("\t")
            taxon = cols[1].replace("_sp", "")
            metabolite = cols[7].replace("_e", "")
            if relabel_nodes is not None and metabolite in relabel_nodes:
                metabolite = relabel_nodes[metabolite]
            direction = cols[8]
            flux = abs(float(cols[5]))
            if hide_taxa is None or taxon not in hide_taxa:
                G.add_node(taxon, bipartite=0)
            if hide_metabolites is None or metabolite not in hide_metabolites:
                if (
                    environmental_carbon_sources is None
                    or metabolite in environmental_carbon_sources
                ):
                    G.add_node(metabolite, bipartite=1)
                elif target_taxon is not None and (
                    taxon == target_taxon or metabolite == target_taxon
                ):
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

    for source in environmental_carbon_sources:
        if source not in G.nodes:
            G.add_node(source, bipartite=1)

    for node in list(G.nodes):
        if (
            "bipartite" in G.nodes[node]
            and (node not in environmental_carbon_sources)
            and G.nodes[node]["bipartite"] == 1
        ):
            if not G.has_edge(node, target_taxon) and (
                keep_metabolites is None or node not in keep_metabolites
            ):
                G.remove_node(node)

    G.remove_nodes_from(list(nx.isolates(G)))

    bipartite_graph = G
    data = nx.readwrite.json_graph.node_link_data(G)

    if output_graph is not None:
        with open(output_graph, "w") as f:
            json.dump(data, f)

    return bipartite_graph


def plot_ecological_interactions(
    ko: DataFrame,
    figsize: tuple = (15, 15),
    node_size: int = 10000,
    node_color: str = "grey",
    arrowsize: int = 20,
) -> None:
    """
    Function to plot an interaction graph based on a knockout analysis dataframe.

    Parameters:
    ko (DataFrame): A pandas DataFrame from a knockout analysis.
    figsize (tuple): A tuple specifying the size of the figure. Default is (15, 15).
    node_size (int): The size of the nodes in the graph. Default is 10000.
    node_color (str): The color of the nodes in the graph. Default is 'grey'.
    arrowsize (int): The size of the arrow heads. Default is 20.

    Returns:
    None
    """
    G = nx.DiGraph()
    for taxon in ko.columns:
        G.add_node(taxon)
    for taxon1 in ko.columns:
        for taxon2 in ko.columns:
            if taxon1 != taxon2:
                growth_change = ko.loc[taxon1, taxon2]
                if growth_change > 0:
                    G.add_edge(taxon1, taxon2, color="red")
                elif growth_change < 0:
                    G.add_edge(taxon1, taxon2, color="blue")
    colors = nx.get_edge_attributes(G, "color").values()
    pos = nx.circular_layout(G)
    plt.figure(figsize=figsize)
    node_labels = {node: node.replace("_sp", "") for node in G.nodes()}
    nx.draw(
        G,
        pos,
        labels=node_labels,
        font_color="white",
        node_size=node_size,
        node_color=node_color,
        edge_color=colors,
        arrowsize=arrowsize,
    )
    plt.show()


def plot_trophic_interactions(
    bipartite_graph: nx.Graph,
    environmental_carbon_sources: list[str],
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
    medium_donors = environmental_carbon_sources
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
