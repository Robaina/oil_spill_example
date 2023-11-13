from __future__ import annotations
from pathlib import Path
import json
import random
import networkx as nx
import matplotlib.pyplot as plt


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
    ko,
    figsize: tuple = (10, 10),
    selected_taxa: list[str] = None,
    node_size: int = 10000,
    node_color: str = "silver",
    edge_width: int = 3,
    cooperative_color: str = "lightblue",
    competitive_color: str = "lightcoral",
    arrowsize: int = 20,
) -> None:
    """
    Function to plot an interaction graph based on a knockout analysis dataframe.

    Parameters:
    ko (DataFrame): A pandas DataFrame from a knockout analysis.
    selected_taxa (List[str]): A list of taxa to be included in the plot. If None, all taxa are included. Default is None.
    figsize (tuple): A tuple specifying the size of the figure. Default is (10, 10).
    node_size (int): The size of the nodes in the graph. Default is 10000.
    node_color (str): The color of the nodes in the graph. Default is 'grey'.
    edge_width (int): The width of the edges in the graph. Default is 3.
    cooperative_color (str): The color of the cooperative edges in the graph. Default is 'blue'.
    competitive_color (str): The color of the competitive edges in the graph. Default is 'red'.
    arrowsize (int): The size of the arrow heads. Default is 20.

    Returns:
    None
    """
    if selected_taxa is not None:
        ko = ko.loc[selected_taxa, selected_taxa]
    fig, ax = plt.subplots(figsize=figsize)

    G = nx.DiGraph()
    for taxon in ko.columns:
        G.add_node(taxon)
    for taxon1 in ko.columns:
        for taxon2 in ko.columns:
            if taxon1 != taxon2:
                growth_change = ko.loc[taxon1, taxon2]
                if growth_change > 0:
                    G.add_edge(taxon1, taxon2, color=competitive_color)
                elif growth_change < 0:
                    G.add_edge(taxon1, taxon2, color=cooperative_color)

    pos = nx.circular_layout(G)
    node_labels = {node: node for node in G.nodes()}
    nx.draw(
        G,
        pos,
        labels=node_labels,
        font_color="black",
        node_size=node_size,
        node_color=node_color,
        edge_color=nx.get_edge_attributes(G, "color").values(),
        width=edge_width,
        arrowsize=arrowsize,
        connectionstyle="arc3, rad = 0.1",
    )
    ax.set_facecolor("black")
    ax.axis("off")
    fig.set_facecolor("black")
    plt.show()


def plot_trophic_interactions(
    bipartite_graph,
    environmental_carbon_sources,
    ax,
    highlight_compounds=None,
    target_taxon="Acinetobacter",
    target_compound="tol",
    color_oil_nodes="violet",
    color_other_edges="lightgrey",
    color_other_nodes="silver",
    large_node_size=6000,
    small_node_size=1000,
    edge_width_target_taxon=3.5,
    edge_width_other=0.8,
    arrow_size_target_taxon=10,
    arrow_size_other=5,
    font_size=8,
    seed: int = None,
):
    if highlight_compounds is None:
        highlight_compounds = []

    if seed is None:
        seed = random.randint(0, 100)

    taxon_nodes_connected_to_highlight = [
        n
        for compound in highlight_compounds
        for n in bipartite_graph.neighbors(compound)
        if bipartite_graph.nodes[n]["bipartite"] == 0
    ]

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
    target_taxon_compounds = [
        n
        for n in bipartite_graph.to_undirected().neighbors(target_taxon)
        if bipartite_graph.nodes[n]["bipartite"] == 1
    ]
    species_connected_to_target_taxon_compounds = [
        n
        for compound in target_taxon_compounds
        for n in bipartite_graph.neighbors(compound)
        if bipartite_graph.nodes[n]["bipartite"] == 0
    ]
    subgraph_nodes = (
        direct_nodes
        + indirect_nodes
        + indirect_nodes_extended
        + target_taxon_compounds
        + species_connected_to_target_taxon_compounds
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
        rotate=seed,
    )

    non_highlight_edges = [
        (u, v)
        for u, v in extended_subgraph.edges()
        if u not in highlight_compounds and v not in highlight_compounds
    ]
    non_highlight_subgraph = extended_subgraph.edge_subgraph(non_highlight_edges)

    highlight_edges = [
        (u, v)
        for u, v in extended_subgraph.edges()
        if u in highlight_compounds or v in highlight_compounds
    ]
    highlight_subgraph = extended_subgraph.edge_subgraph(highlight_edges)

    nx.draw(
        non_highlight_subgraph,
        shell_layout_extended_subgraph,
        with_labels=False,
        node_size=[
            large_node_size if bipartite == 0 else small_node_size
            for bipartite in nx.get_node_attributes(
                non_highlight_subgraph, "bipartite"
            ).values()
        ],
        node_color=[
            color_oil_nodes
            if node in highlight_compounds
            or node == target_taxon
            or node in taxon_nodes_connected_to_highlight
            else color_other_nodes
            for node, bipartite in nx.get_node_attributes(
                non_highlight_subgraph, "bipartite"
            ).items()
        ],
        edge_color=color_other_edges,
        arrowsize=arrow_size_other,
        width=edge_width_other,
        ax=ax,
    )

    nx.draw(
        highlight_subgraph,
        shell_layout_extended_subgraph,
        with_labels=False,
        node_size=[
            large_node_size if bipartite == 0 else 2 * small_node_size
            for bipartite in nx.get_node_attributes(
                highlight_subgraph, "bipartite"
            ).values()
        ],
        node_color=[
            color_oil_nodes
            if node in highlight_compounds
            or node == target_taxon
            or node in taxon_nodes_connected_to_highlight
            else color_other_nodes
            for node, bipartite in nx.get_node_attributes(
                highlight_subgraph, "bipartite"
            ).items()
        ],
        edge_color=color_oil_nodes,
        arrowsize=arrow_size_target_taxon,
        width=edge_width_target_taxon,
        ax=ax,
    )
    # nx.draw_networkx_edge_labels(
    #     extended_subgraph,
    #     shell_layout_extended_subgraph,
    #     edge_labels=nx.get_edge_attributes(extended_subgraph, "weight"),
    #     ax=ax,
    # )
    highlight_nodes = list(highlight_subgraph.nodes)
    non_highlight_nodes = [
        node for node in extended_subgraph.nodes if node not in highlight_nodes
    ]
    nx.draw_networkx_labels(
        extended_subgraph,
        shell_layout_extended_subgraph,
        labels={n: n for n in non_highlight_nodes},
        font_size=font_size,
        ax=ax,
    )
    nx.draw_networkx_labels(
        extended_subgraph,
        shell_layout_extended_subgraph,
        labels={n: n for n in highlight_nodes},
        font_size=1.5 * font_size,
        ax=ax,
    )

    ax.set_facecolor("black")
    ax.axis("off")
