import argparse
from csv import DictReader
from sys import argv

import networkx as nx
from Bio import Phylo
from networkx.drawing.nx_agraph import write_dot
from networkx.relabel import convert_node_labels_to_integers


def _get_arg_parser():
    parser = argparse.ArgumentParser(
        description='Convert newick to map4sci dot format')
    parser.add_argument(
        '--mapping', help='cell => type label mappings in csv format')
    parser.add_argument('--node-weights', help='How to weight nodes',
                        choices=['pagerank', 'betweenness'], required=True)
    parser.add_argument('input', help='Input newick file.')
    parser.add_argument('output', help='Output dot file')
    return parser


def main(args):
    INPUT = args.input
    OUTPUT = args.output
    MAPPING = args.mapping
    WEIGHTS = args.node_weights
    SCALED_MAX_WEIGHT = 1000

    labels = {}
    if MAPPING:
        for row in DictReader(open(MAPPING)):
            label = row['type']
            if label not in [ 'NA', 'other' ]:
                labels[row['cell']] = label

    # Find the rooted tree
    phylo_tree = None
    for t in Phylo.parse(INPUT, "newick"):
        if t.root:
            phylo_tree = t

    root = phylo_tree.root.name
    tree = Phylo.to_networkx(phylo_tree)

    G = nx.Graph()

    for n in tree.nodes():
        G.add_node(n.name, label=labels.get(str(n.name), n.name), weight=0)

    max_weight = max(weight for (_, _, weight) in tree.edges.data('weight'))
    for (source, target, weight) in tree.edges.data('weight'):
        weight = weight / max_weight * SCALED_MAX_WEIGHT

        # Emphasize edges with a root node or custom-labeled nodes
        if source.name in labels or target.name in labels or source.name == root or target.name == root:
            weight += SCALED_MAX_WEIGHT

        G.add_edge(source.name, target.name, weight=weight)

    if WEIGHTS == 'pagerank':
        pagerank_items = nx.pagerank(G, weight='weight').items()
        max_rank = max(rank for (_n, rank) in pagerank_items)
        for (n, rank) in pagerank_items:
            G.nodes[n]['weight'] = rank / max_rank * SCALED_MAX_WEIGHT
    elif WEIGHTS == 'betweenness':
        beweenness_items = nx.betweenness_centrality(G, weight='weight', normalized=False).items()
        max_rank = max(rank for (_n, rank) in beweenness_items)
        for (n, centrality) in beweenness_items:
            G.nodes[n]['weight'] = centrality / max_rank * SCALED_MAX_WEIGHT

    # Emphasize root node and custom-labeled nodes
    for n in G.nodes():
        if n == root or n in labels:
            G.nodes[n]['weight'] += SCALED_MAX_WEIGHT

    G = G.subgraph(max(nx.connected_components(G), key=len)).copy()
    H = convert_node_labels_to_integers(G, 1, 'decreasing degree', 'cell')

    write_dot(H, OUTPUT)

    print(H.number_of_nodes(), H.number_of_edges())

if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)
