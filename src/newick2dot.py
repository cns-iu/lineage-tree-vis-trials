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

    labels = {}
    if MAPPING:
        for row in DictReader(open(MAPPING)):
            labels[row['cell']] = row['type']

    phylo_tree = Phylo.read(INPUT, "newick")
    tree = Phylo.to_networkx(phylo_tree)

    G = nx.Graph()

    for n in tree.nodes():
        G.add_node(n.name, label=labels.get(str(n.name), n.name), weight=0)

    max_weight = max(weight for (_, _, weight) in tree.edges.data('weight'))
    for (source, target, weight) in tree.edges.data('weight'):
        G.add_edge(source.name, target.name, weight=weight / max_weight * 1000)

    if WEIGHTS == 'pagerank':
        pagerank_items = nx.pagerank(G, weight='weight').items()
        max_rank = max(rank for (_n, rank) in pagerank_items)
        for (n, rank) in pagerank_items:
            G.nodes[n]['weight'] = rank / max_rank * 1000
    elif WEIGHTS == 'betweenness':
        for (n, centrality) in nx.betweenness_centrality(G, weight='weight', normalized=False).items():
            G.nodes[n]['weight'] = centrality

    G = G.subgraph(max(nx.connected_components(G), key=len)).copy()
    H = convert_node_labels_to_integers(G, 1, 'decreasing degree', 'cell')

    write_dot(H, OUTPUT)


if __name__ == "__main__":
    parser = _get_arg_parser()
    args = parser.parse_args()
    main(args)
