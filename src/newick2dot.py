from sys import argv

import networkx as nx
from Bio import Phylo
from networkx.drawing.nx_agraph import write_dot
from networkx.relabel import convert_node_labels_to_integers

INPUT=argv[1]
OUTPUT=argv[2]

phylo_tree = Phylo.read(INPUT, "newick")
tree = Phylo.to_networkx(phylo_tree)

G = nx.Graph()

for (n, data) in tree.nodes.data():
  G.add_node(n.name, label=n.name, weight=0)

max_weight = max( weight for (_, _, weight) in tree.edges.data('weight') )
for (source, target, weight) in tree.edges.data('weight'):
  G.add_edge(source.name, target.name, weight=weight / max_weight * 1000)

print(G.number_of_nodes(), G.number_of_edges())

pagerank_items = nx.pagerank(G, weight='weight').items();
max_rank = max( rank for (_n, rank) in pagerank_items)
for (n, rank) in pagerank_items:
  G.nodes[n]['weight'] = rank / max_rank * 1000

# for (n, centrality) in nx.betweenness_centrality(G, weight='weight', normalized=False).items():
#   G.nodes[n]['weight'] = centrality

G = G.subgraph(max(nx.connected_components(G), key=len)).copy()
H = convert_node_labels_to_integers(G, 1, 'decreasing degree', 'label')

write_dot(H, OUTPUT)
