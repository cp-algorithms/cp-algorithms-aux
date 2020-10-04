import networkx as nx
import matplotlib.pyplot as plt
for seed in range(100):
    G = nx.gnp_random_graph(6, 0.6, seed=seed)

    plt.subplot(111)
    nx.draw(G)   # default spring_layout
    print(seed)
    plt.show()
