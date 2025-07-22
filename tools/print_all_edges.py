import pickle

# Load and print the pickle file
file_path = 'data/saved_graphs_noncanonical20250627.pickle'
with open(file_path, 'rb') as file:
    graphs = pickle.load(file)

# Find keys and print all graph info and edges
print("Keys in pickle file:", list(graphs.keys()))
print()

for key in graphs.keys():
    graph = graphs[key]
    print(f"Graph: {key}")
    print(f"Nodes: {len(graph.nodes())}")
    print(f"Edges: {len(graph.edges())}")
    print(f"All nodes: {list(graph.nodes())}")
    print("All edges with attributes:")
    
    for u, v, data in graph.edges(data=True):
        print(f"  {u} -- {v}: {data}")
    
    print("-" * 80) 