import sys
import graph_tool.all as gt
import numpy as np

# https://carlonicolini.github.io/sections/science/2018/09/12/weighted-graph-from-adjacency-matrix-in-graph-tool.html
# To create an undirected, weighted graph with graph-tool
def to_graph_tool(adj):
    """Create an undirected, weighted graph with graph-tool from the adjacency matrix.
    
    The function returns an undirected, weighted graph derived from the adjacency matrix.
    
    Args:
        adj (ndarray): An adjacency matrix in the format [n_voxels, n_voxels].
        
    Returns:
        g (graph_tool.Graph): A graph object.
        
    """
    
    g = gt.Graph(directed=False)
    edge_weights = g.new_edge_property('double')
    g.edge_properties['weight'] = edge_weights
    # Set the lower triangle of the adjacency matrix and the diagonal to 0
    nnz = np.nonzero(np.triu(adj,1))
    # Get the number of edges (i.e. non-zero values)
    nedges = len(nnz[0])
    # Create the edge value list
    g.add_edge_list(
        # Create rows of THREE values
        np.hstack(
        [
        # Transpose nnz so that you have TWO values in each row, where
        # the first is the row index and the second is the column index
        # of this particular edge
        np.transpose(nnz),
        # Get the 3RD values for each row, i.e. the edge weight
        np.reshape(adj[nnz],(nedges,1))
        ]
        ),
    # If given, eprops should specify an iterable containing edge property maps that will be filled with the remaining values at each row, if there are more than two.
    eprops=[edge_weights])
    return g

if __name__ == "__main__":
    adj = sys.argv[1]
    to_graph_tool(adj)


