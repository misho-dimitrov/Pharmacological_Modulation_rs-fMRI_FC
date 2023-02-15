import sys
sys.path.append('/usr/lib/python3/dist-packages')
import graph_tool.all as gt
from to_graph_tool import to_graph_tool 

# Get n of degrees for each voxel
def calc_dc(adjacency_matrix):
    """Get n of degrees for each voxel (first option) or get the weighted sum for each voxel (second option).
    
    The function returns a degree centrality array derived from the graph.
    
    Args:
        adjacency_matrix (ndarray): An adjacency matrix in the format [n_voxels, n_voxels].
        
    Returns:
        dc (ndarray): A 1D voxel-wise degree centrality array.
        
    """
    
    g = to_graph_tool(adjacency_matrix)
    
    #Un-weighted
    #dc = g.get_total_degrees([i for i in range(adjacency_matrix.shape[0])])
    
    #Weighted
    weights = g.degree_property_map("total", g.edge_properties['weight'])
    dc_list = [weights[i] for i in range(adjacency_matrix.shape[0])]
    dc = np.array(dc_list)
    
    return dc

if __name__ == "__main__":
    adjacency_matrix = sys.argv[1]
    calc_dc(adjacency_matrix)


