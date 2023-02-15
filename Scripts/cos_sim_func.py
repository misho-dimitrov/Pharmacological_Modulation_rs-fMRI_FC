import sys
from sklearn.metrics.pairwise import cosine_similarity

# Compute cosine similarity (N.B. same as Pearson correlation when data is centered)
def cos_sim_func(time_series):
    """Compute cosine similarity on each pair of voxel timeseries.
    
    The function returns a voxel-wise functional connectivity matrix as calculated by cosine similarity.
    
    Args:
        time_series (ndarray): A brain timeseries array in the format [n_volumes, n_voxels].
        
    Returns:
        cos_sim (ndarray): A connectivity (i.e. adjacency) matrix of the brain voxels in the format [n_voxels, n_voxels].
        
    """
    
    cos_sim = cosine_similarity(time_series.T, time_series.T)
    return cos_sim

if __name__ == "__main__":
    time_series = sys.argv[1]
    cos_sim_func(time_series)


