import sys
from sklearn.metrics.pairwise import cosine_similarity

# Compute cosine similarity (N.B. same as Pearson correlation when data is centered)
def cos_sim_func(time_series):
    cos_sim = cosine_similarity(time_series.T, time_series.T)
    return cos_sim

if __name__ == "__main__":
    time_series = sys.argv[1]
    cos_sim_func(time_series)


