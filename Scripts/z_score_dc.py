import sys
import numpy
from sklearn.preprocessing import StandardScaler

# Z-score
def z_score_dc(dc_array):
    dc = dc_array.reshape(-1, 1)
    scaler = StandardScaler()
    dc_scaled = scaler.fit_transform(dc)
    return dc_scaled

if __name__ == "__main__":
    dc_array = sys.argv[1]
    z_score_dc(dc_array)


