# Calculate resting-state functional brain connectivity (fMRI)

1) Load fMRI image
2) Load grey matter (GM) mask
3) Extract time series from all voxels using the GM mask
4) Correlate each pair of time series
5) Threshold correlation (adjacency) matrix
6) Create a graph and obtain graph properties (currently limited to weighted degree centrality)
7) Convert back to NIfTI (.nii) image and save to disk

Jupyter notebook (environment file - Project_3_TIA.yml) and scripts included.