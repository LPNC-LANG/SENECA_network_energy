import numpy as np

def brain_energy(conn_matrix):
    # Ensure the input is a square matrix
    if not isinstance(conn_matrix, np.ndarray) or conn_matrix.shape[0] != conn_matrix.shape[1]:
        raise ValueError("Input must be a square matrix")
    
    # Initialize variables
    dim_size = conn_matrix.shape[0]
    
    # Precompute the cube root signs and absolute values
    total_energy = 0
    count_triplets = 0
    
    # Iterate over all unique triplets of nodes using vectorization
    for i in range(dim_size - 2):
        submatrix = conn_matrix[i + 1:, i + 1:]
        vec1 = conn_matrix[i, i + 1:]
        
        # Compute energy for all j, k pairs
        triplet_energies = vec1[:, None] * vec1[None, :] * submatrix
        total_energy += np.sum(np.sign(triplet_energies) * np.abs(triplet_energies)**(1/3))
        
        # Count the triplets
        count_triplets += triplet_energies.size
    
    # Return the network energy (average energy over triplets)
    average_energy = -1 * total_energy / count_triplets
    return average_energy


def network_energy_for_networks(conn_matrix, region_df):
    """
    Compute the energy for specific networks based on a region-to-network mapping.

    Parameters:
    - conn_matrix: np.ndarray, square connection matrix
    - region_df: pd.DataFrame, contains columns 'region' and 'network'

    Returns:
    - energy_per_network: dict, mapping of network to its computed energy
    """
    # Ensure conn_matrix is square
    if not isinstance(conn_matrix, np.ndarray) or conn_matrix.shape[0] != conn_matrix.shape[1]:
        raise ValueError("Connection matrix must be a square numpy array.")
    
    # Ensure region_df has the correct structure
    if not all(col in region_df.columns for col in ['index', 'network']):
        raise ValueError("region_df must have 'index' and 'network' columns.")
    
    # Ensure regions in region_df match the size of conn_matrix
    if len(region_df) != conn_matrix.shape[0]:
        raise ValueError("Region DataFrame must have the same number of rows as the size of conn_matrix.")
    
    # Initialize result dictionary
    energy_per_network = {}
    
    # Group regions by network
    for network, group in region_df.groupby('network'):
        # Get the indices for the regions in this network
        indices = group['index'].values
        
        # Extract the submatrix for this network
        submatrix = conn_matrix[np.ix_(indices, indices)]
        
        # Compute the energy for the submatrix
        if len(indices) < 3:
            # Skip networks with fewer than 3 nodes (no triplets possible)
            energy_per_network[network] = None
        else:
            energy_per_network[network] = brain_energy(submatrix)
    
    return energy_per_network
