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
    average_energy = -1 * total_energy / count_triplets if count_triplets > 0 else 0
    return average_energy

def brain_energy_within_networks(conn_matrix, region_df):
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
    
    within_network_energy = {}

    # Group regions by network
    region_df = region_df.sort_values(by='index').reset_index(drop=True)
    network_groups = {network: group['index'].values for network, group in region_df.groupby('network')}

    # Calculate within-network energy
    for network, indices in network_groups.items():
        # Extract the submatrix for this network
        submatrix = conn_matrix[np.ix_(indices, indices)]
        
        # Compute the energy for the submatrix
        if len(indices) < 3:
            # Skip networks with fewer than 3 nodes (no triplets possible)
            within_network_energy[network] = None
        else:
            within_network_energy[network] = brain_energy(submatrix)
    
    return within_network_energy

def brain_energy_between_networks(conn_matrix, region_df):
    # Ensure conn_matrix is square
    if not isinstance(conn_matrix, np.ndarray) or conn_matrix.shape[0] != conn_matrix.shape[1]:
        raise ValueError("Connection matrix must be a square numpy array.")
    
    # Ensure region_df has the correct structure
    if not all(col in region_df.columns for col in ['index', 'network']):
        raise ValueError("region_df must have 'index' and 'network' columns.")
    
    # Ensure regions in region_df match the size of conn_matrix
    if len(region_df) != conn_matrix.shape[0]:
        raise ValueError("Region DataFrame must have the same number of rows as the size of conn_matrix.")
    
    between_network_energy = {}

    # Group regions by network
    region_df = region_df.sort_values(by='index').reset_index(drop=True)
    network_groups = {network: group['index'].values for network, group in region_df.groupby('network')}

    # Calculate between-network energy
    for network1, indices1 in network_groups.items():
        for network2, indices2 in network_groups.items():
            if network1 >= network2:
                # Avoid duplicate calculations and self-comparisons
                continue
            
            # Compute the energy for the between-network connections
            # Here, we only consider triplets spanning both networks
            if len(indices1) > 0 and len(indices2) > 0:
                triplet_energy = compute_between_network_energy(conn_matrix, indices1, indices2)
                between_network_energy[(network1, network2)] = triplet_energy

    return between_network_energy

def compute_between_network_energy(conn_matrix, indices1, indices2):
    """
    Compute the energy for triplets spanning two networks.

    Parameters:
    - conn_matrix: np.ndarray, square connection matrix
    - indices1: np.ndarray, indices of nodes in network 1
    - indices2: np.ndarray, indices of nodes in network 2

    Returns:
    - between_network_energy: float, energy of the triplets spanning the two networks
    """
    energy_sum = 0
    count = 0

    # Iterate through triplets with one node in indices1 and two in indices2, or vice versa
    for i in indices1:
        for j, k in combinations(indices2, 2):
            triplet_energy = conn_matrix[i, j] * conn_matrix[i, k] * conn_matrix[j, k]
            energy_sum += np.sign(triplet_energy) * np.abs(triplet_energy)**(1/3)
            count += 1

    for i in indices2:
        for j, k in combinations(indices1, 2):
            triplet_energy = conn_matrix[i, j] * conn_matrix[i, k] * conn_matrix[j, k]
            energy_sum += np.sign(triplet_energy) * np.abs(triplet_energy)**(1/3)
            count += 1

    # Normalize by the number of triplets
    return -1 * energy_sum / count if count > 0 else 0

import pandas as pd
from itertools import combinations

def region_energy(conn_matrix, region_df):
    """
    Compute and normalize the contribution of each region to positive and negative triplet energies.
    (Highly optimized)
    Returns:
    - contribution_df: pd.DataFrame with columns:
        'region', 'positive_contribution', 'negative_contribution',
        'normalized_positive', 'normalized_negative', 'region_energy'
    """
    # Ensure the input is a square matrix
    if not isinstance(conn_matrix, np.ndarray) or conn_matrix.shape[0] != conn_matrix.shape[1]:
        raise ValueError("Input must be a square matrix")
    
    dim_size = conn_matrix.shape[0]
    contributions = np.zeros((dim_size, 2))  # Columns: [positive, negative]
    positive_triplet_counts = np.zeros(dim_size)
    negative_triplet_counts = np.zeros(dim_size)
    
    # Get all combinations of triplets (i, j, k)
    triplet_indices = np.array(list(combinations(range(dim_size), 3)))
    i, j, k = triplet_indices[:, 0], triplet_indices[:, 1], triplet_indices[:, 2]
    
    # Compute triplet energies using vectorized operations
    triplet_energies = conn_matrix[i, j] * conn_matrix[i, k] * conn_matrix[k, j]
    triplet_energies = np.sign(triplet_energies) * np.abs(triplet_energies)**(1/3)
    
    # Separate positive and negative energies
    positive_mask = triplet_energies > 0
    negative_mask = triplet_energies < 0
    
    # Update positive contributions and counts
    np.add.at(contributions[:, 0], i[positive_mask], triplet_energies[positive_mask])
    np.add.at(contributions[:, 0], j[positive_mask], triplet_energies[positive_mask])
    np.add.at(contributions[:, 0], k[positive_mask], triplet_energies[positive_mask])
    np.add.at(positive_triplet_counts, i[positive_mask], 1)
    np.add.at(positive_triplet_counts, j[positive_mask], 1)
    np.add.at(positive_triplet_counts, k[positive_mask], 1)
    
    # Update negative contributions and counts
    np.add.at(contributions[:, 1], i[negative_mask], triplet_energies[negative_mask])
    np.add.at(contributions[:, 1], j[negative_mask], triplet_energies[negative_mask])
    np.add.at(contributions[:, 1], k[negative_mask], triplet_energies[negative_mask])
    np.add.at(negative_triplet_counts, i[negative_mask], 1)
    np.add.at(negative_triplet_counts, j[negative_mask], 1)
    np.add.at(negative_triplet_counts, k[negative_mask], 1)
    
    # Normalize contributions
    normalized_positive = np.divide(
        contributions[:, 0], positive_triplet_counts, out=np.zeros_like(contributions[:, 0]), where=positive_triplet_counts != 0
    )
    normalized_negative = np.divide(
        contributions[:, 1], negative_triplet_counts, out=np.zeros_like(contributions[:, 1]), where=negative_triplet_counts != 0
    )
    
    # Compute instability index
    region_energy = np.divide(
        np.abs(normalized_negative),
        normalized_positive,
        out=np.zeros_like(normalized_negative),
        where=(positive_triplet_counts != 0) & (negative_triplet_counts != 0)
    )
    
    # Convert to DataFrame
    contribution_df = pd.DataFrame({
        'region': region_df["Region"],
        'positive_contribution': contributions[:, 0],
        'negative_contribution': contributions[:, 1],
        'normalized_positive': normalized_positive,
        'normalized_negative': normalized_negative,
        'region_energy': region_energy
    })
    
    return contribution_df

def region_energy_between_networks(conn_matrix, region_df, network_a, network_b):
    """
    Compute the instability index of each region considering only triads spanning two specified networks.

    Parameters:
    - conn_matrix: np.ndarray, square connection matrix.
    - region_df: pd.DataFrame, contains columns 'Region' and 'network'.
    - network_a: str, the name of the first network.
    - network_b: str, the name of the second network.

    Returns:
    - between_network_df: pd.DataFrame with columns:
        'region', 'network', 'positive_contribution', 'negative_contribution',
        'normalized_positive', 'normalized_negative', 'region_energy'
    """
    # Ensure the input is a square matrix
    if not isinstance(conn_matrix, np.ndarray) or conn_matrix.shape[0] != conn_matrix.shape[1]:
        raise ValueError("Input must be a square matrix")
    
    region_df = region_df.sort_values(by='index').reset_index(drop=True)
    
    # Extract regions for the two specified networks
    regions_a = region_df[region_df['network'] == network_a].index.to_numpy()
    regions_b = region_df[region_df['network'] == network_b].index.to_numpy()
    if len(regions_a) == 0 or len(regions_b) == 0:
        raise ValueError(f"Networks {network_a} or {network_b} do not contain any regions.")
    
    # Initialize contributions and counts
    dim_size = conn_matrix.shape[0]
    contributions = np.zeros((dim_size, 2))  # Columns: [positive, negative]
    positive_triplet_counts = np.zeros(dim_size)
    negative_triplet_counts = np.zeros(dim_size)
    
    # Get all combinations of triplets (i, j, k)
    triplet_indices = np.array(list(combinations(range(dim_size), 3)))
    i, j, k = triplet_indices[:, 0], triplet_indices[:, 1], triplet_indices[:, 2]
    
    # Compute triplet energies
    triplet_energies = conn_matrix[i, j] * conn_matrix[i, k] * conn_matrix[k, j]
    triplet_energies = np.sign(triplet_energies) * np.abs(triplet_energies) ** (1/3)
    
    # Filter triads spanning the two networks
    valid_triads = (
        ((np.isin(i, regions_a) & np.isin(j, regions_b)) | # region i in network a and region j in network b or vice versa
         (np.isin(i, regions_b) & np.isin(j, regions_a))) &
        (np.isin(k, np.concatenate([regions_a, regions_b]))) # third region can be whatever
    )
    i, j, k = i[valid_triads], j[valid_triads], k[valid_triads]
    triplet_energies = triplet_energies[valid_triads]
    
    # Separate positive and negative energies
    positive_mask = triplet_energies > 0
    negative_mask = triplet_energies < 0
    
    # Update positive contributions and counts
    np.add.at(contributions[:, 0], i[positive_mask], triplet_energies[positive_mask])
    np.add.at(contributions[:, 0], j[positive_mask], triplet_energies[positive_mask])
    np.add.at(contributions[:, 0], k[positive_mask], triplet_energies[positive_mask])
    np.add.at(positive_triplet_counts, i[positive_mask], 1)
    np.add.at(positive_triplet_counts, j[positive_mask], 1)
    np.add.at(positive_triplet_counts, k[positive_mask], 1)
    
    # Update negative contributions and counts
    np.add.at(contributions[:, 1], i[negative_mask], triplet_energies[negative_mask])
    np.add.at(contributions[:, 1], j[negative_mask], triplet_energies[negative_mask])
    np.add.at(contributions[:, 1], k[negative_mask], triplet_energies[negative_mask])
    np.add.at(negative_triplet_counts, i[negative_mask], 1)
    np.add.at(negative_triplet_counts, j[negative_mask], 1)
    np.add.at(negative_triplet_counts, k[negative_mask], 1)
    
    # Normalize contributions
    normalized_positive = np.divide(
        contributions[:, 0], positive_triplet_counts, out=np.zeros_like(contributions[:, 0]), where=positive_triplet_counts != 0
    )
    normalized_negative = np.divide(
        contributions[:, 1], negative_triplet_counts, out=np.zeros_like(contributions[:, 1]), where=negative_triplet_counts != 0
    )
    
    # Compute instability index
    region_energy = np.divide(
        np.abs(normalized_negative), 
        normalized_positive, 
        out=np.zeros_like(normalized_negative), 
        where=(positive_triplet_counts != 0) & (negative_triplet_counts != 0)
    )
    
    # Convert to DataFrame
    between_network_df = pd.DataFrame({
        'region': region_df['Region'],
        'network': region_df['network'],
        'positive_contribution': contributions[:, 0],
        'negative_contribution': contributions[:, 1],
        'normalized_positive': normalized_positive,
        'normalized_negative': normalized_negative,
        'region_energy': region_energy
    })
    
    # Filter for regions belonging to the specified networks
    return between_network_df[between_network_df['network'].isin([network_a, network_b])].reset_index(drop=True)
