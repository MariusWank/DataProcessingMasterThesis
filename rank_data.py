import os

import numpy as np
import pandas as pd
from pathlib import Path
from Visualize_Boxplots_Growth_Branching import main as boxplot_vis

# Function to process a folder and generate plots
def process_folder(folder):
    # Create dataframe from results of the parameter combination
    df = pd.read_csv(folder / "measurements" / "agent-statistics.csv", delimiter=';')

    # Group by 'id' (= run ID)
    groups = df.groupby(['id'])

    # Collect all sums of relevant validation parameters
    hyphal_length_per_conidium_list = []
    number_of_branches_per_conidium_list = []
    highest_branch_level_per_conidium_list = []
    hyphal_growth_unit_per_conidium_list = []
    hyphal_growth_coverage_per_conidium_list = []

    # Iterate over individual runs and process them
    for group_name, group_data in groups:
        # Find the maximum value in column 'agentid'
        number_of_conidia = group_data['agentid'].max() + 1

        # Find all individual spheres of the branches
        all_hyphae_df = group_data[group_data["cellpart"] == "Hyphae"]
        all_individual_spheres_df = all_hyphae_df.drop_duplicates(subset=["x", "y", "z"])

        # Get radius and number of spheres
        radius = float(all_individual_spheres_df["radius"].iloc[0])

        number_of_spheres = all_individual_spheres_df.shape[0]

        # Calculate hyphal_length per conidium by considering step size
        hyphal_length = (2 * radius + (number_of_spheres-1) * 0.95 * radius)
        hyphal_length_per_conidium_list.append(hyphal_length / number_of_conidia)

        # Get number of branches per conidium by finding unique set of cellpart IDs excluding ID for branch points
        pattern = r"\d+_\d+"
        all_branches_df = group_data[group_data['cellpart_id'].str.match(pattern)]
        all_branches_filtered_by_agentid_df = (all_branches_df.drop_duplicates(subset=['agentid', 'cellpart_id']).
                                               reset_index(drop=True))
        number_of_branches = all_branches_filtered_by_agentid_df.shape[0]
        number_of_branches_per_conidium_list.append(number_of_branches / number_of_conidia)

        # Get Hyphal Growth Unit
        hyphal_growth_unit = hyphal_length / number_of_branches
        hyphal_growth_unit_per_conidium_list.append(hyphal_growth_unit)

        # Calculate coverage in µm²
        coverage = hyphal_length * 2*radius
        hyphal_growth_coverage_per_conidium_list.append(coverage/number_of_conidia)

        # Get mean highest branch level
        all_branches_df = all_branches_df.copy()
        all_branches_df['branch_level'] = all_branches_df['cellpart_id'].str.extract(r'(\d+)$').astype(int)
        grouped_by_conidia = all_branches_df.groupby("agentid")
        max_branch_levels_per_agentid = grouped_by_conidia['branch_level'].max()
        highest_branch_level_sum = max_branch_levels_per_agentid.sum()
        average_max_branch_level = highest_branch_level_sum / number_of_conidia
        highest_branch_level_per_conidium_list.append(average_max_branch_level)

    return (hyphal_length_per_conidium_list, number_of_branches_per_conidium_list,
            highest_branch_level_per_conidium_list, hyphal_growth_unit_per_conidium_list,
            hyphal_growth_coverage_per_conidium_list)


def main(mean_length, mean_num_branches, mean_branch_level, mean_HGU):
    # Read the CSV files
    root = Path("Data_to_Rank")
    folders = [file for file in root.iterdir() if file.is_dir() and len(os.listdir(file/"measurements")) != 0]
    faulty_folders = [file for file in root.iterdir() if file.is_dir() and len(os.listdir(file/"measurements")) == 0]

    summed_errors = np.zeros(len(folders))
    data_list = []

    for v, folder in enumerate(folders):
        print("Starting with folder " + str(v + 1))
        data_list.append(process_folder(folder))

    # Process each folder and generate plots for each category
    categories = ['hyphal length', 'number of branches', 'branch level', 'hyphal growth unit']
    for i, category in enumerate(categories):

        # Get mean value of category
        real_mean_value = 0.0
        match category:
            case "hyphal length":
                real_mean_value = mean_length
            case "number of branches":
                real_mean_value = mean_num_branches
            case "branch level":
                real_mean_value = mean_branch_level
            case "hyphal growth unit":
                real_mean_value = mean_HGU

        # Calculate errors of each folder
        for j, folder in enumerate(folders):
            # Get data for category
            category_data = np.array(data_list[j][i])

            # Calculate error and add it to list of errors
            mean_value = np.mean(category_data)
            error = abs(mean_value - real_mean_value)/real_mean_value
            summed_errors[j] += error

    # Select top lowest error folders and put them in Data_to_Vis folder
    top = 3
    indices = np.argpartition(summed_errors, top)
    data_list_transfer = [None] * top
    for i, index in enumerate(indices[:top]):
        src = folders[index]
        os.system('cp -r {} Data_to_Vis'.format(src))
        data_list_transfer[i] = data_list[index]

    print("Following folders are faulty:")
    for folder in faulty_folders:
        print(folder)
    print("End of faulty folders!")

    # boxplot_vis(mean_length, mean_num_branches, mean_branch_level, mean_HGU, data_list_transfer)

if __name__ == "__main__":
    MEAN_HYPHAL_LENGTH_DATA = 427.0
    MEAN_NUM_BRANCHES_DATA = 4
    MEAN_BRANCH_LEVEL_DATA = 2.34
    MEAN_HGU_DATA = 106.0
    main(MEAN_HYPHAL_LENGTH_DATA, MEAN_NUM_BRANCHES_DATA, MEAN_BRANCH_LEVEL_DATA, MEAN_HGU_DATA)
