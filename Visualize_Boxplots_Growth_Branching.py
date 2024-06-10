import math
import re

import imagej
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import os


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
    pierce_percentage_per_conidium_list = []

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
        step_length = 0.95 * radius
        hyphal_length = (2 * radius + (number_of_spheres-1) * step_length)
        hyphal_length_per_conidium_list.append(hyphal_length / number_of_conidia)

        # Calculate coverage in µm² on upper surface
        all_sphere_on_upper_surface_df = all_individual_spheres_df[all_individual_spheres_df["distance_to_upper_surface"] != 0]
        number_of_spheres_on_upper_surface = all_sphere_on_upper_surface_df.shape[0]
        upper_hyphal_length = (2 * radius + (number_of_spheres_on_upper_surface-1) * step_length)
        coverage = upper_hyphal_length * 2*radius
        hyphal_growth_coverage_per_conidium_list.append(coverage/number_of_conidia)

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

        # Get mean highest branch level
        all_branches_df = all_branches_df.copy()
        all_branches_df.loc[:, 'branch_level'] = all_branches_df['cellpart_id'].str.extract(r'(\d+)$').astype(int)
        grouped_by_conidia = all_branches_df.groupby("agentid")
        max_branch_levels_per_agentid = grouped_by_conidia['branch_level'].max()
        highest_branch_level_sum = max_branch_levels_per_agentid.sum()
        average_max_branch_level = highest_branch_level_sum / number_of_conidia
        highest_branch_level_per_conidium_list.append(average_max_branch_level)

        # Get pierce percentage
        all_pierced_branches = all_branches_df[all_branches_df["pierced"] == 1]
        all_pierced_conidia = all_pierced_branches.drop_duplicates(
            subset=["agentid"])
        pierce_percentage = len(all_pierced_conidia) / number_of_conidia
        pierce_percentage_per_conidium_list.append(pierce_percentage)


    return (hyphal_length_per_conidium_list, number_of_branches_per_conidium_list,
            highest_branch_level_per_conidium_list, hyphal_growth_unit_per_conidium_list,
            hyphal_growth_coverage_per_conidium_list, pierce_percentage_per_conidium_list)


def main(mean_length, mean_num_branches, mean_branch_level, mean_HGU, pierce_percentage, data_list_transfer=None):
    # Read the CSV files
    root = Path("Data_to_Vis")
    folders = [file for file in root.iterdir() if file.is_dir()]

    # Create dataframe from results of the parameter combination
    df = pd.read_csv(folders[0] / "measurements" / "agent-statistics.csv", delimiter=';')
    number_of_runs = len(df.groupby(['id']))
    unique_time_df = df.drop_duplicates(subset="time")
    unique_time_array = unique_time_df["time"].to_numpy()
    max_time = unique_time_array.max()
    min_time = unique_time_array.min()
    timestep = unique_time_array[1] - unique_time_array[0]

    data_list = []
    if data_list_transfer is None:
        for v, folder in enumerate(folders):
            print("Starting with folder " + str(v + 1))
            data_list.append(process_folder(folder))
    else:
        data_list = data_list_transfer

    # Process each folder and generate plots for each category
    categories = ['hyphal length', 'number of branches', 'branch level', 'hyphal growth unit', 'coverage', "pierce percentage"]
    for i, category in enumerate(categories):
        # Create Figure
        plt.figure(figsize=(6, 4))
        plt.style.use("seaborn-v0_8-whitegrid")

        # Plot line representing mean value from data
        mean_value = 0.0
        match category:
            case "hyphal length":
                mean_value = mean_length
            case "number of branches":
                mean_value = mean_num_branches
            case "branch level":
                mean_value = mean_branch_level
            case "hyphal growth unit":
                mean_value = mean_HGU
            case "pierce percentage":
                mean_value = pierce_percentage

        if category != "coverage":
            plt.axhline(y=mean_value, color="black", linestyle="--")

        # Plot data in figure
        for j, folder in enumerate(folders):
            data = data_list[j]
            number_of_runs = len(data[0])
            jittered_x = np.random.normal(j + 1, 0.1, size=len(data[i]))
            plt.scatter(jittered_x, data[i], marker='o', alpha=0.5, c="grey", s=5*2**0,edgecolors="white")
            plt.boxplot(data[i], positions=[j + 1], widths=0.5, patch_artist=True, showmeans=True,
                        boxprops=dict(color="black", linewidth=2, facecolor="none", alpha=0.5),
                        medianprops=dict(color="black", linewidth=2),
                        meanprops=dict(color="grey", marker="o", markersize=4, markerfacecolor="grey",
                                       markeredgecolor="black", alpha=1))

            if i == 0:
                print("Set {} with parameters ".format(j+1) + str(folder))

        # Create formatting for all plots
        plt.title('Mean {} after {} minutes for {} runs '.format(category, max_time,number_of_runs))
        plt.xticks(range(1, len(folders) + 1), ["Set " + str(k) for k in range(1, len(folders) + 1)])
        plt.xlabel('Parameter Set')
        y_label_str = category.capitalize()
        if category == "hyphal length":
            y_label_str += " in µm"
        elif category == "hyphal growth unit":
            y_label_str += " in µm/branch"
        elif category == 'coverage':
            y_label_str += " in µm²"
        plt.ylabel(y_label_str)
        if category != "branch level":
            plt.ylim(bottom=0)
        else:
            plt.ylim(bottom=0, top=4)

        # Save the plot as an image
        plt.savefig(root / '{}_comparison.png'.format(category.lower().replace(' ', '_')), bbox_inches='tight')
        plt.close()


if __name__ == "__main__":
    MEAN_HYPHAL_LENGTH_DATA = 427.0
    MEAN_NUM_BRANCHES_DATA = 4
    MEAN_BRANCH_LEVEL_DATA = 2.34
    MEAN_HGU_DATA = 106.0
    PIERCE_PERCENTAGE = 0.12
    main(MEAN_HYPHAL_LENGTH_DATA, MEAN_NUM_BRANCHES_DATA, MEAN_BRANCH_LEVEL_DATA, MEAN_HGU_DATA, PIERCE_PERCENTAGE)
