import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path

# Function to process a folder and generate plots
def process_folder(folder):
    # Create dataframe from results of the parameter combination
    df = pd.read_csv(folder / "measurements" / "agent-statistics.csv", delimiter=';')

    # Group by 'id' (= run ID)
    groups = df.groupby(['id'])

    # Collect all sums of relevant validation parameters
    hyphal_length_per_conidium_list = []
    number_of_branches_per_conidium_list = []
    mean_branch_level_per_conidium_list = []
    hyphal_growth_unit_per_conidium_list = []

    # Iterate over individual runs and process them
    for group_name, group_data in groups:
        # Find the maximum value in column 'A'
        number_of_conidia = group_data['agentid'].max() + 1

        # Find all individual spheres of the branches
        all_hyphae_df = group_data[group_data["cellpart"] == "Hyphae"]
        all_individual_spheres_df = all_hyphae_df.drop_duplicates(subset=["x", "y", "z"])

        # Get radius and number of spheres
        radius = float(all_individual_spheres_df["radius"].iloc[0])
        number_of_spheres = all_individual_spheres_df.shape[0]

        # Calculate hyphal_length per condium by considering step size
        hyphal_length = (number_of_spheres * 0.95 * radius)
        hyphal_length_per_conidium_list.append(hyphal_length / number_of_conidia)

        # Get number of branches per conidium by finding unique set of cellpart IDs excluding ID for branch points
        pattern = r"\d+_\d+"
        all_branches_df = group_data[group_data['cellpart_id'].str.match(pattern)].drop_duplicates(subset=['cellpart_id'])
        number_of_branches = all_branches_df.shape[0]
        number_of_branches_per_conidium_list.append(number_of_branches / number_of_conidia)

        # Get Hyphal Growth Unit
        hyphal_growth_unit = hyphal_length / number_of_branches
        hyphal_growth_unit_per_conidium_list.append(hyphal_growth_unit / number_of_conidia)

        # Get mean branch level
        branch_levels = all_branches_df['cellpart_id'].str.extract(r'(\d+)$').astype(int)
        branch_levels_array = branch_levels.values.flatten()
        mean_branch_level = branch_levels_array.sum() / number_of_branches
        mean_branch_level_per_conidium_list.append(mean_branch_level / number_of_conidia)
        print("Finished run " + str(group_name) + "/1000")

    return hyphal_length_per_conidium_list, number_of_branches_per_conidium_list, mean_branch_level_per_conidium_list, hyphal_growth_unit_per_conidium_list

# Read the CSV files
root = Path("DataFiles")
folders = [file for file in root.iterdir() if file.is_dir()]

number_of_runs = 0

# Process each folder and generate plots for each category
categories = ['hyphal length', 'number of branches', 'branch level', 'hyphal growth unit']
for i, category in enumerate(categories):
    plt.figure(figsize=(6, 4))
    for j, folder in enumerate(folders):
        print("Starting with folder " + str(j+1) + " for " + category)
        folder_name = "parameter_set_" + str(j + 1)
        data = process_folder(folder)
        number_of_runs = len(data[0])
        jittered_x = np.random.normal(j + 1, 0.1, size=len(data[i]))
        plt.scatter(jittered_x, data[i], marker='o', alpha=0.5, c="grey", s=5*2**0,edgecolors="white")
        plt.boxplot(data[i], positions=[j + 1], widths=0.5, patch_artist=True, showmeans=True,
                    boxprops=dict(color="black", linewidth=2, facecolor="none", alpha=0.5),
                    medianprops=dict(color="black", linewidth=2),
                    meanprops=dict(color="grey", marker="o", markersize=4, markerfacecolor="grey",
                                   markeredgecolor="black", alpha=1))


    plt.title('Mean {} over {} runs '.format(category, number_of_runs))
    plt.xticks(range(1, len(folders) + 1), ["Set " + str(k) for k in range(1, len(folders) + 1)])
    plt.xlabel('Parameter Set')
    plt.ylabel(category.capitalize())
    if category != "branch level":
        plt.ylim(bottom=0)
    else:
        plt.ylim(bottom=1.5, top=4)

    # Save the plot as an image
    plt.savefig(root / '{}_comparison.png'.format(category.lower().replace(' ', '_')), bbox_inches='tight')

    plt.close()
