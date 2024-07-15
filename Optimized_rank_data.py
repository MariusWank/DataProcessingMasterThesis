import math
import platform
import re
import warnings
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from itertools import combinations
import Welzl
from joblib import Parallel, delayed
import os


# Process a single group of data
def process_group(group_data, number_of_conidia):
    radius = 1.245

    # Find all individual spheres of the branches
    all_hyphae_df = group_data[group_data["cellpart"] == "Hyphae"]
    all_individual_spheres_df = all_hyphae_df.drop_duplicates(subset=["x", "y", "z"])

    # Calculate hyphal length
    hyphal_length = 0
    grouped_hyphae = all_individual_spheres_df.groupby(['agentid', 'cellpart_id'])
    for (agent_id, branch_id), hypha in grouped_hyphae:
        points = hypha[['x', 'y', 'z']].values
        length = np.sum(np.linalg.norm(points[1:] - points[:-1], axis=1))
        hyphal_length += length + 2*radius  # Add diameter to each hypha

    hyphal_length_per_conidium = hyphal_length / number_of_conidia

    # Get number of branches per conidium
    pattern = r"\d+_\d+"
    all_branches_df = group_data[group_data['cellpart_id'].str.match(pattern)]
    all_branches_filtered_by_agentid_df = all_branches_df.drop_duplicates(subset=['agentid', 'cellpart_id'])
    number_of_branches = all_branches_filtered_by_agentid_df.shape[0]
    number_of_branches_per_conidium = number_of_branches / number_of_conidia

    # Get mean highest branch level
    all_branches_df = all_branches_df.copy()
    all_branches_df['branch_level'] = all_branches_df['cellpart_id'].str.extract(r'(\d+)$').astype(int)
    max_branch_levels_per_agentid = all_branches_df.groupby("agentid")['branch_level'].max()
    highest_branch_level_sum = max_branch_levels_per_agentid.sum()
    average_max_branch_level = highest_branch_level_sum / number_of_conidia

    # Get pierce percentage
    all_pierced_branches = all_branches_df[all_branches_df["pierced"] == 1]
    all_pierced_conidia = all_pierced_branches.drop_duplicates(subset=["agentid"])
    pierce_percentage = len(all_pierced_conidia) / number_of_conidia

    # Calculate smallest circle radius average
    circle_radii = []
    all_spheres_by_agentid = group_data.groupby("agentid")
    for agentid, group in all_spheres_by_agentid:

        """#Code to check correctness
        if agentid == 1:
            # Plotting
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111, projection='3d')
        """

        unique_time_array = group["time"].unique()
        min_time = unique_time_array.min()
        starting_points = group[group["time"] == min_time]

        terminal_points = group[group["cellpart_id"] != "basic"].groupby('cellpart_id').tail(1)
        relevant_points = pd.concat([starting_points, terminal_points], ignore_index=True)
        points_list = list(relevant_points[['x', 'y']].drop_duplicates().itertuples(index=False, name=None))


        if len(points_list) > 1:
            circle = Welzl.smallest_enclosing_circle(points_list)
            circle_radii.append(circle.r)

        """#Code to check correctness
        if agentid == 1:
            # Extract x, y, z positions
            x = group['x'].tolist()
            y = group['y'].tolist()
            z = group['z'].tolist()

            marked_x = relevant_points["x"].to_list()
            marked_y = relevant_points["y"].to_list()

            # Scatter plot
            ax.scatter(x, y, z, c='b', marker='o')
            ax.scatter(marked_x, marked_y, c='r', marker='o')

            # Show plot
            plt.show()
        """

    smallest_circle_radius_per_conidium = np.mean(circle_radii)

    # Calculate confinement ratios of every branch
    confinement_ratios = []
    grouped_branches = all_branches_filtered_by_agentid_df.groupby(['agentid', 'cellpart_id'])
    for (agent_id, branch_id), _ in grouped_branches:

        # Filter branch data for the specific agent_id and branch_id
        branch = all_branches_df[(all_branches_df['cellpart_id'] == branch_id) &
                                 (all_branches_df['agentid'] == agent_id)].drop_duplicates(subset=["x", "y", "z"])
        """# Code to check correctness
        if agent_id == 5 and branch_id == "0_1":
            # Plotting
            fig = plt.figure(figsize=(8, 6))
            ax = fig.add_subplot(111, projection='3d')
        """

        branch_points = branch[['x', 'y', 'z']].values
        if len(branch_points) >= 2:
            start_point, end_point = branch_points[0], branch_points[-1]
            euclidean_distance = np.linalg.norm(start_point - end_point)
            branch_length = np.sum(np.linalg.norm(branch_points[1:] - branch_points[:-1], axis=1))
            if branch_length > 0:
                confinement_ratio = euclidean_distance / branch_length
                confinement_ratios.append(confinement_ratio)

                """#Code to check correctness
                if agent_id == 5 and branch_id == "0_1":
                    # Extract x, y, z positions
                    x = branch['x'].tolist()
                    y = branch['y'].tolist()
                    z = branch['z'].tolist()

                    marked_x = [start_point[0], end_point[0]]
                    marked_y = [start_point[1], end_point[1]]
                    marked_z = [start_point[2], end_point[2]]

                    # Scatter plot
                    ax.scatter(x, y, z, c='b', marker='o')
                    ax.scatter(marked_x, marked_y, marked_z, c='r', marker='o')

                    # Show plot
                    plt.show()
                """

    confinement_ratio_per_conidium = np.mean(confinement_ratios)

    return (hyphal_length_per_conidium, number_of_branches_per_conidium, average_max_branch_level,
            smallest_circle_radius_per_conidium, confinement_ratio_per_conidium, pierce_percentage)


# Function to process a folder and generate plots
def process_folder(df):

    # Drop the last column
    df.drop(df.columns[-1], axis=1, inplace=True)

    df.dropna(inplace=True)

    groups = df.groupby(['id'])

    results = Parallel(n_jobs=-1)(delayed(process_group)(group_data, group_data['agentid'].max() + 1)
                                  for group_name, group_data in groups)

    return list(zip(*results))


def main(mean_length, mean_num_branches, mean_branch_level, mean_confinement_ratio, mean_smallest_circle,
         mean_pierce_percentage, root):

    real_mean_values = {
        "hyphal length": mean_length,
        "number of branches": mean_num_branches,
        "branch level": mean_branch_level,
        "confinement ratio": mean_confinement_ratio,
        "smallest circle radius": mean_smallest_circle,
        "pierce percentage": mean_pierce_percentage
    }

    # Read the CSV files
    folders = sorted([file for file in root.iterdir() if file.is_dir() and len(os.listdir(file/"measurements")) != 0])
    faulty_folders = [file for file in root.iterdir() if file.is_dir() and len(os.listdir(file/"measurements")) == 0]

    data_dict = {}

    # Convert DtypeWarning to an exception
    warnings.simplefilter('error', pd.errors.DtypeWarning)
    folders_to_remove = []

    for v, folder in enumerate(folders):
        print("Starting with folder " + str(v + 1))
        try:
            df = pd.read_csv(folder / "measurements" / "agent-statistics.csv", delimiter=';')
        except pd.errors.DtypeWarning as e:
            print(e)
            folders_to_remove.append(folder)
            continue
        data_dict[folder] = process_folder(df)

    # Remove faulty folders
    for folder in folders_to_remove:
        folders.remove(folder)
        faulty_folders.append(folder)

    # Process each folder and generate plots for each category
    categories = ["hyphal length", "number of branches", "branch level", "smallest circle radius", "confinement ratio",
                  "pierce percentage"]
    error_dict = {folder: {category: {} for category in categories} for folder in folders}
    for i, category in enumerate(categories):

        # Check if the category to match exists in the dictionary
        if category in real_mean_values:
            real_mean_value = real_mean_values[category]
        else:
            print("No mean value to match {} was provided!".format(category))
            continue

        # Calculate errors of each folder
        for folder in folders:
            # Get data for category
            category_data = np.array(data_dict[folder][i])

            # Calculate error and add it to list of errors
            mse = np.mean((category_data - real_mean_value) ** 2)
            error_dict[folder][category]["MSE"] = mse

            mae = np.mean(np.abs(category_data - real_mean_value))
            error_dict[folder][category]["MAE"] = mae

            mean_value = np.mean(category_data)
            relative_error = abs(mean_value - real_mean_value) / real_mean_value
            error_dict[folder][category]["Relative Error"] = relative_error

    # Select top lowest error folders and put them in Data_to_Vis folder
    # Step 1: Calculate the sum of "Relative Error" for each folder
    folder_sums = {folder: sum(data.get('Relative Error', 0) for data in categories.values())
                   for folder, categories in error_dict.items()}

    # Step 2: Sort the folders based on the calculated sums
    sorted_folders = sorted(folder_sums.items(), key=lambda item: item[1])

    top = 3

    # Step 3: Get the top lowest folders
    top_lowest = sorted_folders[:top]

    # Step 4: Get the corresponding keys
    top_lowest_keys = [key for key, value in top_lowest]

    # Check the operating system
    system = platform.system()
    destination = Path("Data_to_Vis")

    for key in top_lowest_keys:
        if system == "Windows":
            if os.path.isdir(key):
                # Use xcopy for directories
                os.system(f'xcopy /E /I "{key}" "{destination / key.name}\\"')
            else:
                # Use copy for files
                os.system(f'copy "{key}" "{destination}"')
        else:
            # Use cp for Unix-based systems
            if os.path.isdir(key):
                os.system(f'cp -r "{key}" "{destination}"')
            else:
                os.system(f'cp "{key}" "{destination}"')

    print("Following folders are faulty:")
    for folder in faulty_folders:
        print(folder)
    print("End of faulty folders!")

if __name__ == "__main__":
    MEAN_HYPHAL_LENGTH_DATA = 427.0
    MEAN_NUM_BRANCHES_DATA = 4
    MEAN_BRANCH_LEVEL_DATA = 2.34
    MEAN_CONFINEMENT_RATIO_DATA = 0.9
    MEAN_SMALLEST_CIRCLE_RADIUS_DATA = 119
    MEAN_PIERCE_PERCENTAGE = 0.1238
    #root = Path("/media/mwank/TOSHIBA EXT/Masterarbeit_Data/Test_Folder")
    root = Path("H:\Masterarbeit_Data\FullScreenExcept_cAng0.025_OnlyApi")
    main(MEAN_HYPHAL_LENGTH_DATA, MEAN_NUM_BRANCHES_DATA, MEAN_BRANCH_LEVEL_DATA, MEAN_CONFINEMENT_RATIO_DATA,
         MEAN_SMALLEST_CIRCLE_RADIUS_DATA, MEAN_PIERCE_PERCENTAGE, root)
