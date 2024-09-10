import os
import warnings

import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as plt
import Welzl
import plotly.graph_objs as go
import plotly.express as px


def distance(P1, P2):
    p1 = np.array(P1)
    p2 = np.array(P2)

    return np.linalg.norm(p1 - p2)


# Function to process a folder and generate plots
def process_folder(df):

    # Drop the last column
    df.drop(df.columns[-1], axis=1, inplace=True)

    # Drop entries with NaN values
    df.dropna(inplace=True)

    # Group by 'id' (= run ID)
    groups = df.groupby(['id'])

    # Collect all sums of relevant validation parameters
    hyphal_length_per_conidium_list = []
    number_of_branches_per_conidium_list = []
    highest_branch_level_per_conidium_list = []
    confinement_ratio_per_conidium_list = []
    smallest_circle_radius_per_conidium_list = []

    # Iterate over individual runs and process them
    for group_name, group_data in groups:
        # Find the maximum value in column 'agentid'
        number_of_conidia = group_data['agentid'].max() + 1

        # Find all individual spheres of the branches
        all_hyphae_df = group_data[group_data["cellpart"] == "Hyphae"]
        all_individual_spheres_df = all_hyphae_df.drop_duplicates(subset=["x", "y", "z"])

        # Get radius and number of spheres
        radius = 1.245

        number_of_spheres = all_individual_spheres_df.shape[0]

        # Calculate hyphal_length per conidium by considering step size
        step_length = 0.95 * radius
        hyphal_length = (2 * radius + (number_of_spheres-1) * step_length)
        hyphal_length_per_conidium_list.append(hyphal_length / number_of_conidia)

        # Get number of branches per conidium by finding unique set of cellpart IDs excluding ID for branch points
        pattern = r"\d+_\d+"
        group_data = group_data.dropna(subset=['cellpart_id'])
        group_data['cellpart_id'] = group_data['cellpart_id'].astype(str)
        all_branches_df = group_data[group_data['cellpart_id'].str.match(pattern)]

        all_branches_filtered_by_agentid_df = (all_branches_df.drop_duplicates(subset=['agentid', 'cellpart_id']).
                                               reset_index(drop=True))
        number_of_branches = all_branches_filtered_by_agentid_df.shape[0]
        number_of_branches_per_conidium_list.append(number_of_branches / number_of_conidia)

        # Get mean highest branch level
        all_branches_df = all_branches_df.copy()
        all_branches_df.loc[:, 'branch_level'] = all_branches_df['cellpart_id'].str.extract(r'(\d+)$').astype(int)
        grouped_by_conidia = all_branches_df.groupby("agentid")
        max_branch_levels_per_agentid = grouped_by_conidia['branch_level'].max()
        highest_branch_level_sum = max_branch_levels_per_agentid.sum()
        average_max_branch_level = highest_branch_level_sum / number_of_conidia
        highest_branch_level_per_conidium_list.append(average_max_branch_level)

        # calculate smallest circle radius average
        all_spheres_by_agentid = group_data.groupby("agentid")
        circle_radii = []
        # Iterate over each agentid
        for agentid, group in all_spheres_by_agentid:

            # get min and max time values
            unique_time_df = group.drop_duplicates(subset="time")
            unique_time_array = unique_time_df["time"].to_numpy()
            min_time = unique_time_array.min()
            max_time = all_hyphae_df[all_hyphae_df["agentid"] == agentid]["time"].max()  # must not include mothercells

            starting_points = group[group["time"] == min_time]
            terminal_points = group[group["time"] == max_time]
            relevant_points = pd.concat([starting_points, terminal_points], ignore_index=True)

            # Extract the points (Pt Position X, Pt Position Y)
            points = relevant_points[['x', 'y']].drop_duplicates()

            # Convert to a list of tuples
            points_list = list(points.itertuples(index=False, name=None))

            # Skip trivial radius 0 circles
            if len(points_list) <= 1:
                continue

            # Calculate the smallest circle for hypha
            circle = Welzl.smallest_enclosing_circle(points_list)

            circle_radii.append(circle.r)

        smallest_circle_radius_per_conidium_list.append(np.mean(circle_radii))

        # calculate confinement ratios of every branch
        grouped_branches = all_branches_filtered_by_agentid_df.groupby(['agentid', 'cellpart_id'])

        confinement_ratio_per_branch_list = []
        for (agent_id, branch_id), _ in grouped_branches:
            # Filter branch data for the specific agent_id and branch_id
            branch = all_branches_df[(all_branches_df['cellpart_id'] == branch_id) &
                                     (all_branches_df['agentid'] == agent_id)]

            if len(branch) >= 2:
                start_point = [branch.iloc[0]["x"], branch.iloc[0]["y"], branch.iloc[0]["z"]]
                end_point = [branch.iloc[-1]["x"], branch.iloc[-1]["y"], branch.iloc[-1]["z"]]
                euclidean_distance = distance(start_point, end_point)

                if euclidean_distance == 0:
                    continue

                number_of_spheres_in_branch = branch.shape[0]
                branch_length = (2 * radius + (number_of_spheres_in_branch - 1) * step_length)

                confinement_ratio = euclidean_distance / branch_length
                confinement_ratio_per_branch_list.append(confinement_ratio)

        confinement_ratio_per_conidium_list.append(np.mean(confinement_ratio_per_branch_list))

    return (hyphal_length_per_conidium_list, number_of_branches_per_conidium_list,
            highest_branch_level_per_conidium_list, confinement_ratio_per_conidium_list,
            smallest_circle_radius_per_conidium_list)


def main(mean_length, mean_num_branches, mean_branch_level, mean_confinement_ratio, mean_smallest_circle, root):

    real_mean_values = {
        "hyphal length": mean_length,
        "number of branches": mean_num_branches,
        "branch level": mean_branch_level,
        "confinement ratio": mean_confinement_ratio,
        "smallest circle radius": mean_smallest_circle
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
        except pd.errors.DtypeWarning:
            folders_to_remove.append(folder)
            continue
        data_dict[folder] = process_folder(df)

    # Remove faulty folders
    for folder in folders_to_remove:
        folders.remove(folder)
        faulty_folders.append(folder)

    # Process each folder and generate plots for each category
    categories = ["hyphal length", "number of branches", "branch level", "confinement ratio", "smallest circle radius"]
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

    for key in top_lowest_keys:
        os.system('cp -r "{}" Data_to_Vis'.format(key))
        print(key)
        print(error_dict[key])

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
    root = Path("/media/mwank/TOSHIBA EXT/Masterarbeit_Data/Test_Folder")
    #root = Path("/media/mwank/TOSHIBA EXT/Masterarbeit_Data/FullScreenExcept_cAng0.025_OnlyApi")
    main(MEAN_HYPHAL_LENGTH_DATA, MEAN_NUM_BRANCHES_DATA, MEAN_BRANCH_LEVEL_DATA, MEAN_CONFINEMENT_RATIO_DATA,
         MEAN_SMALLEST_CIRCLE_RADIUS_DATA, root)
