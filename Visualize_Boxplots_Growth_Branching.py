import math
import re

import imagej
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import os
import Welzl
import scipy.stats as stats
from itertools import combinations


def distance(P1, P2):
    p1 = np.array(P1)
    p2 = np.array(P2)

    return np.linalg.norm(p1 - p2)


# Annotate the boxplot with significance brackets
def add_bracket(x1, x2, y, significance_text):
    line_x = [x1, x1, x2, x2]
    line_y = [y, y+0.01*y, y+0.01*y, y]
    plt.plot(line_x, line_y, lw=1.5, color='black')
    plt.text((x1 + x2) * 0.5, y + 0.02*y, significance_text, ha='center', va='bottom', color='black')


# Function to process a folder and generate plots
def process_folder(folder):

    usecols = ['id', "time", "agentid", "x", "y", "z", "cellpart", "cellpart_id", "distance_to_upper_surface", "pierced"]

    # Create dataframe from results of the parameter combination
    df = pd.read_csv(folder / "measurements" / "agent-statistics.csv", delimiter=';', usecols=usecols)

    # Identify and print rows with NaN values
    nan_rows = df[df.isna().any(axis=1)]
    if not nan_rows.empty:
        print("Rows with NaN values that will be dropped:")
        print(nan_rows)

    # Drop entries with NaN values
    df.dropna(inplace=True)

    df = df.drop_duplicates(subset=["x", "y", "z"])

    # Group by 'id' (= run ID)
    groups = df.groupby(['id'])

    # Collect all sums of relevant validation parameters
    hyphal_length_per_conidium_list = []
    number_of_branches_per_conidium_list = []
    highest_branch_level_per_conidium_list = []
    hyphal_growth_unit_per_conidium_list = []
    hyphal_growth_covered_area_per_conidium_list = []
    pierce_percentage_per_conidium_list = []
    smallest_circle_radius_per_conidium_list = []
    coverage_per_conidium_list = []
    confinement_ratio_per_conidium_list = []

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

        # Calculate covered area in µm²
        covered_area = hyphal_length * 2*radius
        hyphal_growth_covered_area_per_conidium_list.append(covered_area/number_of_conidia)

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

        # calculate smallest circle radius average
        all_spheres_by_agentid = group_data.groupby("agentid")
        circle_radii = []
        # Iterate over each agentid
        for agentid, group in all_spheres_by_agentid:

            """ Code to check correctness
            if group_name == (1,) and agentid == 1:
                # Plotting
                fig = plt.figure(figsize=(8, 6))
                ax = fig.add_subplot(111, projection='3d')
            """

            # get time values
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

            """ Code to check correctness
            if group_name == (1,) and agentid == 1:
                # Extract x, y, z positions
                x = group['x'].tolist()
                y = group['y'].tolist()
                z = group['z'].tolist()

                marked_x = points["x"].to_list()
                marked_y = points["y"].to_list()

                # Scatter plot
                ax.scatter(x, y, z, c='b', marker='o')
                ax.scatter(marked_x, marked_y, c='r', marker='o')

                # Show plot
                plt.show()
            """

        smallest_circle_radius_per_conidium_list.append(np.mean(circle_radii))

        # Calculate coverage in %
        welzl_area = math.pi * np.power(np.mean(circle_radii), 2)
        mean_coverage = (covered_area/number_of_conidia)/welzl_area

        coverage_per_conidium_list.append(mean_coverage)

        # calculate confinement ratios of every branch
        grouped_branches = all_branches_filtered_by_agentid_df.groupby(['agentid', 'cellpart_id'])

        confinement_ratio_per_branch_list = []
        for (agent_id, branch_id), _ in grouped_branches:
            # Filter branch data for the specific agent_id and branch_id
            branch = all_branches_df[(all_branches_df['cellpart_id'] == branch_id) &
                                     (all_branches_df['agentid'] == agent_id)].drop_duplicates(subset=["x", "y", "z"])

            """ Code to check correctness
            if group_name == (1,) and agent_id == 5 and branch_id == "0_1":
                # Plotting
                fig = plt.figure(figsize=(8, 6))
                ax = fig.add_subplot(111, projection='3d')
            """

            if len(branch) >= 2:  #TODO FOR DISCUSSION: Imaris does probably not include branch angles in confinement
                start_point = [branch.iloc[0]["x"], branch.iloc[0]["y"], branch.iloc[0]["z"]]
                end_point = [branch.iloc[-1]["x"], branch.iloc[-1]["y"], branch.iloc[-1]["z"]]
                euclidean_distance = distance(start_point, end_point)

                if euclidean_distance == 0:
                    continue

                # Calculate the actual branch length by summing the distances between consecutive points
                branch_length = 0
                for i in range(1, len(branch)):
                    point1 = [branch.iloc[i - 1]["x"], branch.iloc[i - 1]["y"], branch.iloc[i - 1]["z"]]
                    point2 = [branch.iloc[i]["x"], branch.iloc[i]["y"], branch.iloc[i]["z"]]
                    branch_length += distance(point1, point2)

                confinement_ratio = euclidean_distance / branch_length
                confinement_ratio_per_branch_list.append(confinement_ratio)

                """ Code to check correctness
                if group_name == (1,) and agent_id == 5 and branch_id == "0_1":
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

        confinement_ratio_per_conidium_list.append(np.mean(confinement_ratio_per_branch_list))

    return (hyphal_length_per_conidium_list, number_of_branches_per_conidium_list,
            highest_branch_level_per_conidium_list, hyphal_growth_unit_per_conidium_list,
            hyphal_growth_covered_area_per_conidium_list, pierce_percentage_per_conidium_list,
            confinement_ratio_per_conidium_list, smallest_circle_radius_per_conidium_list,
            coverage_per_conidium_list)


def main(mean_length, mean_num_branches, mean_branch_level, mean_HGU, mean_covered_area, pierce_percentage,
         mean_smallest_circle_radius, mean_coverage, mean_confinement_ratio, root, vis_significance):

    # Read the CSV files
    folders = sorted([file for file in root.iterdir() if file.is_dir()])

    real_mean_values = {
        "hyphal length": mean_length,
        "number of branches": mean_num_branches,
        "branch level": mean_branch_level,
        'hyphal growth unit': mean_HGU,
        'covered area': mean_covered_area,
        "pierce percentage": pierce_percentage,
        "confinement ratio": mean_confinement_ratio,
        "smallest circle radius": mean_smallest_circle_radius,
        "coverage": mean_coverage
    }

    usecols = ["id", "time"]  # Specify only the columns you need

    dtype = {
        "id": "int64",
        'time': 'int64'
    }

    # Create dataframe from results of the parameter combination
    df = pd.read_csv(folders[0] / "measurements" / "agent-statistics.csv", delimiter=';', low_memory=False, dtype=dtype,
                     usecols=usecols)
    number_of_runs = len(df.groupby(['id']))
    unique_time_df = df.drop_duplicates(subset="time")
    unique_time_array = unique_time_df["time"].to_numpy()
    max_time = unique_time_array.max()
    min_time = unique_time_array.min()
    timestep = unique_time_array[1] - unique_time_array[0]

    data_dict = {}
    for v, folder in enumerate(folders):
        print("Starting with folder " + str(v + 1))
        data_dict[folder] = process_folder(folder)

    # Process each folder and generate plots for each category
    categories = ['hyphal length', 'number of branches', 'branch level', 'hyphal growth unit', 'covered area',
                  "pierce percentage", "confinement ratio", "smallest circle radius", "coverage"]

    error_dict = {folder: {category: {} for category in categories} for folder in folders}

    if vis_significance:
        rate_of_change_dict = {category: 0 for category in categories}

    for i, category in enumerate(categories):
        # Create Figure
        plt.figure(figsize=(6, 4))
        plt.style.use("seaborn-v0_8-whitegrid")

        # Plot line representing mean value from data
        # Check if the category to match exists in the dictionary
        if category in real_mean_values:
            mean_value = real_mean_values[category]
        else:
            print("No mean value to match {} was provided!".format(category))
            continue
        if not vis_significance:
            plt.axhline(y=mean_value, color="black", linestyle="--")

        # Create Significance comparison
        if vis_significance:
            p_values = np.ones((len(folders), len(folders)))
            d_values = np.zeros((len(folders), len(folders)))
            for (u, v) in combinations(range(len(folders)), 2):
                values1 = data_dict[folders[u]][i]
                values2 = data_dict[folders[v]][i]
                mean1 = np.mean(values1)
                mean2 = np.mean(values2)
                var1 = np.var(values1)
                var2 = np.var(values2)
                cohens_d = abs((mean1 - mean2) / ((var1 + var2)/2)**0.5)
                t_stat, p_val = stats.ttest_ind(values1, values2, equal_var=False)
                p_values[u, v] = p_val
                p_values[v, u] = p_val
                d_values[u, v] = cohens_d
                d_values[v, u] = cohens_d

        # Plot data in figure
        for j, folder in enumerate(folders):
            category_data = np.array(data_dict[folder][i])
            jittered_x = np.random.normal(j + 1, 0.1, size=len(category_data))
            plt.scatter(jittered_x, category_data, marker='o', alpha=0.5, c="grey", s=5*2**0, edgecolors="white")
            plt.boxplot(category_data, positions=[j + 1], widths=0.5, patch_artist=True, showmeans=True,
                        boxprops=dict(color="black", linewidth=2, facecolor="none", alpha=0.5),
                        medianprops=dict(color="black", linewidth=2),
                        meanprops=dict(color="grey", marker="o", markersize=4, markerfacecolor="grey",
                                       markeredgecolor="black", alpha=1))

            if i == 0:
                print("Set {} with parameters ".format(j+1) + str(folder))

            # Calculate error and add it to list of errors
            mse = np.mean((category_data - mean_value) ** 2)
            error_dict[folder][category]["MSE"] = mse

            mae = np.mean(np.abs(category_data - mean_value))
            error_dict[folder][category]["MAE"] = mae

            in_silico_mean = np.mean(category_data)
            relative_error = abs(in_silico_mean - mean_value) / mean_value
            error_dict[folder][category]["Relative Error"] = relative_error

        # Create formatting for all plots
        y_max = max([max(data_dict[folder][i]) for folder in folders])
        plt.title('Mean {} after {} minutes for {} runs '.format(category, max_time, number_of_runs))
        plt.xticks(range(1, len(folders) + 1), ["Set " + str(k) for k in range(1, len(folders) + 1)])
        plt.xlabel('Parameter Set')
        y_label_str = category.capitalize()
        if category == "hyphal length":
            y_label_str += " in µm"
        elif category == "hyphal growth unit":
            y_label_str += " in µm/branch"
        elif category == 'covered area':
            y_label_str += " in µm²"
        elif category == "smallest circle radius":
            y_label_str += " in µm"
        elif category == "coverage":
            y_label_str += " in %"
        plt.ylabel(y_label_str)
        plt.ylim(bottom=0, top=1.5*y_max)

        if vis_significance:
            # Add significance brackets
            bracket_height = y_max * 0.1
            count = 1
            for (u, v) in combinations(range(len(folders)), 2):
                if p_values[u, v] < 0.05:
                    cohens_d = d_values[u, v]
                    significance = '*'
                    if p_values[u, v] < 0.001:
                        significance = '***'
                    elif p_values[u, v] < 0.01:
                        significance = '**'
                    add_bracket(u + 1, v + 1, y_max + count * bracket_height,
                                significance + " with d=" + str(round(cohens_d, 2)))
                    count += 1

        # Save the plot as an image
        plt.savefig(root / '{}_comparison.png'.format(category.lower().replace(' ', '_')), bbox_inches='tight')
        plt.close()

    if not vis_significance:
        for key in error_dict:
            print(key)
            print(error_dict[key])


if __name__ == "__main__":
    MEAN_HYPHAL_LENGTH_DATA = 427.0
    MEAN_NUM_BRANCHES_DATA = 4
    MEAN_BRANCH_LEVEL_DATA = 2.34
    MEAN_HGU_DATA = 106.0
    PIERCE_PERCENTAGE_DATA = 0.12
    MEAN_COVERED_AREA_DATA = 1062
    MEAN_SMALLEST_CIRCLE_RADIUS_DATA = 119
    MEAN_COVERAGE_DATA = 0.0238
    MEAN_CONFINEMENT_RATIO_DATA = 0.9

    visualize_significance = False

    root = Path("Data_to_Rank")
    #root = Path("/media/mwank/TOSHIBA EXT/Masterarbeit_Data/Exporative_NEW/bDD_k1_k2_interaction/bDD_60_k1")
    main(MEAN_HYPHAL_LENGTH_DATA, MEAN_NUM_BRANCHES_DATA, MEAN_BRANCH_LEVEL_DATA, MEAN_HGU_DATA, MEAN_COVERED_AREA_DATA,
         PIERCE_PERCENTAGE_DATA, MEAN_SMALLEST_CIRCLE_RADIUS_DATA, MEAN_COVERAGE_DATA, MEAN_CONFINEMENT_RATIO_DATA, root, visualize_significance)
