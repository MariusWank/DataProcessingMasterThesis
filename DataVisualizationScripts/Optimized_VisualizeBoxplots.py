import math
import re
from pathlib import Path
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
from itertools import combinations
import Welzl
from joblib import Parallel, delayed
from scipy.stats import t


def t_cdf(t_stat, df):
    # Using the scipy.stats t-distribution CDF function
    return t.cdf(t_stat, df)


# Function to decrement the integer part
def decrement_match(match):
    number = int(match.group(1)) - 1
    return f'null_{number}'


# Function to split based on 'null' entries
def split_by_null(df):
    split_dfs = []
    current_split = []
    for idx, row in df.iterrows():
        if 'null' in row['cellpart_id']:
            if current_split:
                split_dfs.append(pd.DataFrame(current_split))
                current_split = []
        else:
            current_split.append(row)
    if current_split:
        split_dfs.append(pd.DataFrame(current_split))
    return split_dfs


def distance(P1, P2):
    p1 = np.array(P1)
    p2 = np.array(P2)
    return np.linalg.norm(p1 - p2)


# Annotate the boxplot with significance brackets
def add_bracket(x1, x2, y, significance_text):
    line_x = [x1, x1, x2, x2]
    line_y = [y, y + 0.01 * y, y + 0.01 * y, y]
    plt.plot(line_x, line_y, lw=1.5, color='black')
    plt.text((x1 + x2) * 0.5, y + 0.02 * y, significance_text, ha='center', va='bottom', color='black')


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

    # Calculate covered area in µm²
    covered_area = hyphal_length * 2 * radius
    hyphal_growth_covered_area_per_conidium = covered_area / number_of_conidia

    # Get number of branches per conidium
    pattern = r"\d+_\d+"
    all_branches_df = group_data[group_data['cellpart_id'].str.match(pattern)]
    all_branches_filtered_by_agentid_df = all_branches_df.drop_duplicates(subset=['agentid', 'cellpart_id'])
    number_of_branches = all_branches_filtered_by_agentid_df.shape[0]
    number_of_branches_per_conidium = number_of_branches / number_of_conidia

    # Get Hyphal Growth Unit
    hyphal_growth_unit = hyphal_length / number_of_branches

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

    # Calculate coverage in %
    welzl_area = math.pi * np.power(smallest_circle_radius_per_conidium, 2)
    mean_coverage = hyphal_growth_covered_area_per_conidium / welzl_area

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

    return (hyphal_length_per_conidium, number_of_branches_per_conidium,
            average_max_branch_level, hyphal_growth_unit,
            hyphal_growth_covered_area_per_conidium, pierce_percentage,
            confinement_ratio_per_conidium, smallest_circle_radius_per_conidium,
            mean_coverage)


# Function to process a folder and generate plots
def process_folder(folder):
    usecols = ['id', "time", "agentid", "x", "y", "z", "cellpart", "cellpart_id", "distance_to_upper_surface",
               "pierced"]
    df = pd.read_csv(folder / "measurements" / "agent-statistics.csv", delimiter=';', usecols=usecols)
    df.dropna(inplace=True)

    groups = df.groupby(['id'])
    results = Parallel(n_jobs=-1)(delayed(process_group)(group_data, group_data['agentid'].max() + 1)
                                  for group_name, group_data in groups)

    return list(zip(*results))


def main(real_mean_values_dict, real_std_values_dict, root, vis_significance, vis_in_vitro_significance):

    # Read the CSV files
    folders = sorted([file for file in root.iterdir() if file.is_dir()])

    # Create dataframe from results of the parameter combination
    df = pd.read_csv(folders[0] / "measurements" / "agent-statistics.csv", delimiter=';')
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
        print("Finished folder " + str(v + 1))

    # Process each folder and generate plots for each category
    categories = ['hyphal length', 'number of branches', 'branch level', 'hyphal growth unit', 'covered area',
                  "pierce percentage", "confinement ratio", "smallest circle radius", "coverage"]

    error_dict = {folder: {category: {} for category in categories} for folder in folders}

    if vis_in_vitro_significance:
        in_vitro_significance_dict = {folder: {category: {} for category in categories} for folder in folders}

    if vis_significance:  #TODO: IMPLEMENT
        rate_of_change_dict = {category: 0 for category in categories}

    for i, category in enumerate(categories):
        # Create Figure
        plt.figure(figsize=(6, 4))
        plt.style.use("seaborn-v0_8-whitegrid")

        # Plot line representing mean value from data
        # Check if the category to match exists in the dictionary
        if category in real_mean_values_dict:
            mean_value = real_mean_values_dict[category]
        else:
            print("No mean value to match {} was provided!".format(category))
            continue
        if not vis_significance:
            data_std_value = real_std_values_dict[category]
            plt.axhline(y=mean_value, color="black", linestyle="--")
            plt.axhline(y=mean_value - data_std_value, color="black", linestyle="dotted")
            plt.axhline(y=mean_value + data_std_value, color="black", linestyle="dotted")

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
                cohens_d = abs((mean1 - mean2) / ((var1 + var2) / 2) ** 0.5)
                t_stat, p_val = stats.ttest_ind(values1, values2, equal_var=False)
                p_values[u, v] = p_val
                p_values[v, u] = p_val
                d_values[u, v] = cohens_d
                d_values[v, u] = cohens_d

        if vis_in_vitro_significance:
            data_size = 60

            if category == "confinement ratio":
                data_size = 250

            for j, folder in enumerate(folders):
                model_values = data_dict[folders[j]][i]
                model_size = len(model_values)
                model_mean = np.mean(model_values)
                model_var = np.var(model_values)
                data_mean = real_mean_values_dict[category]
                data_var = real_std_values_dict[category]**2
                cohens_d = abs((model_mean - data_mean) / ((model_var + data_var) / 2) ** 0.5)

                # Calculate t-statistic
                t_stat = (model_mean - data_mean) / math.sqrt(model_var / model_size + data_var / data_size)

                # Calculate degrees of freedom
                numerator = (model_var / model_size + data_var / data_size) ** 2
                denominator = ((model_var / model_size) ** 2 / (model_size - 1)) + (
                            (data_var / data_size) ** 2 / (data_size - 1))
                df = numerator / denominator

                # Calculate the p-value (two-tailed test)
                p_val = 2 * (1 - t_cdf(abs(t_stat), df))

                # Add values to dictionary
                in_vitro_significance_dict[folder][category]["p_value"] = p_val
                in_vitro_significance_dict[folder][category]["df"] = df
                in_vitro_significance_dict[folder][category]["t_stat"] = t_stat
                in_vitro_significance_dict[folder][category]["Cohens_d"] = cohens_d

        # Plot data in figure
        for j, folder in enumerate(folders):
            category_data = np.array(data_dict[folder][i])
            jittered_x = np.random.normal(j + 1, 0.1, size=len(category_data))
            plt.scatter(jittered_x, category_data, marker='o', alpha=0.5, c="grey", s=5 * 2 ** 0, edgecolors="white")
            plt.boxplot(category_data, positions=[j + 1], widths=0.5, patch_artist=True, showmeans=True,
                        boxprops=dict(color="black", linewidth=2, facecolor="none", alpha=0.5),
                        medianprops=dict(color="black", linewidth=2),
                        meanprops=dict(color="grey", marker="o", markersize=4, markerfacecolor="grey",
                                       markeredgecolor="black", alpha=1))

            if i == 0:
                print("Set {} with parameters ".format(j + 1) + str(folder))

            # Add mean and standard deviation to error dict
            error_dict[folder][category]["Mean"] = np.mean(category_data)
            error_dict[folder][category]["Std"] = np.std(category_data)

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
        plt.xticks(np.arange(1, 3), ["Evasion", "No evasion"])
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
        plt.ylabel(y_label_str)
        plt.ylim(bottom=0, top=mean_value + real_std_values_dict[category] + 0.5 * y_max)

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

    # Select top lowest error folders for the categories to be ranked and put them in Data_to_Vis folder
    rank_categories = ["branch level", "hyphal length", "smallest circle radius", "number of branches",
                       "confinement ratio"]
    # Step 1: Calculate the sum of "Relative Error" for each folder
    folder_sums = {
        folder: sum(data.get('Relative Error', 0)
                    for category, data in categories.items()
                    if category in rank_categories)
        for folder, categories in error_dict.items()
    }

    # Step 2: Sort the folders based on the calculated sums
    sorted_folders = sorted(folder_sums.items(), key=lambda item: item[1])

    top = 3

    # Step 3: Get the top lowest folders
    top_lowest = sorted_folders[:top]

    # Step 4: Get the corresponding keys
    top_lowest_keys = [key for key, value in top_lowest]

    for key in top_lowest_keys:
        print(key)
        print(error_dict[key])
        if vis_in_vitro_significance:
            print(in_vitro_significance_dict[key])


if __name__ == "__main__":
    MEAN_HYPHAL_LENGTH_DATA = 427.0
    MEAN_NUM_BRANCHES_DATA = 4
    MEAN_BRANCH_LEVEL_DATA = 2.34
    MEAN_HGU_DATA = 106.75
    PIERCE_PERCENTAGE_DATA = 0.1238
    MEAN_COVERED_AREA_DATA = 1062
    MEAN_SMALLEST_CIRCLE_RADIUS_DATA = 119
    MEAN_COVERAGE_DATA = 0.0238
    MEAN_CONFINEMENT_RATIO_DATA = 0.9

    real_mean_values = {
        "hyphal length":  427.0,
        "number of branches":  4,
        "branch level": 2.34,
        'hyphal growth unit': 106.75,
        'covered area': 1062,
        "pierce percentage": 0.1238,
        "confinement ratio": 0.9,
        "smallest circle radius": 119,
        "coverage": 0.0238
    }

    real_std_values = {
        "hyphal length": 358,
        "number of branches": 3,
        "branch level": 1.13,
        'hyphal growth unit': 120,
        'covered area': 967,
        "pierce percentage": 0.0515,
        "confinement ratio": 0.1,
        "smallest circle radius": 74,
        "coverage": 0.0368
    }

    visualize_significance = True  # Shows differences between model data
    vis_in_vitro_significance = not visualize_significance  # Shows difference between model data and in_vitro data

    root = Path("Data_to_Vis")
    #root = Path("/media/mwank/TOSHIBA EXT/Masterarbeit_Data/Results/Exporative_NEW/bAng_influence")
    #root = Path("/media/mwank/TOSHIBA EXT/Masterarbeit_Data/Results/Exporative_NEW/bDD_influence")
    #root = Path("/media/mwank/TOSHIBA EXT/Masterarbeit_Data/Results/Exporative_NEW/bLen_influence")
    #root = Path("/media/mwank/TOSHIBA EXT/Masterarbeit_Data/Results/Exporative_NEW/cAng_influence")
    #root = Path(r"H:\Masterarbeit_Data\Results\NewerThanNewest_Results\cAng_Screen_OnlyApi_Winners")
    #root = Path("/media/mwank/TOSHIBA EXT/Masterarbeit_Data/Results/Exporative_NEW/k2_influence")
    #root = Path("/media/mwank/TOSHIBA EXT/Masterarbeit_Data/Results/Exporative_NEW/lSat_influence")
    #root = Path("/media/mwank/TOSHIBA EXT/Masterarbeit_Data/Results/Exporative_NEW/pSym_influence")
    #root = Path("/media/mwank/TOSHIBA EXT/Masterarbeit_Data/Results/Exporative_NEW/confinement_ratio_circle_test/Eva,Surface,cAng")
    #root = Path("/media/mwank/TOSHIBA EXT/Masterarbeit_Data/Results/Exporative_NEW/confinement_ratio_circle_test/Lat_vs_Apical/Apical_Vs_Lateral")
    #root = Path("/media/mwank/TOSHIBA EXT/Masterarbeit_Data/Results/Exporative_NEW/confinement_ratio_circle_test/Lat_vs_Apical/lateral_delay")
    #root = Path("/media/mwank/TOSHIBA EXT/Masterarbeit_Data/Results/Exporative_NEW/confinement_ratio_circle_test/Lat_vs_Apical/Symmetry")
    main(real_mean_values, real_std_values, root, visualize_significance, vis_in_vitro_significance)