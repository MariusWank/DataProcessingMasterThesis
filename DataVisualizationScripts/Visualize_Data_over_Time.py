import os.path

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path


# Function to process a folder and generate plots
def process_folder(folder, num_of_timesteps):
    # Create dataframe from results of the parameter combination
    df = pd.read_csv(folder / "measurements" / "agent-statistics.csv", delimiter=';')

    if df["time"].max() % 30 != 0:
        num_of_timesteps -= 1

    # Group by 'id' (= run ID)
    run_groups = df.groupby(['id'])

    # Collect all relevant validation parameters
    summed_mean_length_over_time = np.zeros(num_of_timesteps)
    summed_mean_branches_over_time = np.zeros(num_of_timesteps)
    summed_pierce_percent_over_time = np.zeros(num_of_timesteps)

    mean_branches_per_timestep = [np.zeros(len(run_groups)) for _ in range(num_of_timesteps)]
    mean_length_per_timestep = [np.zeros(len(run_groups)) for _ in range(num_of_timesteps)]
    mean_pierce_percent_per_timestep = [np.zeros(len(run_groups)) for _ in range(num_of_timesteps)]

    # Iterate over individual runs and process them
    for i, (run_group_name, run_group_data) in enumerate(run_groups):
        # Find the maximum value in column 'agentid'
        number_of_conidia = run_group_data['agentid'].max() + 1

        # Collect all relevant validation parameters
        mean_length_over_time = np.zeros(num_of_timesteps)
        mean_branches_over_time = np.zeros(num_of_timesteps)
        pierce_percent_over_time = np.zeros(num_of_timesteps)

        # Group by timesteps
        time_groups = run_group_data.groupby(["time"])
        for j, (time_group_name, time_group_data) in enumerate(time_groups):
            # skip last timestep if not divisible by timestep
            time = time_group_data["time"].max()
            if time != 0 and time % 30 != 0:
                continue

            # Find all individual spheres of the branches
            all_hyphae_df = time_group_data[time_group_data["cellpart"] == "Hyphae"]
            all_individual_spheres_df = all_hyphae_df.drop_duplicates(subset=["x", "y", "z"])

            if all_individual_spheres_df.empty:
                continue

            # Get radius and number of spheres
            radius = float(all_individual_spheres_df["radius"].iloc[0])
            number_of_spheres = all_individual_spheres_df.shape[0]

            # Calculate hyphal_length per condium by considering step size
            hyphal_length = (2 * radius + (number_of_spheres - 1) * 0.95 * radius)
            mean_length_over_time[j:] += (hyphal_length / number_of_conidia)
            for k in range(j, num_of_timesteps):
                mean_length_per_timestep[k][i] += hyphal_length / number_of_conidia

            # Get number of branches per conidium by finding unique set of cellpart IDs excluding ID for branch points
            pattern = r"\d+_\d+"
            all_branches_df = run_group_data[run_group_data["time"] <= time]  # Get all branches up to this point
            all_branches_df = all_branches_df[all_branches_df['cellpart_id'].str.match(pattern)]
            all_branches_filtered_by_agentid_df = (all_branches_df.drop_duplicates(subset=['agentid', 'cellpart_id']).
                                                   reset_index(drop=True))
            number_of_branches = all_branches_filtered_by_agentid_df.shape[0]
            mean_branches_over_time[j] += (number_of_branches / number_of_conidia)

            # Get all pierced branches
            all_pierced_branches = all_branches_df[all_branches_df["pierced"] == 1]
            all_pierced_conidia = all_pierced_branches.drop_duplicates(
                subset=["agentid"])
            pierce_in_percent = len(all_pierced_conidia) / number_of_conidia
            pierce_percent_over_time[j] += pierce_in_percent

            mean_branches_per_timestep[j][i] += (number_of_branches / number_of_conidia)
            mean_pierce_percent_per_timestep[j][i] += pierce_in_percent

        summed_mean_branches_over_time += mean_branches_over_time
        summed_mean_length_over_time += mean_length_over_time
        summed_pierce_percent_over_time += pierce_percent_over_time

    summed_mean_length_over_time /= len(run_groups)
    summed_mean_branches_over_time /= len(run_groups)
    summed_pierce_percent_over_time /= len(run_groups)

    return (summed_mean_length_over_time, summed_mean_branches_over_time, summed_pierce_percent_over_time,
            mean_length_per_timestep, mean_branches_per_timestep, mean_pierce_percent_per_timestep)


def main(label_list, root):
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

    timesteps = np.arange(min_time, max_time + 1, timestep)
    data_list = []

    for v, folder in enumerate(folders):
        print("Starting with folder " + str(v + 1))
        data_list.append(process_folder(folder, len(unique_time_array)))

    # Process each folder and generate plots for each category
    categories = ['hyphal length', 'number of branches']
    colors = plt.colormaps["tab10"]

    for i, category in enumerate(categories):
        # Create Figure
        plt.figure(figsize=(6, 4))
        plt.style.use("seaborn-v0_8-whitegrid")

        # Plot data in figure
        for j, folder in enumerate(folders):
            mean_data_over_time = data_list[j][i]
            mean_data_per_run = data_list[j][3+i]
            color = colors(j)  # Get the color for this parameter set

            # Plot line over time
            plt.plot(timesteps, mean_data_over_time, label=label_list[j], marker=' ', markersize=5,
                     linestyle='-', color=color)
            plt.errorbar(timesteps, mean_data_over_time, yerr=np.std(mean_data_per_run, axis=1), elinewidth=1,
                         capsize=3, color=color, marker=".")

        # Create formatting for all plots
        plt.title('Mean {} over time for {} runs '.format(category, number_of_runs))
        plt.xlabel('Time in minutes')
        plt.xticks(np.arange(min_time, max_time+1, 60))
        y_label_str = category.capitalize()
        if category == "hyphal length":
            y_label_str += " in µm"
        if category == "pierce percentage":
            y_label_str = "Pierced conidia in %"
        plt.ylabel(y_label_str)
        plt.legend()

        # Save the plot as an image
        plt.savefig(
            root / '{}_over_time.png'.format(category.lower().replace(" ", "_")),
            bbox_inches='tight')
        plt.close()


if __name__ == "__main__":
    #labels = ["$b_{delay}=0 min$", "$b_{delay}=1min$", "$b_{delay}=100µm$"]
    root = Path("/home/mwank/Desktop/ExplorativeVis/confinement_ratio_circle_test/Eva,Surface,cAng")
    labels = ["Evasion", "Nothing", "RandomSurface"]
    main(labels, root)
