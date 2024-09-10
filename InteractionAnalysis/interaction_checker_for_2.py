import re
from pathlib import Path

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.lines import Line2D

from Optimized_VisualizeBoxplots import process_folder, process_group


# Helper function to extract the bDD value from folder name
def extract_bdd_value(folder_name):
    for label in x_divider:
        if label in folder_name:
            return label
    return None


# Sorting function for folders based on x_ticks value
def sort_folders(folders, x_divider_names):
    sorted_folders = sorted(folders, key=lambda f: x_divider_names.index(extract_bdd_value(str(f))))
    return sorted_folders


def main(root, x_divider, x_ticks_labels, red_blue_divider, linestyle_divider, categories, enable_legend):

    """
    # Error check
    if len(red_blue_divider) != 2:
        print("Error: red_blue divider must have length of 2")
        return None
    if len(linestyle_divider) != 3:
        print("Error: linestyle divider must have length of 3")
        return None
    if len(x_divider) != len(x_ticks_labels):
        print("Error: x_labels and  x_divider must have the same length")
        return None
    """

    # Read the CSV files
    folders = sorted([file for file in root.iterdir() if file.is_dir()])

    # Create a data dict for both groups
    red_data_dict = {}
    blue_data_dict = {}

    for v, folder in enumerate(folders):
        if all(name in str(folder) for name in red_blue_divider) and any(name in str(folder) for name in x_divider):
            print("Added folder " + str(v + 1) + " to group red and blue!")
            data = process_folder(folder)
            red_data_dict[folder] = data
            blue_data_dict[folder] = data
        elif red_blue_divider[0] in str(folder) and any(name in str(folder) for name in x_divider):
            print("Added folder " + str(v + 1) + " to group red!")
            blue_data_dict[folder] = process_folder(folder)
        elif red_blue_divider[1] in str(folder) and any(name in str(folder) for name in x_divider):
            print("Added folder " + str(v + 1) + " to group blue!")
            red_data_dict[folder] = process_folder(folder)
        else:
            print("Added folder " + str(v + 1) + " to no group!")

    line_styles = ["solid", "dashed", "dotted"]
    data_dicts = [blue_data_dict, red_data_dict]
    for i, category in enumerate(categories):
        plt.figure(figsize=(6, 4))
        plt.style.use("seaborn-v0_8-whitegrid")
        plt.xticks(range(1, len(x_ticks_labels) + 1), x_ticks_labels)

        for u, name1 in enumerate(red_blue_divider):
            # Set color and data_dict
            color = "#74B6CF"
            data_dict = data_dicts[u]
            lines = []
            for v, name2 in enumerate(linestyle_divider):
                category_data = [np.array(data_dict[folder][i]) for folder in sort_folders(data_dict, x_divider)
                                 if name1 in str(folder) and name2 in str(folder)]
                lines.append(category_data)

            for linestyle_index, line_data in enumerate(lines):
                mean_data1 = [np.mean(data) for data in line_data]
                std_data1 = [np.std(data) for data in line_data]
                plt.plot(range(1, len(x_divider) + 1), mean_data1, marker=' ', markersize=5,
                         linestyle=line_styles[linestyle_index], color=color)
                plt.errorbar(range(1, len(x_divider) + 1), mean_data1, yerr=std_data1,
                             elinewidth=1, capsize=3, color=color, marker=".", linestyle=line_styles[linestyle_index])

        # Create formatting for all plots
        plt.title('Interaction of parameters')
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
        plt.xlabel("Apical dominance range in µm")
        plt.ylim(bottom=0)
        if category == "confinement ratio":
            plt.ylim(bottom=0, top=1)

        # Creating custom legend entries for line styles
        legend_elements_styles = [
            Line2D([0], [0], color="#74B6CF", lw=2, linestyle="solid", label=r'$p_{A} = 0$'),
            Line2D([0], [0], color="#74B6CF", lw=2, linestyle="dashed", label=r'$p_{A} = 0.5$'),
            Line2D([0], [0], color="#74B6CF", lw=2, linestyle="dotted", label=r'$p_{A} = 1$')
        ]

        if enable_legend:
            # Shrink current axis's height by 10% on the bottom
            box = plt.gca().get_position()
            plt.gca().set_position([box.x0, box.y0 + box.height * 0.2,
                                    box.width, box.height * 0.8])

            # Put the first legend below current axis
            plt.legend(handles=legend_elements_styles, loc='upper center', bbox_to_anchor=(0.5, -0.1),
                       fancybox=True, shadow=True, ncol=3)

        # Save the plot as an image
        plt.savefig(
            root / '{}_parameter_interaction.png'.format(category.lower().replace(" ", "_")),
            bbox_inches='tight')
        plt.close()


if __name__ == "__main__":

    root = Path("Data_To_Vis_Interaction")
    enable_legend = False
    x_divider = ["bDelay0", "bDelay20", "bDelay40", "bDelay60", "bDelay80", "bDelay100"]
    x_ticks_labels = ["0", "20", "40", "60", "80", "100"]
    red_blue_divider = ["pA", "pA"]
    linestyle_divider = ["pA0_", "pA0.5", "pA1_"]
    categories = ['hyphal length', 'number of branches', 'branch level', 'hyphal growth unit', 'covered area',
                  "pierce percentage", "confinement ratio", "smallest circle radius", "coverage"]
    main(root, x_divider, x_ticks_labels, red_blue_divider, linestyle_divider, categories, enable_legend)
