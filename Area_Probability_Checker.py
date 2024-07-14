import random

# Define the properties of the areas
areas = [
    {'percentage': 38.5, 'branch_level': 1, 'branch_number': 1},
    {'percentage': 29.2, 'branch_level': 2, 'branch_number': 2.5},  # Average of 2 to 3 branches
    {'percentage': 19.86, 'branch_level': 3, 'branch_number': 4},
    {'percentage': 12.34, 'branch_level': 3, 'branch_number': 5}
]

def assign_area(areas):
    """Assign an area to an object based on the area percentages."""
    rand_value = random.uniform(0, 100)
    cumulative_percentage = 0
    for area in areas:
        cumulative_percentage += area['percentage']
        if rand_value <= cumulative_percentage:
            return area
    return areas[-1]  # In case of any rounding errors

def calculate_averages(n, areas):
    total_branch_level = 0
    total_branch_number = 0

    for _ in range(n):
        area = assign_area(areas)
        total_branch_level += area['branch_level']
        total_branch_number += area['branch_number']

    avg_branch_level = total_branch_level / n
    avg_branch_number = total_branch_number / n

    return avg_branch_level, avg_branch_number

# Number of objects to place in the square
n = 200*32  # You can change this value as needed

avg_branch_level, avg_branch_number = calculate_averages(n, areas)

print(f"Average Branch Level: {avg_branch_level}")
print(f"Average Number of Branches: {avg_branch_number}")
