import random
import math
import pandas as pd
import numpy as np


class Circle:
    def __init__(self, x, y, r):
        self.x = x
        self.y = y
        self.r = r


def dist(p1, p2):
    return math.hypot(p1[0] - p2[0], p1[1] - p2[1])


def is_in_circle(c, p):
    return dist((c.x, c.y), p) <= c.r


def circle_from_2_points(p1, p2):
    center_x = (p1[0] + p2[0]) / 2
    center_y = (p1[1] + p2[1]) / 2
    radius = dist(p1, p2) / 2
    return Circle(center_x, center_y, radius)


def circle_from_3_points(p1, p2, p3):
    ax, ay = p1
    bx, by = p2
    cx, cy = p3
    d = 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    if d == 0:
        return None
    ux = ((ax**2 + ay**2) * (by - cy) + (bx**2 + by**2) * (cy - ay) + (cx**2 + cy**2) * (ay - by)) / d
    uy = ((ax**2 + ay**2) * (cx - bx) + (bx**2 + by**2) * (ax - cx) + (cx**2 + cy**2) * (bx - ax)) / d
    center = (ux, uy)
    radius = dist(center, p1)
    return Circle(ux, uy, radius)


def trivial_circle(points):
    if len(points) == 0:
        return Circle(0, 0, 0)
    elif len(points) == 1:
        return Circle(points[0][0], points[0][1], 0)
    elif len(points) == 2:
        return circle_from_2_points(points[0], points[1])
    elif len(points) == 3:
        return circle_from_3_points(points[0], points[1], points[2])
    else:
        raise ValueError("trivial_circle called with more than 3 points")


def welzl(P, R):
    if len(P) == 0 or len(R) == 3:
        return trivial_circle(R)

    p = random.choice(P)
    P.remove(p)

    D = welzl(P, R)

    if is_in_circle(D, p):
        P.append(p)
        return D

    R.append(p)
    result = welzl(P, R)
    P.append(p)
    R.pop()
    return result


def smallest_enclosing_circle(points):
    P = points.copy()
    random.shuffle(P)
    return welzl(P, [])


# Define file path of Excel file
file_path = 'DataFiles/AufgereinigteDaten.xlsx'

# Specify the sheet name
sheet_name = 'Welzl'

# Read the Excel file and select the specified sheet
df = pd.read_excel(file_path, sheet_name=sheet_name)

# Group by 'Exp' (Experiment)
groups = df.groupby(['Exp'])
circle_radii = []

for experiment_name, experiment_data in groups:

    # Group by FilamentID (individual hyphae)
    hyphae_groups = experiment_data.groupby(["FilamentID"])

    for hypha_name, hypha_data in hyphae_groups:

        # Extract the points (Pt Position X, Pt Position Y, Pt Position Z)
        points = hypha_data[['Pt Position X', 'Pt Position Y']].drop_duplicates()

        # Convert to a list of tuples
        points_list = list(points.itertuples(index=False, name=None))

        # Skip trivial radius 0 circles
        if len(points_list) <= 1:
            continue

        # Calculate the smallest circle for hypha
        circle = smallest_enclosing_circle(points_list)

        # Add radius of smallest circle
        circle_radii.append(circle.r)

mean_circle_radius = np.mean(circle_radii)
deviation_circle_radius = np.std(circle_radii)
print(mean_circle_radius)
print(deviation_circle_radius)
