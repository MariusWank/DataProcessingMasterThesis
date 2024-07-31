import matplotlib.pyplot as plt
import numpy as np
import random
import cv2
import os


# Helper function to compute the distance between two points
def distance(p1, p2):
    return np.sqrt((p1[0] - p2[0]) ** 2 + (p1[1] - p2[1]) ** 2)


# Function to check if a point is inside a circle
def is_in_circle(circle, p):
    center, radius = circle
    return distance(center, p) <= radius


# Function to compute the smallest circle from up to three points
def trivial_circle(R):
    if len(R) == 0:
        return (0, 0), 0
    elif len(R) == 1:
        return R[0], 0
    elif len(R) == 2:
        center = ((R[0][0] + R[1][0]) / 2, (R[0][1] + R[1][1]) / 2)
        radius = distance(R[0], center)
        return center, radius
    else:
        A, B, C = R
        ox = (min(A[0], B[0], C[0]) + max(A[0], B[0], C[0])) / 2
        oy = (min(A[1], B[1], C[1]) + max(A[1], B[1], C[1])) / 2
        ax, ay = A[0] - ox, A[1] - oy
        bx, by = B[0] - ox, B[1] - oy
        cx, cy = C[0] - ox, C[1] - oy
        d = (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by)) * 2
        if d == 0:
            return (ox, oy), max(distance((ox, oy), A), distance((ox, oy), B), distance((ox, oy), C))
        x = ox + ((ax * ax + ay * ay) * (by - cy) + (bx * bx + by * by) * (cy - ay) + (cx * cx + cy * cy) * (
                    ay - by)) / d
        y = oy + ((ax * ax + ay * ay) * (cx - bx) + (bx * bx + by * by) * (ax - cx) + (cx * cx + cy * cy) * (
                    bx - ax)) / d
        center = (x, y)
        radius = distance(center, A)
        return center, radius


# Welzl's algorithm
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


# Function to plot the current state of the algorithm and save it as an image
def plot_welzl(P, R, circle, step):
    plt.figure(figsize=(8, 8))
    if P:
        plt.scatter(*zip(*P), color='blue', label='Points')
    if R:
        plt.scatter(*zip(*R), color='red', label='Boundary Points')

    if circle is not None:
        center, radius = circle
        circle_plot = plt.Circle(center, radius, color='green', fill=False, linestyle='--', linewidth=2)
        plt.gca().add_patch(circle_plot)
        plt.scatter(*center, color='green', marker='x', label='Circle Center')

    plt.title(f'Step {step}')
    plt.legend()
    plt.gca().set_aspect('equal', adjustable='box')
    plt.xlim(-10, 10)
    plt.ylim(-10, 10)
    plt.savefig(f'frames/frame_{step:03d}.png')
    plt.close()


# Generating a set of random points and ensure three specific points are included
np.random.seed(0)
points = np.random.uniform(-2, 2, (7, 2)).tolist()

# Ensure these three points form a non-collinear triangle outside the other points
special_points = [(8, 0), (-4, 7), (-4, -7)]
points.extend(special_points)

# Running the Welzl algorithm and plotting each step
R = []
circle = None
step = 0


def welzl_with_plot(P, R, step):
    global circle
    plot_welzl(P, R, circle, step)

    if len(P) == 0 or len(R) == 3:
        circle = trivial_circle(R)
        plot_welzl(P, R, circle, step + 1)
        return circle

    p = random.choice(P)
    P.remove(p)

    D = welzl_with_plot(P, R, step + 1)

    if is_in_circle(D, p):
        P.append(p)
        circle = D
        plot_welzl(P, R, circle, step + 1)
        return D

    R.append(p)
    result = welzl_with_plot(P, R, step + 1)
    P.append(p)
    R.pop()
    circle = result
    plot_welzl(P, R, circle, step + 1)
    return result


# Ensure frames directory exists
os.makedirs('frames', exist_ok=True)

smallest_circle = welzl_with_plot(points, R, step)

print(smallest_circle)

# Create a video from the saved images
image_files = [f'frames/frame_{i:03d}.png' for i in range(len(os.listdir('frames')))]
frame = cv2.imread(image_files[0])
height, width, layers = frame.shape
video = cv2.VideoWriter('welzl_algorithm.avi', cv2.VideoWriter_fourcc(*'DIVX'), 2, (width, height))  # 2 fps

# Duplicate each frame 10 times to slow down the video
for image_file in image_files:
    frame = cv2.imread(image_file)
    for _ in range(2):  # Adjust this value to change the duration of each frame
        video.write(frame)

cv2.destroyAllWindows()
video.release()

# Cleaning up the images
for image_file in image_files:
    os.remove(image_file)
os.rmdir('frames')
