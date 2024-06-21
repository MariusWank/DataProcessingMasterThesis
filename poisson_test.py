import random
import numpy as np

rate = 0.0002130556191
num = 50
events = np.zeros(num)

for t in range(0, 600):
    for i, hyphae in enumerate(range(0, num)):
        if random.random() < rate:
            events[i] += 1

percent_breached = np.count_nonzero(events) / events.size

print(percent_breached)