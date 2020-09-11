import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.animation as anim
from IPython.display import HTML

# initialising variables
print("How many cores did you use?")
cores = int(input())
settings = np.loadtxt("settings.txt", dtype=int)
period = settings[2]
rows = settings[0]
cols = settings[1]
p_rows = settings[3]
p_cols = settings[4]

p_data = []
dimensions = []
for i in range(cores):
    # reading in the game of life grids for every core
    file_name = "output_" + str(i) + ".txt"
    p_data.append(np.loadtxt(file_name, dtype=int))
    # reading in the processor dimensions
    file_name2 = "dimensions_" + str(i) + ".txt"
    dimensions.append(np.loadtxt(file_name2, dtype=int))

whole = []
imgs = []
# combining the grids of each processor for multiple generations
for data in range(0,period):
    # initialising how big each array is
    sub_row = dimensions[0][0]
    ind_low = sub_row*data
    ind_mid = sub_row+sub_row*data
    ind_up = (sub_row*2)+sub_row*data

    initial = p_data[0][ind_low:ind_mid]
    initial2 = p_data[p_cols][ind_low:ind_mid]
    total = 0
    for j in range(p_cols-1):
        for i in range(p_rows-1):
            sub_row = dimensions[i+j][0]
            ind_low = sub_row*data
            ind_mid = sub_row+sub_row*data
            ind_up = (sub_row*2)+sub_row*data
            add = p_data[j+1][ind_mid:ind_up]
            TESTCOL = np.concatenate((initial, add), axis=1)
            initial = TESTCOL
            under = p_data[p_cols+i+j+1][ind_mid:ind_up]
            TESTROW = np.concatenate((initial2, under), axis=1)
            initial2 = TESTROW
    total = np.concatenate((initial, initial2), axis=0)
    whole.append(total)

# creating an animation
fig = plt.figure(figsize=(10,10))
for j in range(period):
    plt.title('Game of Life')
    img = plt.imshow(whole[j], cmap='gray')
    imgs.append([img])

print()
print('Building animation...')

ani = anim.ArtistAnimation(fig, imgs, interval=100, blit=True)
ani.save('life.gif', fps=60)
