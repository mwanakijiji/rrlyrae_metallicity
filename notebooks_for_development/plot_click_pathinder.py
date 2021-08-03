import time
import numpy as np
import matplotlib.pyplot as plt


def tellme(s):
    print(s)
    plt.title(s, fontsize=16)
    plt.draw()

plt.clf()
#plt.setp(plt.gca(), autoscale_on=False)

# plot junk data
x_test = np.arange(100)
y_test = np.random.rand(100)
plt.plot(x_test,y_test)

num_rays = input("Enter number of cosmic rays: ")
tellme('Click to begin')

#plt.waitforbuttonpress()

row_index = 0

for n_ray in range(0,num_rays):
    pts = []
    while len(pts) < 2:
        tellme('Select 2 limits with mouse')
        pts = np.asarray(plt.ginput(2, timeout=-1))
        if len(pts) < 2:
            tellme('Too few points, starting over')
            time.sleep(1)  # Wait a second

    ph = plt.fill(pts[:, 0], pts[:, 1], 'r', lw=2)

    tellme('Happy? KEY CLICK = yes and done with spectrum; MOUSE CLICK = no; ')

    print(pts[:,0])
    # append to dataframe

    if plt.waitforbuttonpress():
      print("button press")
      break

    # Get rid of fill
    for p in ph:
        p.remove()
