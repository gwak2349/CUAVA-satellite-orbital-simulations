import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np





### Code below was adapted from https://www.youtube.com/watch?v=NXnLWnTwnFI&list=PLOIRBaljOV8ghaN9XcC4ubv-QLPOQByYQ&index=2





fig, ax = plt.subplots()
# def plot_world_map():
coastal_coordinates = np.genfromtxt('coastlines.csv', delimiter=',')
plt.plot(coastal_coordinates[:,0], coastal_coordinates[:,1], 'mo', markersize=0.3, color='fuchsia')
ax.patch.set_facecolor('black')
plt.grid()
plt.show()


###inspired from https://www.youtube.com/watch?v=3c4Vk5wDIsk

import cartopy.crs as cy
import cartopy.feature as cf

plt.figure(figsize=(10,8))
ax=plt.axes(projection=cy.PlateCarree())
ax.add_feature(cf.COASTLINE, alpha=0.8)
ax.add_feature(cf.BORDERS, alpha=0.8, linestyle='--')
ax.add_feature(cf.LAND)
ax.add_feature(cf.LAKES)
ax.add_feature(cf.RIVERS)
ax.add_feature(cf.OCEAN)

ax.stock_img()

ax.gridlines(draw_labels=True, color='black', alpha=0.6, linestyle='--')


plt.show()
