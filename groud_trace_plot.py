import cartopy.crs as cy
import cartopy.feature as cf
import matplotlib.pyplot as plt

def ground_trace_plot(alpha, delta):

    ## Plotting the ground trace plot

    #utilising cartopy module as a background image for ground trace plot
    fig = plt.figure(figsize=(10,8))
    ax=plt.axes(projection=cy.PlateCarree())
    ax.add_feature(cf.COASTLINE, alpha=0.8)
    ax.add_feature(cf.BORDERS, alpha=0.8, linestyle='--')
    ax.add_feature(cf.LAND)
    ax.add_feature(cf.LAKES)
    ax.add_feature(cf.RIVERS)
    ax.add_feature(cf.OCEAN)
    ax.stock_img()
    ax.gridlines(draw_labels=True, color='black', alpha=0.6, linestyle='--')

    po = [] #longitude vector
    ho = [] #latitude vector

    last_i = 0 #counter

## The for loop below matches the corresponding longitude with its latitude
    for i in range(len(alpha)-1):    
        if alpha[i+1]>alpha[i]:
            po = alpha[last_i:i]
            ho = delta[last_i:i]
            plt.plot(po, ho, color='black')
            po = []
            ho = []
            last_i = i+1


    po = alpha[last_i:-1]
    ho = delta[last_i:-1]

    plt.plot(po,ho,'black')
    plt.show()