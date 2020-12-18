# import scorpio
import matplotlib.pyplot as plt
# from astropy.io import fits

# from matplotlib.testing.decorators import check_figures_equal

def plotear_linea(arr, ax=None):
    if ax is None:
        ax = plt.gca()
        fig = plt.gcf()
        fig.set_size_inches(10,10)

    ax.plot(arr)
    ax.set_title("Un titulo")
    
    return ax

