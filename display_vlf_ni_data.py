from nilearn.plotting import plot_anat, plot_img, plot_stat_map
from nibabel.viewers import OrthoSlicer3D
def plot_anatomy_nifti(im):   #  Needs a nifti
    plot_anat(im, colorbar=True, cbar_tick_format="%i")

def plot_anatomy_raw(im):
    OrthoSlicer3D(im).show()
