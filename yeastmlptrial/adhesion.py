#Import packages
import numpy as np
import matplotlib.pyplot as plt
import statistics as st
import scipy.stats as stats
import scipy
import pandas as pd
from skimage.morphology import square
from skimage.measure import label, regionprops
from matplotlib.collections import PatchCollection
import matplotlib.patches as mpatches
from matplotlib.axes import Axes
from skimage.filters import threshold_otsu
import skimage.io
import skimage.filters
import skimage.measure
import seaborn as sns
import os
import glob
import re
import warnings


#Show the actual image from an image matrix
def show_image(image, title='Image', cmap_type='viridis'):
    """Displays an image.
    image: numpy array, representing an image.
    title: string
    cmap_type: string, the type of colormap to be used in the display,
               for example 'gray' or 'jet'.
    """
    
    plt.imshow(image, cmap=cmap_type)
    plt.title(title)
    plt.axis('off')


def crop(image, x1=60, x2=60, y1=60, y2=60):
    """Crop an image.
    image: numpy array, representing an image.
    x1: int, number of pixels to remove from the left side.
    x2: int, number of pixels to remove from the right side.
    y1: int, number of pixels to remove from the top.
    y2: int, number of pixels to remove from the bottom."""


    # Remove x1 pixels from the left side
    image = np.delete(image, range(x1), axis=1)
    # Remove x2 pixels from the right side
    image = np.delete(image, range(image.shape[1]-x2, image.shape[1]), axis=1)
    # Remove y1 pixels from the top
    image = np.delete(image, range(y1), axis=0)
    # Remove y2 pixels from the bottom
    image = np.delete(image, range(image.shape[0]-y2, image.shape[0]), axis=0)
    return image

#Define each colony as a distinct "object"
def connected_components(filename, connectivity=2,x1=60,x2=60,y1=60,y2=60):
    """Find and label objects in an image.
    filename: string, the name of an image file.
    connectivity: int, the connectivity of the objects, either 1 or 2.
    x1: int, number of pixels to remove from the left side, goes into crop().
    x2: int, number of pixels to remove from the right side, goes into crop().
    y1: int, number of pixels to remove from the top, goes into crop().
    y2: int, number of pixels to remove from the bottom, goes into crop().
    
    Returns: list of region properties and the mask of the image."""

    # load the image
    image = plt.imread(filename)
    #Crop the image
    image=crop(image,x1=x1,x2=x2,y1=y1,y2=y2)
    #blur the image
    image = skimage.filters.gaussian(image, sigma=7, truncate=8, channel_axis=True)
    # make an inverted binary
    thresh = threshold_otsu(image)
    binary = image > thresh
    binary_mask=np.invert(binary)

    mask = skimage.morphology.remove_small_objects(binary_mask,1000)

    labeled_image, count = skimage.measure.label(mask,connectivity=connectivity, return_num=True)

    #regions group the individual object properties
    regions=regionprops(labeled_image)
    return regions, mask


#Name each object after its position (A1-H12)
def create_raw_layout(regions,mask):
    """Create a raw layout of the plate.
    regions: list of region properties, as returned by connected_components().
    mask: numpy array, the mask returned by connected_components().
    
    Returns: The raw layout of the plate."""


    x_array=[]
    y_array=[]
    for r in regions:
        x,y=r["centroid"]
        x_array.append(x)
        y_array.append(y)
    rows=pd.cut(x_array,8, labels=["A","B","C","D","E","F","G","H"])
    cols=pd.cut(y_array,12, labels=["1","2","3","4","5","6","7","8","9","10","11","12"])
    index=[]
    for i in range(0,len(rows)):
        index.append(rows[i]+cols[i])
    d = {'x': x_array, 'y': y_array, "row": rows, "col":cols, "index": index}
    df = pd.DataFrame(data=d)


    m=np.zeros((len(mask),len(mask[0])))

    for ind in df["index"]:
        i=df[df['index']==ind].index[0]
        r=df[df["row"]==df.iat[i,2]]
        c=df[df["col"]==df.iat[i,3]]
        x=st.mean([df.iat[r[r["col"]=="12"].index[0],0],df.iat[r[r["col"]=="1"].index[0],0]])
        y=st.mean([df.iat[c[c["row"]=="A"].index[0],1],df.iat[c[c["row"]=="H"].index[0],1]])
        mean_region_area = st.mean([r.area for r in regions])*1.1
        sqrt_area = int(np.sqrt(mean_region_area))
        minr=int(round(x-sqrt_area/2))
        minc=int(round(y-sqrt_area/2))
        maxr=int(round(x+sqrt_area/2))
        maxc=int(round(y+sqrt_area/2))
        m[minr:maxr,minc:maxc]=square(sqrt_area)

    #Find and label objects in the new binary
    labeled_image, count = skimage.measure.label(m,connectivity=2, return_num=True)
    return labeled_image

def im_to_raw_layout(image_fname,x1=60,x2=60,y1=60,y2=60):
    """Create a raw layout from an image. Chains together
      connected_components() and create_raw_layout().
    image_fname: string, the name of an image file.
    x1: int, number of pixels to remove from the left side, goes into crop().
    x2: int, number of pixels to remove from the right side, goes into crop().
    y1: int, number of pixels to remove from the top, goes into crop().
    y2: int, number of pixels to remove from the bottom, goes into crop().

    Returns: The raw layout of the plate.
    """
    regions,mask=connected_components(image_fname, x1=x1,x2=x2,y1=y1,y2=y2)
    labeled_image=create_raw_layout(regions,mask)
    return labeled_image

def pad_layout(image, raw_layout):
    """Pad the raw_layout with zeros to make it the same shape as the image.
    image: numpy array, representing an image.
    raw_layout: numpy array, representing a layout.

    Returns: The padded layout.
    """
    # Get the difference in shape
    x_diff = image.shape[0] - raw_layout.shape[0]
    y_diff = image.shape[1] - raw_layout.shape[1]
    # Pad the layout with zeros
    layout = np.pad(raw_layout, (int(x_diff/2), int(y_diff/2)), 'constant')
    
    
    return layout

def file_to_initial_layout(image_for_raw_layout,target_image,x1=60,x2=60,y1=60,y2=60):
    """Create an initial layout from an image. Chains together
        im_to_raw_layout() and pad_layout().
    image_for_raw_layout: string, the name of an image file.
    target_image: numpy array, representing an image.
    x1: int, number of pixels to remove from the left side, goes into crop().
    x2: int, number of pixels to remove from the right side, goes into crop().
    y1: int, number of pixels to remove from the top, goes into crop().
    y2: int, number of pixels to remove from the bottom, goes into crop().

    Returns: The initial layout.
    """

    raw_layout = im_to_raw_layout(image_for_raw_layout,x1=x1,x2=x2,y1=y1,y2=y2)
    initial_layout = pad_layout(target_image, raw_layout)
    return initial_layout

def calculate_sums(image, layout):
    """Calculate the sum of squares for each square in the layout.
    This will be used to optimize the layout.
    image: numpy array, representing an image.
    layout: numpy array, representing a layout.

    Returns: A list of the sum of squares for each square in the layout."""

    
    #make binary from image
    image = image > threshold_otsu(image)
    # Calculate the sum of each square
    sums = np.sum(np.square(image[layout != 0]))  
    
    return sums

def optimize_layout(image, layout):
    """
    Optimize the layout by testing different offsets and finding the one with the lowest sum of squares.
    image: numpy array, representing an image.
    layout: numpy array, representing a layout.

    Returns: The first optimized layout by roughly fitting the layout to the image.
    """


    # Define the search space for the offsets
    x_offsets = range(-250, 250, 12)
    y_offsets = range(-250, 250, 12)
    
    # Initialize the minimum sum of squares and the optimal layout
    min_sum_of_squares = float('inf')
    optimal_layout = None
    
    # Iterate over all possible offsets and find the one with the lowest sum of squares
    for x_offset in x_offsets:
        for y_offset in y_offsets:
            # Shift the layout by the current offset
            shifted_layout = np.roll(layout, (y_offset, x_offset), axis=(0, 1))
            
            # Calculate the sum of squares between the image and the shifted layout
            
            sum_of_squares = calculate_sums(image, shifted_layout)
            
            
            # Update the minimum sum of squares and the optimal layout if necessary
            if sum_of_squares < min_sum_of_squares:
                min_sum_of_squares = sum_of_squares
                optimal_layout = shifted_layout
    
    # Return the optimal layout and the minimum sum of squares
    return min_sum_of_squares, optimal_layout
    

def optimize_layout_2(image,optimal_layout):
    """Optimize the layout of the image.
    image: numpy array, representing an image.
    optimal_layout: numpy array, representing a layout, output from optimize_layout().

    Returns: The final optimized layout.
    """
    
    # Define the size of the image
    height, width = image.shape
    #max offset cannot be more so long that the layout is outside the image, using layout_original_dims
    

    # Define the initial guess for the size of each square and border
    x0 = [0,0]
    # Define the function to minimize
    
    
    # Define the bounds for the optimization
    bounds = ( (-15, 15), (-15, 15) )
    def objective(x):
        #calculate max offset as the difference between the size of the image and the size of the layout
        new_layout = optimal_layout.copy()
        new_layout = np.roll(new_layout, int(x[0]), axis=1)
        new_layout = np.roll(new_layout, int(x[1]), axis=0)
        sums = calculate_sums(image, new_layout)
        
        # Return the sum of the sums
        return sums
    # Perform the optimization using basin hopping
    result = scipy.optimize.basinhopping(objective, x0, niter=100, minimizer_kwargs={"method":"L-BFGS-B", "bounds":bounds},stepsize=1)


    
    
    final_layout = optimal_layout.copy()
    final_layout = np.roll(final_layout, int(result.x[0]), axis=1)
    final_layout = np.roll(final_layout, int(result.x[1]), axis=0)
    return final_layout, result.x

def show_layout_over_im(image_before, image_after,initial_layout, final_layout):
    """Show the layouts over the before after images.
    image_before: numpy array, representing an image.
    image_after: numpy array, representing an image.
    initial_layout: numpy array, representing an initial layout.
    final_layout: numpy array, representing a final layout.

    Returns: None, just prints the images.
    """

    # Create a new figure
    fig, ax = plt.subplots(1, 3, figsize=(15, 15))
    
    # Show the image twice side by side
    ax[0].imshow(image_before, cmap='gray')
    ax[1].imshow(image_before, cmap='gray')
    ax[2].imshow(image_after, cmap='gray')
    
    
    # Show the initial layout with alpha=0.3
    ax[0].imshow(initial_layout, alpha=0.35)
    ax[1].imshow(final_layout, alpha=0.35)
    ax[2].imshow(final_layout, alpha=0.35)


def layout_optimization(image_before, image_after, initial_layout,show_layout_fits=False):
    """Optimize the layout of the image. Chains together optimize_layout() and
    optimize_layout_2() with the option to show the layout fits.

    image_before: numpy array, representing an image.
    image_after: numpy array, representing an image.
    initial_layout: numpy array, representing an initial layout.
    show_layout_fits: bool, whether to show the layout fits.

    Returns: The final optimized layout.
    """

    min_sum_of_squares, optimal_layout = optimize_layout(image_before, initial_layout)
    final_layout, result = optimize_layout_2(image_before,optimal_layout)
    if show_layout_fits:

        show_layout_over_im(image_before, image_after, initial_layout, final_layout)
    return final_layout, result


def create_strain_layout(map_filename):
    """Create a list of strains in the layout.
    map_filename: str, the filename of the map.

    Returns: A list of strains in the layout."""

    try:
        map_df = pd.read_excel(map_filename,header=None)
    except:
        map_df = pd.read_csv(map_filename,header=None)
    # Transpose the dataframe to get the strains in each of the 8 columns
    strains_df = pd.DataFrame.transpose(map_df)
    # Create a list of strains in the layout
    strain_layout = []
    for i in range(0, 8):
        strain_layout.extend(strains_df[i])
    return strain_layout


def mean_int(image,labeled_image, mean_or_median="mean"):
    """Calculate the mean intensity of each square of the image."""


    mi=[]
    for i in range(96):
        #take the plate value except the first letter
        if mean_or_median == "mean":
            m=np.mean(image[np.where(labeled_image==i+1)])
        elif mean_or_median == "median":
            m=np.median(image[np.where(labeled_image==i+1)])
        else:
            raise ValueError("mean_or_median must be mean or median")
        mi.append(m)
    maxval=np.max(image)
    return mi, maxval 

def create_plate():
    """Create a list of positions in the layout."""
    plate=[]
    for i in range(8):
        for j in range(12):
            plate.append(chr(i+65)+str(j+1))
    return plate


def processing(image, image_w, plate, strain_layout, labeled_image, keep_all=False, mean_or_median="mean"):
    """Take a before and after image and compare the mean intensity of each square.
    image: numpy array, representing an image.
    image_w: numpy array, representing an image.
    plate: list, representing the positions in the layout.
    strain_layout: list, representing the strains in the layout.
    labeled_image: numpy array, representing the labeled image.
    keep_all: bool, whether to keep all the values.
    mean_or_median: str, whether to use the mean or median intensity.

    Returns: A dataframe with the position, strain, before wash, after wash, and ratio.
    """
    
    mean_intensities, max_val = mean_int(image, labeled_image, mean_or_median)
    mean_intensities_w, max_val_w = mean_int(image_w, labeled_image, mean_or_median)

    before_vals = mean_intensities / max_val
    before_vals = before_vals - np.min(before_vals)

    after_vals = mean_intensities_w / max_val_w
    # after_vals- np.min(after_vals or if less than 0, then 0
    after_vals = np.maximum(after_vals - np.min(after_vals), 0)

    d = {
        "position": plate,
        "strain": strain_layout,
        "before_wash": before_vals,
        "after_wash": after_vals,
    }
    df = pd.DataFrame(data=d)
    df["ratio"] = df["after_wash"] / df["before_wash"]
    # remove entries where position starts with A or H,
    if keep_all == False:
        df = df[~df["position"].str.contains("A|H")]
        # for position, put a _ after the letter and the number like "B_1_"
        
        
        warnings.filterwarnings("ignore")
        df["position"] = df["position"].str.replace("([A-Z])([0-9])", r"\1_\2")
        # add a _ to the end of position entries
        df["position"] = df["position"] + "_"
        # remove entries where position has format letter + 1_end, letter + 12
        df = df[~df["position"].str.contains("_1_|_12_")]
    else:
        warnings.filterwarnings("ignore")
        df["position"] = df["position"].str.replace("([A-Z])([0-9])", r"\1_\2")
        df["position"] = df["position"] + "_"

    return df

def get_filenames(folder):
    """Get the filenames of the images in a folder.
    folder: str, the folder containing the images.

    Returns: A list of filenames.
    """
    filenames = []
    for filename in os.listdir(folder):
        #if ends with jpg or png
        if filename.endswith(".jpg") or filename.endswith(".png"):
            filenames.append(os.path.join(folder, filename))
    # find last integer in name using regex
    index = re.compile(r"\d+")
    # sort filenames by the last integer in the name
    filenames.sort(key=lambda f: int(index.findall(f)[-1]))
    return filenames


#create a function that takes two lists of before and after wash images and returns a list of dataframes


def result_from_of_ims(image_before,image_after,strain_layout_file, image_for_raw_layout,show_layout_fits = False,x1=60,x2=60,y1=60,y2=60, keep_all=False,mean_or_median = "mean"):
    """Take a before and after image and return a dataframe
    image_before: str, the filename of the before image.
    image_after: str, the filename of the after image.
    strain_layout_file: str, the filename of the strain layout.
    image_for_raw_layout: str, the filename of the image used to create the initial layout.
    show_layout_fits: bool, whether to show the layout fits.
    x1: int, the number of pixels to crop from the left.
    x2: int, the number of pixels to crop from the right.
    y1: int, the number of pixels to crop from the top.
    y2: int, the number of pixels to crop from the bottom.
    keep_all: bool, whether to keep all the values.
    mean_or_median: str, whether to use the mean or median intensity.

    Returns: A dataframe with the position, strain, before wash, after wash, and ratio.


    """

    strain_layout = create_strain_layout(strain_layout_file)
    plate = create_plate()
    
    #initial layout will be shared between images, as we assume that the image sizes and square sizes are the same
    initial_layout = file_to_initial_layout(image_for_raw_layout, plt.imread(image_before),x1=x1,
                                            x2=x2,y1=y1,y2=y2)

   
    before=plt.imread(image_before)
    after=plt.imread(image_after)
       

        #Assume that before and after pictures have the same layout: no drastic change in position of plate.
    final_layout, result = layout_optimization(before, after, initial_layout,show_layout_fits = show_layout_fits)
        
        #invert images
    before=np.invert(before)
    after=np.invert(after)
    df=processing(before,after,plate,strain_layout,final_layout, keep_all=keep_all,mean_or_median=mean_or_median)
    
    return df


def result_from_folders_of_ims(images_before_folder,images_after_folder,strain_layout_file, image_for_raw_layout,show_layout_fits = True,
                               x1=60,x2=60,y1=60,y2=60, keep_all=False,mean_or_median = "mean",
                               filtered_ratios=True,filtered_strains=True):
    """Take a folder of before and after images and return a dataframe
    of ratios with summary statistics for each strain.
    images_before_folder: str, the folder containing the before images.
    images_after_folder: str, the folder containing the after images.
    strain_layout_file: str, the filename of the strain layout.
    image_for_raw_layout: str, the filename of the image used to create the initial layout.
    show_layout_fits: bool, whether to show the layout fits.
    x1: int, the number of pixels to crop from the left.
    x2: int, the number of pixels to crop from the right.
    y1: int, the number of pixels to crop from the top.
    y2: int, the number of pixels to crop from the bottom.
    keep_all: bool, whether to keep all the values.
    mean_or_median: str, whether to use the mean or median intensity.
    filtered_ratios: bool, whether to filter the ratios.
    filtered_strains: bool, whether to filter the strains.

    Returns: A dataframe with the position, strain, before wash, after wash, and ratio.
    """

    images_before = get_filenames(images_before_folder)
    images_after = get_filenames(images_after_folder)

    strain_layout = create_strain_layout(strain_layout_file)
    plate = create_plate()
    
    #initial layout will be shared between images, as we assume that the image sizes and square sizes are the same
    initial_layout = file_to_initial_layout(image_for_raw_layout, plt.imread(images_before[0]),x1=x1,
                                            x2=x2,y1=y1,y2=y2)

    for i in range(len(images_before)):
        before=plt.imread(images_before[i])
        after=plt.imread(images_after[i])
       

        #Assume that before and after pictures have the same layout: no drastic change in position of plate.
        final_layout, result = layout_optimization(before, after, initial_layout,show_layout_fits = show_layout_fits)
        
        #invert images
        before=np.invert(before)
        after=np.invert(after)
        df=processing(before,after,plate,strain_layout,final_layout, keep_all=keep_all,mean_or_median=mean_or_median)
        
        if i==0:
            dfs=[df]
        else:
            dfs.append(df)

        

        merged_df = pd.concat(dfs)
        merged_df_pre_filter= merged_df.copy()
        # Set any ratio values greater than 1 to 1
        merged_df.loc[merged_df["ratio"] > 1, "ratio"] = 1
        
        # Remove entries with before_wash less than 0.1
        if filtered_ratios:
            merged_df = merged_df[merged_df["before_wash"] > 0.1]
        
        # Group the dataframe by strain and calculate the mean and SEM of the ratio and before_wash columns
        grouped_df = merged_df.groupby("strain", as_index=False).agg(
            {"ratio": "mean", "before_wash": "mean", "after_wash": "mean"}
        )
        grouped_df["sem"] = merged_df.groupby("strain", as_index=False)["ratio"].sem()[
            "ratio"
        ]
        # Remove entries with before_wash less than 0.1 and strain equal to "TYPO"
        if filtered_strains:
            grouped_df = grouped_df[grouped_df["before_wash"] > 0.1]
        grouped_df = grouped_df[grouped_df["strain"] != "TYPO"]
        grouped_df = grouped_df[grouped_df["strain"] != "EMPTY"]
        #from merged df take every single ratio value for a given strain and put it in a list
        if filtered_ratios:
            grouped_df["filtered_all_ratio_vals"] = grouped_df["strain"].apply(lambda x: merged_df[merged_df["strain"] == x]["ratio"].tolist())
            grouped_df["n"] = grouped_df["strain"].apply(lambda x: len(merged_df[merged_df["strain"] == x]["ratio"].tolist()))
        else:
            grouped_df["unfiltered_all_ratio_vals"] = grouped_df["strain"].apply(lambda x: merged_df_pre_filter[merged_df_pre_filter["strain"] == x]["ratio"].tolist())

    return grouped_df

def barchart_from_res(df):
    """Take a dataframe of ratios and plot a bar chart of the ratios.
    df: dataframe, the dataframe of ratios.

    Returns: A bar chart of the ratios.
    """
    plt.figure(figsize=(20,10))
    plt.bar(df["strain"],df["ratio"],yerr=df["sem"],capsize=5)
    plt.xticks(rotation=90)
    plt.ylabel("Ratio of after/before wash intensity")
    plt.title("Ratio of after/before wash intensity for each strain")
    plt.show()