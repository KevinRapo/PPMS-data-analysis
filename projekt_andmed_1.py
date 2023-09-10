import os
import statistics as stat
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
import matplotlib.colors as mcolors
from tkinter import filedialog
from scipy.interpolate import interp1d
from scipy.signal import argrelextrema


USER_PATH = os.getcwd()


#-------------- OPENING THE FILE AND INDEXING IT -------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
root = tk.Tk()
root.wm_attributes('-topmost', 1)
root.withdraw()

file_path = filedialog.askopenfilename()
save_to_path = os.path.dirname(file_path)
#Opens the selected data file
with open(file_path, 'r') as f:
    i = 0 #Rowindex, how many rows will get_header read in
    
    for line in f:
        i += 1
        if '[Data]' in line:
            break

#--------------READING IN THE HEADER DATA-----------------------------------------------------------------------------------


def get_header(dataFile):
    #reads in i amount of rows from the data file, i indicating the number
    #of lines from the header to the data section
    header = pd.read_csv(dataFile, nrows = i , index_col = 2, encoding = 'ANSI',names = ['1','2','3','4'], on_bad_lines = 'skip')
    
    return header

#Lines 0-i from the file
data_header = get_header(file_path)


def data_type(line1):
    #data_type checks from the first line under the header whether the file is VSM or ACMS and returns a token for further use
    x = line1['1'][1]
    
    if "VSM" in x:
        print("\nThis is a VSM data file \n")
        VSM_token = "VSM"
        return VSM_token
    
    elif "ACMS" in x:
        print("\nThis is an ACMS data file \n")
        ACMS_token = "ACMS"
        
        return ACMS_token

#Prints out the type of data and saves the corresponding token
token = data_type(data_header)

#Sample mass from data header
SAMPLE_MASS = data_header['2']['SAMPLE_MASS']
SAMPLE_MASS_NAME = "SAMPLE_MASS"

#Sample size from data header
SAMPLE_SIZE = data_header['2']['SAMPLE_SIZE']
SAMPLE_SIZE_NAME = "SAMPLE_SIZE"

def extract_float_with_unit(string):
    #extract_float_with_unit function uses regex to index the input string
    #into a (float, unit) format if it has units, if no units (float, None) format, if no value (None)
    regex = r'^([\d.]+)\s*([a-zA-Z]+(\d)?)?$'# regular expression to match float value and unit
    match = re.search(regex, string)
    
    if match:
        float_str = match.group(1)
        unit = match.group(2)
        float_val = float(float_str)
        return (float_val, unit)
    else:
        return None

def header_value_check(info, name):
    #header_value_check checks the value of interest, is the value nan or is it just a string in 
    #the wrong format, otherwise outputs the float value
    print(f"Checking:{name}, {info}")#Hetkel jääb
    
    if not isinstance(info, str) and np.isnan(info): #Checks whether the value is not a string and if it is a nan value
        print(f"NO VALID VALUE FOR {name}, value is NAN \n")
        return None
    match = extract_float_with_unit(info)
    
    if match is None: #condition for extract_float_with_unit when it didn't find a match for a desired format therefore being a string
        print(f"{name} had an input put it is not in the valid format: {info} \n")
        return None
    float_val = match[0]
    #print("float value:", float_val)#Hetkel jääb
    unit = match[1]
    #print("units:", unit)#Hetkel jääb
    
    return float_val, unit #Ühikutega peaks veel midagi tegema

sample_mass_with_unit = header_value_check(SAMPLE_MASS, SAMPLE_MASS_NAME)

def get_mass(mass):
    #If the mass > 1, indicating that it's input is in mg, divides it by a 1000 to get g
    if mass is None:
        return None
    val = mass[0] #Mass value
    if val > 1:
        val= val/1000
        print(f'Sample mass is {val:.5f} grams \n')
        
    return val

sample_mass = get_mass(sample_mass_with_unit)

header_size = header_value_check(SAMPLE_SIZE, SAMPLE_SIZE_NAME)

#Prints out the size if it is valid/ HETKEL AINULT ÜKS CASE =mm2
def get_size(size):
    if size is None:
        return None
    val = size[0]
    unit = size[1]
    if unit == "mm2":
        val = val/100
        unit = "cm2"
        print(f'Sample size is {val:.5f} {unit} \n')
        
        return val

size = get_size(header_size)

def get_thickness(data):
    #Checks whether the title contains sample thickness in nm units: e.g. "25nm" and outputs 25
    thickness = data.iloc[3,1] #Title index in table
    pattern = r"(\d+)\s*(nm)"
    match = re.search(pattern,thickness)
    if match:
        float_str = match.group(1)
        print(f"Checking thickness: {float_str} nm")
        #unit = match.group(2)
        float_val = float(float_str)*10**-7 # nm to cm conversion
        print(f"Sample thickness is: {float_val} cm \n" )
        return float_val
    else:
        print("Sample thickness not found in title \n")

thickness = get_thickness(data_header)

#----------------------READING IN AND PREPROCESSING THE MEASUREMENT DATA-------------------------------------------------

def get_measurements(dataFile):
    #get_measurements function reads in the rest of the data under the [Data] line and formats it as a table
    measurements = pd.read_csv(dataFile, skiprows = i, header = 0, encoding = "ANSI")
    return measurements

#Forms a table from the measurement data and excludes all lines that contain 2 in the "Transport Action" column
MEASUREMENT_TABLE = get_measurements(file_path)
MEASUREMENT_TABLE = MEASUREMENT_TABLE[MEASUREMENT_TABLE["Transport Action"] != 2]
MEASUREMENT_TABLE = MEASUREMENT_TABLE[["Time Stamp (sec)", "Temperature (K)", "Magnetic Field (Oe)", "Moment (emu)"]] #HETKEL ACMS FAILIGA EI TÖÖTA
MEASUREMENT_TABLE = MEASUREMENT_TABLE.dropna()
MEASUREMENT_TABLE = MEASUREMENT_TABLE.reset_index(drop=True)

def get_temp(data):
    temp = data['Temperature (K)']
    return temp

def get_H(data):
    field = data['Magnetic Field (Oe)']
    return field

def get_M(data, token):
    if token == "VSM":
        moment = data['Moment (emu)']
        return moment
    moment = data["DC Moment (emu)"]
    return moment

TEMPERATURE = get_temp(MEASUREMENT_TABLE)
MAGNETIC_FIELD = get_H(MEASUREMENT_TABLE)
MOMENT = get_M(MEASUREMENT_TABLE, token)


# KUNI SIIANI SAMA MIS EELMISES VERSIOONIS


#--------------------------------------------- MEASUREMENT ANALYSIS----------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

#-----------------------------------------------MvsH functions-------------------------------------
#-----------------This part finds the constant temperatures from the data (MvsH)-------------------

def min_max_range(data, step = 0.15):
    #min and max value for temperature range
    #divides the temperature range into intervals for temperature distribution in the data
    start = min(data)
    stop = max(data)+0.5
    # start = 1.7
    # stop = 1000
    #step = 0.15 #0.2 oli error selle failiga DataToTest/8_DyMn0p2Ga3_VSM_DC_HYS_300K_50K_2K_26102022, kaks kõrvuti vahemikku ületasid nivoo hiljem
    values = []
    current_value = start
    
    while current_value < stop:
        values.append(current_value)
        current_value += step
        
    return values

#Hetkel tegin väikese muutusega lihtsalt uue, et see ilusaks tegemine tuleb niikuinii, aga for both on tsükli jaoks kus mõlemad jne
def mood_temp_for_both(intervals):
    #determines at which temperatures the MvsH measurement was made
    #outputs the interval which satisfies 20*mean
    mean = stat.mean(intervals)
    const_temp = []
    print("Mean temp:", mean)
    for index, value in intervals.items():
        if value > 20*mean: # MINGI ILUSAM TINGIMUS!!
            print("Üle mean:", value)
            const_temp.append(index)
            
    return const_temp

def mood_temp_for_MvsH(intervals):
    #determines at which temperatures the MvsH measurement was made
    #outputs the interval which satisfies mean
    mean = stat.mean(intervals)
    const_temp = []
    print("Mean temp:", mean)
    for index, value in intervals.items():
        if value > mean: # MINGI ILUSAM TINGIMUS!!
            print("Üle mean:", value)
            const_temp.append(index)
            print("Millised temp vahemikud:", index)
            
    return const_temp


def get_const_temp(const):
    #outputs the integer values for the constant temperatures
    const_temp = []
    x = 0
    for val in const:
        
        if int(val.left) == int(val.right): # if the temp is low the interval might shift too much to group accurately so instead we round the left value to its first digit

            rounded = round(val.left , 1)
            const_temp.append(rounded)
            print(const_temp[x],"K mõõdetud väljasõltuvust")
            x += 1
            
        else:
            const_temp.append(int(val.right))
            print(const_temp[x],"K mõõdetud väljasõltuvust")
            x += 1
        
    return const_temp

#SIIANI ON MvsH UNIVERSAALNE, SEDA OSA KASUTAB KÕIKIDEL JUHTUDEL

#------------Filtering const temp points from the data for MvsH------

def get_measurement_MvsH(const_T_values, bound = 0.05): #!!!
    #Saves all the indices of the points that fall between the bound
    table = MEASUREMENT_TABLE['Temperature (K)']
    filtered_dfs = []
    all_indices = []
    
    for value in const_T_values:
        lower = value - bound
        upper = value + bound
        filtered_df = table[(table >= lower) & (table <= upper)]
        indices = filtered_df.index.tolist()
        filtered_dfs.append(filtered_df)
        all_indices.append(indices)
        
    return filtered_dfs, all_indices


def filter_measurement_MvsH(unfiltered_MvsH_indices):
    #Filters the points based on the difference between the indices,
    #if it's 1 then the points are considered of interest.
    global table
    filtered_table = []
    j = 0
    for table in unfiltered_MvsH_indices:
        filtered = [table[i] for i in range(len(table)-1) if table[i+1] - table[i] == 1] # see tingimus siin ei tundu väga kindel

        if table[-1] - table[-2] == 1: #Extra condition for the last element
            filtered.append(table[-1])
        filtered_table.append(filtered)
        j += 1
        
    return filtered_table


#---------------------Separating MvsH measurements----------------

#MÕTLE ÄKKI ON PAREM MÕLEMAD ÜHEKS KOKKU PANNA IKKA
def separate_MvsH_index(MvsH_indices):
    transition_indices = []
    for indices_range in MvsH_indices:
        series = MEASUREMENT_TABLE["Magnetic Field (Oe)"].loc[indices_range]
        transition_indices.append(series.idxmin())
    return transition_indices

def separate_MvsH(separation_index, MvsH_indices):
    #Separates the points based on the index and returns the separated series in pairs
    separated_pair = []
    i = 0 #index for MvsH_indices elements
    
    for index in separation_index:
        separated = []
        min_index = min(MvsH_indices[i])
        max_index = max(MvsH_indices[i])
        sliced1 = MEASUREMENT_TABLE[["Magnetic Field (Oe)","Moment (emu)"]].loc[min_index:index-1]
        separated.append(sliced1)
        sliced2 = MEASUREMENT_TABLE[["Magnetic Field (Oe)","Moment (emu)"]].loc[index:max_index]
        separated.append(sliced2)
        i += 1
        separated_pair.append(separated)
        
    return separated_pair

#-----------------------MvsH correction-----------------------------------------------------------
def round_H_for_correction(MvsH_indices):
    #rounds the magnetic field to a nice number (100041.221 -> 100000)
    values = []
    for indices in MvsH_indices:
        MvsH = MEASUREMENT_TABLE["Magnetic Field (Oe)"].loc[indices]
        max_range = max(MvsH)
        if max_range >= 10000:
            nr = round(max_range / 1000) * 1000
        elif max_range >= 1000:
            nr = round(max_range / 100) * 100
        elif max_range >= 100:
            nr = round(max_range / 10) * 10
        else:
            nr = round(max_range)
        values.append(nr)
        
    return values


#-------Applying the corresponding correction table----------

def search_file(folder_path, number):
    #searches for the right table based on the value in the title and returns the filepath
    
    for filename in os.listdir(folder_path):
        # Split the filename on a specific character, e.g., underscore
        parts = filename.split('_') # Modify this based on your file naming convention
        # Check if the number is an exact match in any part of the filename
        
        if str(number) in parts:
            file_path = os.path.join(folder_path, filename)
            return file_path
        
    return None

folder_path = os.path.join(USER_PATH,'PdCorrection tables')
#folder_path = 'C:/Users/Kevin/Desktop/Andmetöötlus/Projekt_andmed1/PdCorrection tables' #PC
#folder_path = "C:/Users/kevin/OneDrive/Desktop/Andmetöötlus/Projekt_andmed1/PdCorrection tables" #Laptop

def nr_to_dict(numbers_to_search):
    #returns the corresponding amount of error tables in a dictionary form with the key being the value the table is for
    error_tables = {}
    
    for name in numbers_to_search:
        error_table_path = search_file(folder_path, name)
        
        if error_table_path is None:
            print(f"{name} Oe jaoks ei leia parandustabelit ")
            continue
        error_table_measurements = pd.read_csv(error_table_path, index_col=0, encoding="ANSI")
        error_tables[name] = error_table_measurements
        
    return error_tables


#---------Interpolation with the correction table--------------

def interpolate_true_field(magnetic_field_values, error_table_measurements):
    x = error_table_measurements["Magnetic Field (Oe)"]
    y = error_table_measurements["True Field (Oe)"]

    # Create an interpolation function
    interp_func = interp1d(x, y, kind='linear', fill_value='extrapolate')

    # Call the interpolation function with magnetic_field_values to get the interpolated values
    true_field_interpolated = interp_func(magnetic_field_values)
    return true_field_interpolated

def interpolate_MvsH(separated_MvsH, error_tables):
    #replaces the old field values with the correct interpolated value
    
    interpolated = []
    interpolated_dict = {}
    for i in range(len(separated_MvsH)):
        max_val = max(separated_MvsH[i][0]["Magnetic Field (Oe)"])
        
        for key in error_tables:
            if max_val - 100 <= key <= max_val + 100:
                print(f"{key} kukub {max_val} vahemikku")
                
                for j in range(len(separated_MvsH[i])):
                    magnetic_field_values = separated_MvsH[i][j]["Magnetic Field (Oe)"]
                    true_field_interpolated = interpolate_true_field(magnetic_field_values, error_tables[key])
                    
                    true_field_pd = pd.DataFrame({"True Field (Oe)": true_field_interpolated})
                    interpolated.append(true_field_pd)
                    interpolated_dict[key] = interpolated
            else:
                print(f"{key} ei kuku {max_val} vahemikku")
                
    return interpolated_dict

#------------------------------------------------------------------------------------------------

#----------------------------------------MvsT functions----------------------------------------------------
#-------------This part finds the constant magnetic field from the data (MvsT)------------------


def get_const_M_for_both(count):
    #outputs the float values for the constant magnetic fields for MvsT
    mean = stat.mean(count)
    print("\nMean H:", mean)
    const_H = []
    for index, value in count.items():
        if value > 2*mean: # muutsin nivood, kui hiljeb peaks probleeme tekkima 28.08
            print("Üle mean:", value)
            const_H.append(index)
            print(f"{index} tulbas viitab const H")
            
    return const_H

def get_const_M_for_MvsT(count):
    #outputs the float values for the constant magnetic fields for MvsT
    mean = stat.mean(count)
    print("\nMean H:", mean/2)
    const_H = []
    for index, value in count.items():
        if value > mean/2:
            print("Üle mean:", value)
            const_H.append(index)
            print(f"{index} tulbas viitab const H")
            
    return const_H

#--------------------Filtering const H points from the data for MvsT---------

def get_measurement_MvsT(const_H_values): #!!!
    #Saves all the indices of the points that are equal to the predetermined H value
    row_indices = []
    table = MEASUREMENT_TABLE['Magnetic Field (Oe)']
    
    for value in const_H_values:

        # filtered_df = table[(table == value)]
        # indices = filtered_df.index.tolist() vana versioon mis mõlema puhul korraga töötas
        # row_indices.append(indices)
        
        indices = table.index[table == value].tolist()
        #print(f"Appending indices for {value}: {indices}")
        row_indices.append(indices)

    return row_indices


#------------------Separating MvsT measurements--------------------


#SEE ON SEE ÜLDINE FUNKTSIOON MIS OSKAB ÜKSKÕIK MIS RIDA TÜKELDADA NII, KUI EKSTREEMUMID TULEVAD ÕIGESSE KOHTA JA SENI PAISTAB ET TÖÖTAB HÄSTI

#Ma tean et selle funktsiooni comment on kindlasti imeline ja küll lõpuks võib kõigi jaoks midagi sarnast teha, see muidu GPT comment

def separation_index_for_clean_series(data, column_name, n = 50): # https://stackoverflow.com/questions/48023982/pandas-finding-local-max-and-min
    """
    Find local peaks indices (maxima and minima) in a DataFrame or Series.

    Parameters:
    - data: DataFrame or Series.
    - column_name: Name of the column to analyze.
    - n: Number of points to be checked before and after.

    Returns:
    - DataFrame with 'min' and 'max' columns indicating local minima and maxima.
    """
    global max_indices, min_indices, local_peaks, katse
    katse = data
    if isinstance(data, pd.Series):
        # Convert a Series to a DataFrame with a specified column name
        data = pd.DataFrame({column_name: data})
    
    index = data.index
    
    # Find local peaks
    relative_min_indices = argrelextrema(data[column_name].values, np.less_equal, order=n)[0]
    relative_max_indices = argrelextrema(data[column_name].values, np.greater_equal, order=n)[0]
    # min_indices = argrelextrema(data[column_name].values, np.less_equal, order=n)[0]
    # max_indices = argrelextrema(data[column_name].values, np.greater_equal, order=n)[0] vana
    
    min_indices = index[relative_min_indices]
    max_indices = index[relative_max_indices]
    
    # Create a DataFrame to store results
    local_peaks = pd.DataFrame(index=index)
    local_peaks['min'] = np.nan
    local_peaks['max'] = np.nan
    
    print("min:", min_indices)
    print("max:", max_indices)
    
    # Fill in the 'min' and 'max' columns with peak values 
    local_peaks.loc[min_indices, 'min'] = data.loc[min_indices, column_name]
    local_peaks.loc[max_indices, 'max'] = data.loc[max_indices, column_name]
    # local_peaks.loc[min_indices, 'min'] = data.iloc[min_indices][column_name]
    # local_peaks.loc[max_indices, 'max'] = data.iloc[max_indices][column_name] vana versioon
    
    # Extract max min indices
    mask = local_peaks.notna()
    max_indices = local_peaks["max"].index[mask["max"]]
    min_indices = local_peaks["min"].index[mask["min"]]
    
    # Plot results
    plt.scatter(data.index, local_peaks['min'], c='r', label='Minima')
    plt.scatter(data.index, local_peaks['max'], c='g', label='Maxima')
    plt.plot(data.index, data[column_name], label=column_name)
    plt.legend()
    plt.show()
        
    return min_indices, max_indices

def separation_index_for_mixed_series(indices, column_name):
    tester = []
    for index in indices:
        
        #print("indeksid:", index)
        
        filtered = MEASUREMENT_TABLE[column_name].loc[index]
        #print("filtered:", filtered)
        
        separation_index = separation_index_for_clean_series(filtered, column_name)
        
        tester.append(separation_index)
        print("check")
        
    return tester
    

def separate_MvsT(separation_index, MvsT_indices):
    #Separates the points based on the index and returns the separated series in pairs
    separated_pair = []
    
    min_index_list = separation_index[0].tolist() #Siit jätka, probleem et enne üli üks tuple nüüd list kus min/max indeksid
    max_index_list = separation_index[1].tolist()

    i = 0
    j = 0
    
    for indices in MvsT_indices:

        tabel = MEASUREMENT_TABLE[["Temperature (K)","Moment (emu)"]].loc[indices]
    
        for max_index in max_index_list:
            
            separated = []

            if max_index in indices:
                
                sliced1 = tabel.loc[min_index_list[j]:max_index]
                separated.append(sliced1)
                
                if j == len(min_index_list) - 1:
                    sliced2 = tabel.loc[max_index+1:min_index_list[j]]
                    
                else:
                    sliced2 = tabel.loc[max_index+1:min_index_list[j+1]]
            else:
                max_index_list = max_index_list[j:]
                
                break
            
            separated.append(sliced2)
            separated_pair.append(separated)
            
            j += 1
            
        i += 1
        
    return separated_pair


def separate_MvsT_index_for_both(MvsT_indices): #!!! saab vist kokku panna indekseerimise
    #Finds the index where to separate the measurement into multiple series,
    #each series being the temperature increase from lower to higher.
    #if the previous value is bigger than the current one, then it's considered the breaking point
    transition_indices = []
    
    for indices_range in MvsT_indices:
        series = MEASUREMENT_TABLE["Temperature (K)"].loc[indices_range]
        previous_index = indices_range[0]  # Initialize with the starting index
        #print(series)
        for index in indices_range:
            if series.loc[index] < series.loc[previous_index]: #loc for index value based indexing, not iloc for label
                transition_indices.append(index)
            previous_index = index
            
    # if len(transition_indices) == 0: # the case where there was only one up/down cycle measured at that H value
    #     return MvsT_indices
    
    return transition_indices

def separate_MvsT_for_both(separation_index, MvsT_indices): 
    #Separates the points based on the index and returns the separated series in pairs
    separated_pair = []
    
    i = 0 #index for MvsT_indices elements
    
    if len(separation_index) == 0: # the case where there was only only up cycle measured at that H
        for indices in MvsT_indices:
            whole = [MEASUREMENT_TABLE[["Temperature (K)","Moment (emu)"]].loc[indices]]
            separated_pair.append(whole)
        return separated_pair
    
    for index in separation_index:
        print("Mitmes separate_MvsT tsükkel: ",i)
        separated = []
        min_index = min(MvsT_indices[i])
        max_index = max(MvsT_indices[i])
        sliced1 = MEASUREMENT_TABLE[["Temperature (K)","Moment (emu)"]].loc[min_index:index-1]
        separated.append(sliced1)
        sliced2 = MEASUREMENT_TABLE[["Temperature (K)","Moment (emu)"]].loc[index:max_index]
        separated.append(sliced2)
        i += 1
        separated_pair.append(separated)
    
    return separated_pair
#-------------------------------------------------------------------------------------------------------------------------------------
def separated_into_dict_pair(separated_pairs, const_val, column):
    # for interval:value pairs
    raamat = {}
    
    for const in const_val:
        
        for val in separated_pairs:
            
            index_to_check = val[0].index[0]
            val_to_check = MEASUREMENT_TABLE[column].iloc[index_to_check]
            
            # print("val",val_to_check)
            # print("int val", round(val_to_check,0))
            # print("const", const)
            # print("const int", int(const))
            
            if val_to_check == const or round(val_to_check, 1) == round(const,1):
                key = const
                # print("millal:", key, "\n")
                if key not in raamat:
                    raamat[key] = []  # Create an empty list for this key if it doesn't exist
                    
                raamat[key].append(val)  # Append the value to the list associated with the key

    return raamat


#-----------------------PLOTTING THE DATA---------------------------------------------------------
COLORS = list(mcolors.BASE_COLORS.keys())
MEASUREMENT_TABLE["color"] = "black"


def plot_MvsT(raamat, const_H_values, MvsT_indices):
    #Plots the MvsT measurement pair with different colors
    
    for index in MvsT_indices:
        
        MEASUREMENT_TABLE.loc[index, "color"] = COLORS[0]
        COLORS.pop(0)
    
    for key in raamat:
        
        for df in raamat[key]:

            fig, ax = plt.subplots()
            T1 = df[0]["Temperature (K)"]
            M1 = df[0]["Moment (emu)"]
            T2 = df[1]["Temperature (K)"] if len(df) > 1 else None
            M2 = df[1]["Moment (emu)"] if len(df) > 1 else None
            
            ax.plot(T1,M1,color = "green", label = "Ascending")
            ax.plot(T2,M2,color = "red", label = "Descending") if len(df) > 1 else None #, marker = "o") #descending ei pea paika kui on alt üle > alt üles mõõtmine
            ax.set_title(f"M vs T at {key} Oe")
            ax.set_xlabel("Temperature (K)")
            ax.set_ylabel("Moment (emu)")
            ax.legend() #Hetkel legend nimetab selle järgi et esimene tsükkel on kasvav ja teine kahanev ehk eeldus et mõõtmisel temp algas kasvamisest
            ax.grid(True)
            fig.savefig(os.path.join(save_to_path,f'MvsT_graph_at_{key}K.png'),bbox_inches = "tight", dpi = 200)
        
    return None

def plot_MvsH(raamat, const_T_values, interpolated_MvsH, MvsH_indices):
    #Plots the MvsH measurement pair with different colors

    j = 0
    
    for index in MvsH_indices:
        
        MEASUREMENT_TABLE.loc[index, "color"] = COLORS[0]
        COLORS.pop(0)
    
    for key1 in raamat:
        
        for df in raamat[key1]:
            
            fig, ax = plt.subplots()
            H1 = df[0]["Magnetic Field (Oe)"]
            M1 = df[0]["Moment (emu)"]
            H2 = df[1]["Magnetic Field (Oe)"]
            M2 = df[1]["Moment (emu)"] 
            
            for key2 in interpolated_MvsH:

                # See on interpoleerimise jaoks
                H3 = interpolated_MvsH[key2][j]["True Field (Oe)"]
                H4 = interpolated_MvsH[key2][j+1]["True Field (Oe)"]
                
                if len(H3) != len(M1) or len(H4) != len(M2):
                    #Siin ongi see, et ei ole interpoleeritud punkte, siis x/y erineva pikkusega ja ei saa joonistada
                    print(f"Warning: Length mismatch between True Field and Moment data for {key1} K.") #error selles,et pole correction tabelit iga välja tugevuse jaoks
                    continue  # Continue to the next iteration
                    
                j += 2
                ax.plot(H3, M1, color = "red", label = "True Field")
                ax.plot(H4, M2, color = "red")
                ax.legend()
                
            ax.plot(H1,M1,color = "blue", label = "Ascending")
            ax.plot(H2,M2,color = "orange", label = "Descending")
            
            val = int(key1) #!!! siin peaks ka probleem tekkima, kui mitu möötmist samal tempil
            
            plot_title = f"M vs H at {val} K"
            ax.set_title(plot_title)
            ax.set_xlabel("Magnetic field (Oe)")
            ax.set_ylabel("Moment (emu)")
            ax.legend() #Hetkel legend nimetab selle järgi et esimene tsükkel on kasvav ja teine kahanev ehk eeldus et mõõtmisel temp algas kasvamisest
            ax.grid(True)
            #fig.savefig(f"C:/Users/kevin/OneDrive/Desktop/Andmetöötlus/Projekt_andmed1/MvsH_graph_at_{val}K.png",bbox_inches = "tight", dpi = 200) #laptop
            fig.savefig(os.path.join(save_to_path,f'MvsH_graph_at_{val}K.png'),bbox_inches = "tight", dpi = 200) #PC
        
    return None

def plot_timeseries(): 

    # Create subplots with shared x-axis
    fig, axes = plt.subplots(nrows=3, sharex=True)
    
    # Plot data on each subplot
    time_axis = "Time Stamp (sec)"
    MEASUREMENT_TABLE.plot(x=time_axis, y="Temperature (K)", kind="scatter", c=MEASUREMENT_TABLE["color"].values, ax=axes[0])
    MEASUREMENT_TABLE.plot(x=time_axis, y="Moment (emu)", kind="scatter", c=MEASUREMENT_TABLE["color"].values, ax=axes[1])
    MEASUREMENT_TABLE.plot(x=time_axis, y="Magnetic Field (Oe)", kind="scatter", c=MEASUREMENT_TABLE["color"].values, ax=axes[2])
    
    # Customize axes labels and other properties
    
    axes[0].set_ylabel("Temperature (K)")
    axes[1].set_ylabel("Moment (emu)")
    axes[2].set_ylabel("Magnetic Field (Oe)")
    axes[-1].set_xlabel("Time Stamp (sec)")
    
    # Show the plot
    plt.show()
    
    return None



#-------------Determining the measurement type from the data (MvsH/MvsT or both)-----------------


def measurement_type_check(RATIO_LEVEL = 0.07):
    #Checks what measurements the file contains
    
    global codes_T, uniques_T, codes_H, uniques_H
    
    type_token = {"Temperature": None, "Field": None}
    
    dataset_T = MEASUREMENT_TABLE["Temperature (K)"]
    dataset_H = MEASUREMENT_TABLE["Magnetic Field (Oe)"]
    
    smoothed_dataset_T = dataset_T.round(decimals=1)
    smoothed_dataset_H = dataset_H.round(decimals=1)
    
    
    codes_T, uniques_T = pd.factorize(smoothed_dataset_T)
    codes_H, uniques_H = pd.factorize(smoothed_dataset_H)
    
    ratio_T = uniques_T.size/smoothed_dataset_T.size
    ratio_H = uniques_H.size/smoothed_dataset_H.size
    
    if ratio_T < RATIO_LEVEL: #discrete
    
        if ratio_H < RATIO_LEVEL: #discrete
        
            type_token["Temperature"] = "discrete"
            type_token["Field"] = "discrete"
            print("T discrete, H discrete = error \n")
            
            return type_token
        
        else: #continous
        
            type_token["Temperature"] = "discrete"
            type_token["Field"] = "continous"
            print("T discrete, H continous = MvsH \n")
            
            return type_token
            
    else: #continous
    
        if ratio_H < RATIO_LEVEL: #discrete
        
            type_token["Temperature"] = "continous"
            type_token["Field"] = "discrete"
            print("T continous, H discrete = MvsT \n")
            
            return type_token
        
        else: #continous
        
            type_token["Temperature"] = "continous"
            type_token["Field"] = "continous"
            print("T continous, H continous = both \n")
            
            return type_token

measurement_type_token = measurement_type_check()

def MvsH():
    
    global ranges_temp, intervals_temp, const_T_interval, const_T_values, unfiltered_MvsH_T_values, unfiltered_MvsH_indices,\
        MvsH_indices, separated_MvsH_indices, separated_MvsH
        
    ranges_temp = min_max_range(TEMPERATURE)
    intervals_temp = TEMPERATURE.groupby(pd.cut(TEMPERATURE, ranges_temp)).count()#Groups the temperatures based on the range and returns the count of each range
    #print(intervals)
    const_T_interval = mood_temp_for_MvsH(intervals_temp)
    const_T_values = get_const_temp(const_T_interval)
    
    #SIIN VAATA filtered_values üle, äkki saad hiljem kasutada
    unfiltered_MvsH_T_values, unfiltered_MvsH_indices = get_measurement_MvsH(const_T_values)
    MvsH_indices = filter_measurement_MvsH(unfiltered_MvsH_indices)
    
    #!!! vb siia ka separation_index_for_series panna
    separated_MvsH_indices = separate_MvsH_index(MvsH_indices) #the indices where the separation is going to be done
    separated_MvsH = separate_MvsH(separated_MvsH_indices, MvsH_indices)
    
    return separated_MvsH, MvsH_indices, const_T_values

def MvsT():
    global H_count, const_H_values, MvsT_indices, separated_MvsT_indices, separated_MvsT
    
    H_count = MAGNETIC_FIELD.value_counts()
    const_H_values = get_const_M_for_MvsT(H_count)
    
    MvsT_indices = get_measurement_MvsT(const_H_values)
    separated_MvsT_indices = separation_index_for_clean_series(MEASUREMENT_TABLE, "Temperature (K)")# the indices where the separation is going to be done
    separated_MvsT = separate_MvsT(separated_MvsT_indices, MvsT_indices)
    
    return separated_MvsT, MvsT_indices, const_H_values

def what_path_main(measurement_type_token):
    
    global ranges_temp, intervals_temp, const_T_interval, const_T_values, H_count, const_H_values, unfiltered_MvsH_T_values, unfiltered_MvsH_indices,\
        MvsH_indices, MvsT_indices, separated_MvsT_indices, separated_MvsT, separated_MvsH_indices, separated_MvsH, min_max_MvsH_val, error_tables,\
        interpolated_MvsH, dict_MvsT, dict_MvsH
        
    if measurement_type_token["Temperature"] == "discrete" and measurement_type_token["Field"] == "continous":
        
        separated_MvsH, MvsH_indices, const_T_values = MvsH()
        
        dict_MvsH = separated_into_dict_pair(separated_MvsH, const_T_values, "Temperature (K)")
        
        min_max_MvsH_val = round_H_for_correction(MvsH_indices) #Magnetic field value indicating the +- field value
        error_tables = nr_to_dict(min_max_MvsH_val)
        interpolated_MvsH = interpolate_MvsH(separated_MvsH, error_tables)
        
        plot_MvsH(dict_MvsH, const_T_values, interpolated_MvsH, MvsH_indices)
        
        print("\n Kas jõudis MvsH lõppu?")
        
        return separated_MvsH, MvsH_indices, const_T_values #HETKEL OUTPUT SELLEKS, ET SAAKS MUUTUJAID JÄLGIDA
    
    elif measurement_type_token["Temperature"] == "continous" and measurement_type_token["Field"] == "discrete":
        
        separated_MvsT, MvsT_indices, const_H_values = MvsT()
        
        dict_MvsT = separated_into_dict_pair(separated_MvsT, const_H_values, "Magnetic Field (Oe)")
        
        plot_MvsT(dict_MvsT, const_H_values, MvsT_indices)
        
        print("\n Kas jõudis MvsT lõppu?")
        
        return separated_MvsT, MvsT_indices, const_H_values #HETKEL OUTPUT SELLEKS, ET SAAKS MUUTUJAID JÄLGIDA
                
    
    elif measurement_type_token["Temperature"] == "continous" and measurement_type_token["Field"] == "continous":

        #This part finds the constant temperatures from the data (MvsH)
        ranges_temp = min_max_range(TEMPERATURE)
        intervals_temp = TEMPERATURE.groupby(pd.cut(TEMPERATURE, ranges_temp)).count()#Groups the temperatures based on the range and returns the count of each range
        const_T_interval = mood_temp_for_both(intervals_temp)
        const_T_values = get_const_temp(const_T_interval)
        
        #This part finds the constant magnetic field from the data (MvsT)
        H_count = MAGNETIC_FIELD.value_counts()
        const_H_values = get_const_M_for_both(H_count)
        
        #Filtering const temp points from the data for MvsH
        unfiltered_MvsH_T_values, unfiltered_MvsH_indices = get_measurement_MvsH(const_T_values)
        MvsH_indices = filter_measurement_MvsH(unfiltered_MvsH_indices)
        
        #Filtering const H points from the data for MvsT
        MvsT_indices = get_measurement_MvsT(const_H_values)
        
        #Separating MvsT measurements
        separated_MvsT_indices = separate_MvsT_index_for_both(MvsT_indices)# the indices where the separation is going to be done
        separated_MvsT = separate_MvsT_for_both(separated_MvsT_indices, MvsT_indices)
        
        dict_MvsT = separated_into_dict_pair(separated_MvsT, const_H_values, "Magnetic Field (Oe)")
        
        #Separating MvsH measurements
        separated_MvsH_indices = separate_MvsH_index(MvsH_indices) #the indices where the separation is going to be done
        separated_MvsH = separate_MvsH(separated_MvsH_indices, MvsH_indices)
        dict_MvsH = separated_into_dict_pair(separated_MvsH, const_T_values, "Temperature (K)")
        
        #MvsH correction
        min_max_MvsH_val = round_H_for_correction(MvsH_indices) #Magnetic field value indicating the +- field value
        error_tables = nr_to_dict(min_max_MvsH_val)
        interpolated_MvsH = interpolate_MvsH(separated_MvsH, error_tables)
        
        #Plots
        plot_MvsH(dict_MvsH, const_T_values, interpolated_MvsH, MvsH_indices)
        plot_MvsT(dict_MvsT, const_H_values, MvsT_indices)
        
        print("\n Kas jõudis duo lõppu?")

        return separated_MvsH, MvsH_indices, const_T_values, separated_MvsT, MvsT_indices, const_H_values
    
what_path_main(measurement_type_token)
plot_timeseries()

tester = separation_index_for_mixed_series(MvsT_indices, "Temperature (K)")



def plot_to_csv(indices, columns, folder_path, dType):
    count = 0
    for index in indices:

        # Select the desired rows and columns
        selected_data = MEASUREMENT_TABLE.loc[index][columns]
        
        # Define a unique file name based on the index
        file_name = f'{dType}_data_{count}.csv'
        
        # Create the full path to the CSV file
        file_path = os.path.join(folder_path, file_name)
        
        # Save the selected data to the CSV file
        selected_data.to_csv(file_path, index=False)
        count += 1
    
    return None

columns_to_save_MvsH = ["Magnetic Field (Oe)", "Moment (emu)"]
columns_to_save_MvsT = ["Temperature (K)", "Moment (emu)"]

#folder_path = 'C:/Users/kevin/OneDrive/Desktop/Andmetöötlus'  # laptop
folder_path = "C:/Users/Kevin/Desktop/Andmetöötlus/Projekt_andmed1" #PC

if "MvsH_indices" not in globals():
    plot_to_csv(MvsT_indices, columns_to_save_MvsT, folder_path, "MvsT")
elif "MvsT_indices" not in globals():
    plot_to_csv(MvsH_indices, columns_to_save_MvsH, folder_path, "MvsH")
else:
    plot_to_csv(MvsT_indices, columns_to_save_MvsT, folder_path, "MvsT")
    plot_to_csv(MvsH_indices, columns_to_save_MvsH, folder_path, "MvsH")



# class PlotSetup:
#     COLORS = list(mcolors.BASE_COLORS.keys())
#     def __init__(self, Dict, const):
#         self.Dict = Dict
#         self.const = const
#         self.interpolated = None
#         self.indices = None
        
#     def setColorIndex(self, indices):
#         self.indices = indices
        
#         for index in indices:
            
#             MEASUREMENT_TABLE.loc[index, "color"] = self.COLORS[0]
#             self.COLORS.pop(0)
#tester = PlotSetup(dict_MvsT, const_T_values)
#Parameetrite test

#Dataframe.to_excel(os.path.join(save_to_path,"parameterTest.xlsx"), index = False) sellega saab exceli faili näiteks teha parameetrite jaoks

#-----------------------DIVIDING MEASUREMENTS BY SAMPLE SIZE--------------------------------------

# def calc_volume(size, thickness):
#     if size is None or thickness is None:
#         return None
#     volume = size*thickness
#     return volume

# volume = calc_volume(size, thickness)

# def Moment_divided_mass(moment, mass):
#     if mass is None:
#         return None
#     M_divided_mass = moment/mass
#     return M_divided_mass

# def Moment_divided_size(moment, size):
#     if size is None:
#         return None
#     M_divided_size = moment/size
#     return M_divided_size

# def Moment_divided_volume(moment, volume):
#     if volume is None:
#         return None
#     M_divided_volume = moment/volume
#     return M_divided_volume

# M_div_mass = Moment_divided_mass(moment,mass)
# M_div_size = Moment_divided_size(moment, size)
# M_div_volume = Moment_divided_volume(moment, volume)