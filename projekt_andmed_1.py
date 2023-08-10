# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 13:00:22 2023

@author: kevin
"""
import os
import statistics as stat
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog
from tkinter import simpledialog
from scipy.interpolate import interp1d

#-------------- OPENING THE FILE AND INDEXING IT -------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
root = tk.Tk()
root.wm_attributes('-topmost', 1)
root.withdraw()

file_path = filedialog.askopenfilename()

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
    print(f"Checking: {name}, {info}")#Hetkel jääb
    
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

header_mass = header_value_check(SAMPLE_MASS, SAMPLE_MASS_NAME)

def get_mass(mass):
    #If the mass > 1, indicating that it's input is in mg, divides it by a 1000 to get g
    if mass is None:
        return None
    val = mass[0] #Mass value
    if val > 1:
        val= val/1000
        print(f'Sample mass is {val:.5f} grams \n')
        
    return val

mass = get_mass(header_mass)

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
data_measurements = get_measurements(file_path)[get_measurements(file_path)["Transport Action"] != 2]

#HETKEL ACMS FAILIGA EI TÖÖTA
measurement_table = data_measurements[["Time Stamp (sec)", "Temperature (K)", "Magnetic Field (Oe)",
                                       "Moment (emu)", "M. Std. Err. (emu)"]]

def get_Temp(data):
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

temperature = get_Temp(data_measurements)
magnetic_field = get_H(data_measurements)
moment = get_M(data_measurements, token)
#--------------------------------------------- MEASUREMENT ANALYSIS----------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------


#-------------Determening the measurement type from the data (MvsH/MvsT)-----------------

#-----------------This part finds the constant temperatures from the data (MvsH)-------------------

#ON VAJA VEEL TEISI KONTROLLE SELLEKS KUI ANDMETES ON AINULT ÜHTE TÜÜPI MÕÕTMISED

def min_max_range(data):
    #min and max value for temperature range with step 0.2
    #divides the temperature range into intervals for temperature distribution in the data
    start = min(data)
    stop = max(data)+0.5
    step = 0.2
    values = []
    current_value = start
    
    while current_value < stop:
        values.append(current_value)
        current_value += step
        
    return values

ranges_temp = min_max_range(temperature)
intervals = temperature.groupby(pd.cut(temperature, ranges_temp)).count()#Groups the temperatures based on the range and returns the count of each range
index_val = intervals.index.values#values of the indices

def mood_temp(intervals):
    #determines at which temperatures the MvsH measurement was made
    #outputs the interval which satisfies 20*mean
    mean = stat.mean(intervals)
    const_temp = []
    print("Mean temp:", mean)
    for index, value in intervals.items():
        if value > 20*mean:
            const_temp.append(index)
            
    return const_temp

const_T_interval = mood_temp(intervals)

def get_const_temp(const):
    #outputs the integer values for the constant temperatures
    const_temp = []
    x = 0
    for val in const:
        const_temp.append(int(val.right))
        print(const_temp[x],"K mõõdetud väljasõltuvust")
        x +=1
        
    return const_temp

const_T_values = get_const_temp(const_T_interval)

#-------------This part finds the constant magnetic field from the data (MvsT)------------------

M_count = magnetic_field.value_counts()


def get_const_M(count):
    #outputs the float values for the constant magnetic fields for MvsT
    mean = stat.mean(count)
    print("Mean H:", mean)
    const_H = []
    for index, value in count.items():
        if value > mean:
            const_H.append(index)
            print(f"{index} tulbas viitab const H")
            
    return const_H

const_H_values = get_const_M(M_count)

#------------Filtering const temp points from the data for MvsH------


def get_measurement_MvsH(const_T_values):
    #Saves all the indices of the points that are close enough to the constant temperature
    table = measurement_table['Temperature (K)']
    filtered_dfs = []
    all_indices = []
    
    for value in const_T_values:
        lower = value - 0.05
        upper = value + 0.05
        filtered_df = table[(table >= lower) & (table <= upper)]
        indices = filtered_df.index.tolist()
        filtered_dfs.append(filtered_df)
        all_indices.append(indices)
        
    return filtered_dfs, all_indices

#SIIN VAATA filtered_values üle, äkki saad hiljem kasutada
unfiltered_MvsH_T_values, unfiltered_MvsH_indices = get_measurement_MvsH(const_T_values)

def filter_measurement_MvsH(unfiltered_MvsH_indices):
    #Filters the points based on the difference between the indices,
    #if it's 1 then the points are considered of interest.
    filtered_table = []
    
    for table in unfiltered_MvsH_indices:
        filtered = [table[i] for i in range(len(table)-1) if table[i+1] - table[i] == 1]
        if table[-1] - table[-2] == 1: #Extra condition for the last element
            filtered.append(table[-1])
        filtered_table.append(filtered)
        
    return filtered_table

MvsH_indices = filter_measurement_MvsH(unfiltered_MvsH_indices)

#--------------------Filtering const H points from the data for MvsT---------

def get_measurement_MvsT(const_H_values):
    #Saves all the indices of the points that are equal to the predetermined H value
    row_indices = []
    table = measurement_table['Magnetic Field (Oe)']
    
    for value in const_H_values:
        filtered_df = table[(table == value)]
        indices = filtered_df.index.tolist()
        row_indices.append(indices)
        
    return row_indices

MvsT_indices = get_measurement_MvsT(const_H_values)

#------------------Separating MvsT measurements--------------------

#MÕELDA KA SELLELE, ET HETKEL TUNNEB ÄRA AINULT MÕÕTMISI KUS TEMP LIIGUB ALT -> ÜLES

#HETKEL OSKAB VÕRRELDA AINULT MÕÕTMIST KUS KÕIK PUNKTID KASVAVAD MONOTOONSELT
def separate_MvsT_index(measurement_table, MvsT_indices):
    #Finds the index where to separate the measurement into multiple series,
    #each series being the temperature increase from lower to higher.
    #if the previous value is bigger than the current one, then it's considered the breaking point
    transition_indices = []
    
    for indices_range in MvsT_indices:
        series = measurement_table["Temperature (K)"].loc[indices_range]
        previous_index = indices_range[0]  # Initialize with the starting index
        
        for index in indices_range:
            if series.loc[index] < series.loc[previous_index]: #loc for index value based indexing, not iloc for label
                transition_indices.append(index)
            previous_index = index
            
    return transition_indices

separated_MvsT_indices = separate_MvsT_index(measurement_table, MvsT_indices) #the indices where the separation is going to be done

def separate_MvsT(separation_index, MvsT_indices):
    #Separates the points based on the index and returns the separated series in pairs
    separated_pair = []
    i = 0 #index for MvsT_indices elements
    
    for index in separation_index:
        separated = []
        min_index = min(MvsT_indices[i])
        max_index = max(MvsT_indices[i])
        sliced1 = measurement_table[["Temperature (K)","Moment (emu)"]].loc[min_index:index-1]
        separated.append(sliced1)
        sliced2 = measurement_table[["Temperature (K)","Moment (emu)"]].loc[index:max_index]
        separated.append(sliced2)
        i += 1
        separated_pair.append(separated)
        
    return separated_pair

separated_MvsT = separate_MvsT(separated_MvsT_indices, MvsT_indices)

#---------------------Separating MvsH measurements----------------

#MÕTLE ÄKKI ON PAREM MÕLEMAD ÜHEKS KOKKU PANNA IKKA
def separate_MvsH_index(measurement_table, MvsH_indices):
    transition_indices = []
    for indices_range in MvsH_indices:
        series = measurement_table["Magnetic Field (Oe)"].loc[indices_range]
        transition_indices.append(series.idxmin())
    return transition_indices

separated_MvsH_indices = separate_MvsH_index(measurement_table, MvsH_indices) #the indices where the separation is going to be done


def separate_MvsH(separation_index, MvsH_indices):
    #Separates the points based on the index and returns the separated series in pairs
    separated_pair = []
    i = 0 #index for MvsH_indices elements
    
    for index in separation_index:
        separated = []
        min_index = min(MvsH_indices[i])
        max_index = max(MvsH_indices[i])
        sliced1 = measurement_table[["Magnetic Field (Oe)","Moment (emu)"]].loc[min_index:index-1]
        separated.append(sliced1)
        sliced2 = measurement_table[["Magnetic Field (Oe)","Moment (emu)"]].loc[index:max_index]
        separated.append(sliced2)
        i += 1
        separated_pair.append(separated)
        
    return separated_pair

separated_MvsH = separate_MvsH(separated_MvsH_indices, MvsH_indices)



#-----------------------DIVIDING MEASUREMENTS BY SAMPLE SIZE--------------------------------------

#thickness = float(input("Insert sample thickness in number format: "))
#thickness = simpledialog.askfloat("Enter a number", "Please enter a number:")
#thickness = 1

def calc_volume(size, thickness):
    if size is None or thickness is None:
        return None
    volume = size*thickness
    return volume

volume = calc_volume(size, thickness)

def Moment_divided_mass(moment, mass):
    if mass is None:
        return None
    M_divided_mass = moment/mass
    return M_divided_mass

def Moment_divided_size(moment, size):
    if size is None:
        return None
    M_divided_size = moment/size
    return M_divided_size

def Moment_divided_volume(moment, volume):
    if volume is None:
        return None
    M_divided_volume = moment/volume
    return M_divided_volume

M_div_mass = Moment_divided_mass(moment,mass)
M_div_size = Moment_divided_size(moment, size)
M_div_volume = Moment_divided_volume(moment, volume)

#-----------------------MvsH correction-----------------------------------------------------------

#------Finding the magnetic field range------------------

def round_H(MvsH_indices):
    #rounds the magnetic field to a nice number (100041.221 -> 100000)
    values = []
    for indices in MvsH_indices:
        MvsH = measurement_table["Magnetic Field (Oe)"].loc[indices]
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

min_max_MvsH_val = round_H(MvsH_indices) #Magnetic field value indicating the +- field value

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

#folder_path = 'C:/Users/Kevin/Desktop/Andmetöötlus/Projekt_andmed1/PdCorrection tables' #PC
folder_path = "C:/Users/kevin/OneDrive/Desktop/Andmetöötlus/Projekt_andmed1/PdCorrection tables" #Laptop

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

error_tables = nr_to_dict(min_max_MvsH_val)


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

interpolated_MvsH = interpolate_MvsH(separated_MvsH, error_tables)
        

#-----------------------PLOTTING THE DATA---------------------------------------------------------

def plot_MvsT(separated_MvsT, MvsT_indices, const_H_values):
    #Plots the MvsT measurement pair with different colors
    for i in range(len(MvsT_indices)):
        fig, ax = plt.subplots()
        T1 = separated_MvsT[i][0]["Temperature (K)"]
        M1 = separated_MvsT[i][0]["Moment (emu)"]
        T2 = separated_MvsT[i][1]["Temperature (K)"]
        M2 = separated_MvsT[i][1]["Moment (emu)"]
        
        ax.plot(T1,M1,color = "red")
        ax.plot(T2,M2,color = "green")
        val = const_H_values[i]
        ax.set_title(f"M vs T at {val} Oe") #hetkel
        ax.set_xlabel("Temperature (K)")
        ax.set_ylabel("Moment (emu)")
        #ax.legend() #KAS SIIN KA LEGENDI VAJA?
        ax.grid(True)
        fig.savefig(f"C:/Users/kevin/OneDrive/Desktop/Andmetöötlus/Projekt_andmed1/MvsT_graph_at_{val}K.png",bbox_inches = "tight", dpi = 200)
        
    return None

plot_MvsT(separated_MvsT, MvsT_indices, const_H_values)

def plot_MvsH(separated_MvsH, const_T_values, interpolated_MvsH):
    #Plots the MvsH measurement pair with different colors
    i = 0
    
    for i in range(len(separated_MvsH)):
        fig, ax = plt.subplots()
        H1 = separated_MvsH[i][0]["Magnetic Field (Oe)"]
        M1 = separated_MvsH[i][0]["Moment (emu)"]
        H2 = separated_MvsH[i][1]["Magnetic Field (Oe)"]
        M2 = separated_MvsH[i][1]["Moment (emu)"] 
        
        if i == 1:
            for key in interpolated_MvsH:
                H3 = interpolated_MvsH[key][0]["True Field (Oe)"]
                H4 = interpolated_MvsH[key][1]["True Field (Oe)"]
                ax.plot(H3, M1, color = "red", label = "True Field")
                ax.plot(H4, M2, color = "red")
                ax.legend()

        ax.plot(H1,M1,color = "blue", label = "Ascending")
        ax.plot(H2,M2,color = "orange", label = "Descending")
        val = const_T_values[i]
        ax.set_title(f"M vs H at {const_T_values[i]} K")
        ax.set_xlabel("Magnetic field (Oe)")
        ax.set_ylabel("Moment (emu)")
        ax.legend() #Hetkel loodab lihtsalt sellele, et algav tsükkel on kasvav, KÜSI ÜLE!
        ax.grid(True)
        fig.savefig(f"C:/Users/kevin/OneDrive/Desktop/Andmetöötlus/Projekt_andmed1/MvsH_graph_at_{val}K.png",bbox_inches = "tight", dpi = 200)
        i += 1
        
    return None

plot_MvsH(separated_MvsH, const_T_values, interpolated_MvsH)

test1

