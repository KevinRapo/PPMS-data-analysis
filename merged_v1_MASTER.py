# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 15:25:54 2023

@author: Kevin
"""

import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import messagebox
from tkinter import filedialog
from scipy.interpolate import interp1d
from scipy.signal import argrelextrema
import copy



USER_PATH = os.getcwd()

#-------------- OPENING THE FILE AND INDEXING IT -------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
root = tk.Tk()
root.wm_attributes('-topmost', 1)
root.withdraw()

        
def askNewDatafilePath():
    """
    Ask for user input for a datafile.

    Returns
    -------
    file_path : STRING
        File path for the measurement file

    """
    
    file_path = filedialog.askopenfilename()
    
    return file_path


def readDatafile(file_path):
    '''
    Reads in the data file and separates the header and measurement info into pandas dataframes.
    
    Parameters
    ----------
    file_path : STRING
        file path for the measurement file

    Returns
    -------
    header : PANDAS DATAFRAME
        header info   
    data : PANDAS DATAFRAME
        measurement data

    '''
    #Opens the selected data file
    with open(file_path, 'r') as f:
        i = 0 #Rowindex, how many rows until header 
        
        for line in f:
            i += 1
            if '[Data]' in line:
                f.close()
                break
    
         
    header = pd.read_csv(file_path, nrows = i , index_col = 2, encoding = 'ANSI',names = ['1','2','3','4'], on_bad_lines = 'skip')
    data = pd.read_csv(file_path, skiprows = i, header = 0, encoding = "ANSI")
    data = data[data["Transport Action"] == 1]
    
    return header, data

 
def determineDatafileType(header):
    """
    Determine if its VSM or ACMS datafile.
    
    Headers of VSM and ACMS files are similiar. DATA columns of those files have a difference in the Moment column. 
    In VSM the column is named Moment (emu), while in ACMS its named DC Moment (emu)

    Parameters
    ----------
    header : PANDAS DATAFRAME
        Measurement file header

    Returns
    -------
    token : STRING
        Data file type, "VSM" or "ACMS"

    """

    token = "error - unknown datafile format"
    
    option_specific_line = header.iloc[1, 0]
    
    if "VSM" in option_specific_line:
        print("\nThis is a VSM data file \n")
        token = "VSM"
        
    elif "ACMS" in option_specific_line:
        print("\nThis is an ACMS data file \n")
        token = "ACMS"
        
    return token
    
    
def headerValueCheck(header, sample_property):
    """
    Text parsing function to return sample property value and unit from the header if they exist.
    
    Uses a helper function extractFloatWithUnit to parse the property.
    
    Parameters
    ----------
    header : PANDAS DATAFRAME
        measurement file header.
    sample_property : STRING
        PROPERTY TO CHECK.

    Returns
    -------
    float_val : TYPE
        description.
    unit : TYPE
        DESCRIPTION.
        
        
    None : Retuns None if no match found
    """
    
    header_column = header['2']

    print(f"Checking:{sample_property}, {header_column[sample_property]}")
    
    if not isinstance(header_column[sample_property], str) and np.isnan(header_column[sample_property]): #Checks whether the value is not a string and if it is a nan value
        print(f"NO VALID VALUE FOR {sample_property}, value is NAN \n")
        return None
    
    def extractFloatWithUnit(string):
        """
        Text parsing helper function to help with sample properties, uses regex to separate the input string into
        a (float, unit) format if it has units, if no units (float, None) format, if no value (None)
        
        Parameters
        ----------
        string : STRING
            the parameter string to examine.

        Returns
        -------
        float_val : FLOAT
            sample parameter value.
        unit : STRING
            sample parameter unit.
            
            
        None : Returns None if no match found
        """

        regex = r'^([\d.]+)\s*([a-zA-Z]+(\d)?)?$'
        match = re.search(regex, string)
        
        if match:
            float_str = match.group(1)
            unit = match.group(2)
            float_val = float(float_str)
            print(f"SAMPLE MASS :{float_val}, {unit}")
            return (float_val, unit)
        else:
            return None
        
    match = extractFloatWithUnit(header_column[sample_property])
    
    if match is None: #condition for extract_float_with_unit when it didn't find a match for a desired format therefore being a string
        print(f"{sample_property} had an input put it is not in the valid format: {header_column[sample_property]} \n")
        return None
    
    float_val = match[0]

    unit = match[1]

    
    return float_val, unit
  

def getMassInGrams(header):
    """
    Returns parsed MASS from Header

    Parameters
    ----------
    header : PANDAS DATAFRAME
        Measurement file header

    Returns
    -------
    mass : FLOAT
        Mass in grams.

    """
    
    #If the mass > 1, indicating that it's input is in mg, divides it by a 1000 to get g
    parameter = "SAMPLE_MASS"
    
    mass_unit = headerValueCheck(header, parameter)

    
    if mass_unit is None:
        return None
    
    mass = mass_unit[0]
    unit = mass_unit[1]
    
    if unit == 'g':
        return mass
    
    if unit == 'mg':
        return mass/1000
    
    if unit is None and mass > 1:
        mass= mass/1000
        print(f'Sample mass is {mass:.5f} grams \n')
        
    return mass


#Parsed sample size
def getAreaCM2(header):
    """
    Returns parsed area if it's in the header

    Parameters
    ----------
    header : PANDAS DATAFRAME
        Measurement file header.

    Returns
    -------
    area : FLOAT
        Area value
    
    None : Returns None if no match
    """
    #If the mass > 1, indicating that it's input is in mg, divides it by a 1000 to get g
    parameter = "SAMPLE_SIZE"
    
    area_unit = headerValueCheck(header, parameter)

    if area_unit is None:
        return None
    
    area = area_unit[0]
    unit = area_unit[1]
    
    if unit == "mm2":
        area = area/100
        unit = "cm2"
        
    print(f'Sample size is {area:.5f} {unit} \n') 
    return area
    

# ei tea kas päris nii ikka teha
#Parsed thickness
def getThickness(data):
    """
    Checks if the datafile title (from the header) contains the sample thickness,
    usually not there but just in case does the control.

    Parameters
    ----------
    data : PANDAS DATAFRAME
        Measurement file header

    Returns
    -------
    float_val : FLOAT
        Thickness value
    None : Return None if no match

    """
    #Checks whether the title contains sample thickness in nm units: e.g. "25nm" and outputs 25
    try:
        thickness = data.iloc[3,1] #Title index in table
        pattern = r"(\d+)\s*(nm)"
        match = re.search(pattern,thickness)
        
        if match:
            float_str = match.group(1)
            print(f"Checking thickness: {float_str} nm")
            float_val = float(float_str)*10**-7 # nm to cm conversion
            print(f"Sample thickness is: {float_val} cm \n" )
            
            return float_val
        
        else:
            print("Sample thickness not found in title \n")
            return None
    except:
        print("Sample thickness unknown Exeption \n")
        return None
        


#---------------------------------------------------------------------------------------------------------------------------------------------------

def checkMeasurementType2(measurement_table, discrete_detection_ration = 0.02, min_count = 5):
    """
    Recognizes what measurement types are present in the file.
    
    Returns the temperatures and magnetic fields where the MvsH and MvsT measurements were made

    Parameters
    ----------
    measurement_table : PANDAS DATAFRAME
        Measurement file data.
    discrete_detection_ration : FLOAT, optional
        Discrete detection ratio. The default is 0.02.
    min_count : INT, optional
        Min amount of instances of a measurement to pass the level!?. The default is 5.

    Returns
    -------
    temperatures_of_interest : PANDAS SERIES
        MvsH measurement temperatures.
    magnetic_fields_of_interest : PANDAS SERIES
        MvsT measurement fields.

    """
    global tempCount
    #Checks what measurements the file contains
    
    rounded_dataset_T = measurement_table["Temperature (K)"].round(decimals=1)
    rounded_dataset_H = measurement_table["Magnetic Field (Oe)"]#.round(decimals=0)
    
    #returnitavad Seried
    magnetic_fields_of_interest = pd.Series(dtype=float)
    temperatures_of_interest = pd.Series(dtype=float)
    
    
    tempCount = rounded_dataset_T.value_counts()
    fieldCount = rounded_dataset_H.value_counts()
    
    # tempCount = tempCount[tempCount > min_count]
    # fieldCount = fieldCount[fieldCount > min_count]
    
    codes_T, uniques_T = pd.factorize(rounded_dataset_T)
    codes_H, uniques_H = pd.factorize(rounded_dataset_H)
    
    ratio_T = uniques_T.size/rounded_dataset_T.size
    ratio_H = uniques_H.size/rounded_dataset_H.size
    print(f"T : {ratio_T}, H : {ratio_H }")
    
    
    if ratio_T < discrete_detection_ration: #discrete
    
        if ratio_H < discrete_detection_ration: #discrete

            print("T discrete, H discrete = error \n")
        
        else: #continous
        
            print("T discrete, H continous = MvsH \n")
            temperatures_of_interest = pd.concat([temperatures_of_interest,pd.Series(tempCount.index.values)], ignore_index = True)
            
    else: #continous
    
        if ratio_H < discrete_detection_ration: #discrete

            print("T continous, H discrete = MvsT \n")
            fieldCount = fieldCount[fieldCount > min_count]
            magnetic_fields_of_interest = pd.concat([magnetic_fields_of_interest,pd.Series(fieldCount.index.values)], ignore_index = True)
        else: #continous

            print("T continous, H continous = both \n")
            
            meanTempCount=tempCount.mean()                     #see osa siin annab segafailide puhul mõistlikud numbrid. Aga mõne erilisema faili korral ei saa ta aru et sega fail on, aga need read annavad ka siis mõistlikud numbrid. Äkki saame kuidagi hoopis neid prinditavaid liste kasutada faili ära tundmiseks.
            muutuja = rounded_dataset_T.value_counts() > meanTempCount*10
            temperatures_of_interest = pd.Series(rounded_dataset_T.value_counts()[muutuja].index.values)
            print("meanTempCount",meanTempCount)
            meanFieldCount=fieldCount.mean()
            muutuja2 = rounded_dataset_H.value_counts() > meanFieldCount*10
            magnetic_fields_of_interest = pd.Series(rounded_dataset_H.value_counts()[muutuja2].index.values)
            
    return temperatures_of_interest, magnetic_fields_of_interest

#----------------------------------------MvsH specific functions-------------------------

def getMeasurementMvsH(const_T_values, bound = 1):
    """
    Saves the initial row indices of data point values that fall between the MvsH bound from the constant temperature points for further filtering.
    Returns them in a nested list.

    Parameters
    ----------
    const_T_values : PANDAS SERIES OF FLOATS
        Measurement temperatures.
    bound : FLOAT, optional
        The plus/min bound from the constant value to choose elements by. The default is 1.

    Returns
    -------
    row_indices : LIST OF LIST OF INT
        Nested list with int indices at each const temp bound.

    """
    #Saves all the indices of the points that fall between the bound
    table = ORIGINAL_DATAFRAME[['Temperature (K)', "color"]]
    filtered_dfs = []
    row_indices = []
    
    for value in const_T_values:
        lower = value - bound
        upper = value + bound
        filtered_df = table[(table['Temperature (K)'] >= lower) & (table['Temperature (K)'] <= upper) & (table['color'] == "black") ]
        indices = filtered_df.index.tolist()
        filtered_dfs.append(filtered_df)
        row_indices.append(indices)
        
    return row_indices #filtered_dfs


def removeBleedingElement(pairs):
    """
    This function ensures that there is no value "bleeding" from the next pair due to the way
    separationIndexForSingleSeries function slices multiple measurements made on the same const value.
    
    Special case when the first element from the next measurement "bleeds" into the previous measurement,
    the function uses a ratio = first_part_max/second_part_max to determine if the first and second part
    have the same max value and removes the second part last element if the ratio is out of bounds.
    
    Parameters
    ----------
    pairs : LIST OF LIST OF DATAFRAMES
        Nested list with dataframe pairs for the measurements.

    Returns
    -------
    new_pairs : LIST
        "pairs" list with removed elements, if there were any to remove.

    """

    new_pairs = []

    for pair in pairs:
        
        first = pair[0]
        
        second = pair[1]
        
        try:
            first_max = max(first["Magnetic Field (Oe)"])
            second_max = max(second["Magnetic Field (Oe)"])
            ratio = first_max/second_max
        except ValueError:
            error_message = "Change separationIndexForSingleSeries argument n for correct extrema points"
            showExtremaError(error_message)
        
        new_list = []

        while not 0.9 < ratio < 1.1:
            if ratio < 0.9:

                second = second[:-1]
                ratio = first_max/max(second["Magnetic Field (Oe)"])
                
            elif ratio > 1.1:

                second = second[:-1]
                break

        new_list.append(first)
        new_list.append(second)
        new_pairs.append(new_list)
        
    return new_pairs

def roundFieldForCorrection(pairs):
    """
    Returns the max field value for each individual MvsH measurement pair for correction fit

    Parameters
    ----------
    pairs : LIST OF LIST OF DATAFRAMES
        Nested list with dataframe pairs for the measurements.

    Returns
    -------
    values : LIST OF INT
        List with max field values for each measurement.

    """
    #rounds the magnetic field to a nice number (100041.221 -> 100000)
    values = []
    for pair in pairs:
        first_max = max(pair[0]["Magnetic Field (Oe)"])
        second_max = max(pair[1]["Magnetic Field (Oe)"])
        max_range = None
        
        # print(f"{first_max=}")
        # print(f"{second_max=}")
        if first_max > second_max:
            max_range = first_max
        elif first_max < second_max:
            max_range = second_max
            
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

def searchCorrectionTable(folder_path, number):
    """
    Searches for the correction table that is closest to the measurement max field value,
    returns the file path.

    Parameters
    ----------
    folder_path : STRING
        Folder path for the correction tables folder.
    number : INT
        Measurement max field value.

    Returns
    -------
            STRING
        File path for the closest correction table.

    """
    closest_match = None
    min_difference = float('inf')  # Initialize with positive infinity
    
    for filename in os.listdir(folder_path):
        parts = filename.split('_')  # Modify this based on your file naming convention
        try:
            value = int(parts[1])
            difference = abs(number - value)
            if difference < min_difference:
                min_difference = difference
                closest_match = filename
        except (ValueError, IndexError):
            # Ignore filenames that don't have a numeric part or invalid parts
            pass
        
    if closest_match:
        return os.path.join(folder_path, closest_match)
    else:
        return None

    
CORRECTION_FOLDER_PATH = os.path.join(USER_PATH,'PdCorrection tables')

def CorrectionTableToDict(correctionField_values):
    """
    Returns the corresponding amount of correction tables for each unique min/max field value
    in a dictionary form with the key being the value the table is for.

    Parameters
    ----------
    correctionField_values : INT
        Correction field max values.

    Returns
    -------
    correction_tables : DICT OF {int : Dataframe}
        Dictionary with max_field_value : values_at_that_field format.

    """
    #returns the corresponding amount of correction tables for each unique min/max measurement value in a dictionary form with the key being the value the table is for
    correction_tables = {}
    
    for nr in correctionField_values:
        correction_table_path = searchCorrectionTable(CORRECTION_FOLDER_PATH, nr)
        
        if correction_table_path is None:
            print(f"{nr} Oe jaoks ei leia parandustabelit ")
            continue
        correction_table_measurements = pd.read_csv(correction_table_path, index_col=0, encoding="ANSI")
        correction_tables[nr] = correction_table_measurements
        
    return correction_tables

def interpolateTrueField(magnetic_field_values, correction_tables):
    """
    Uses the correction table to interpolate the correct values for the measurement

    Parameters
    ----------
    magnetic_field_values : SERIES OF FLOAT
        One part of the measured MvsH field values.
    correction_tables : DATAFRAME
        Correction table in dataframe format.

    Returns
    -------
    true_field_interpolated : NUMPY.NDARRAY
        Interpolated correction values based on the measurement points.

    """
    x = correction_tables["Magnetic Field (Oe)"]
    y = correction_tables["True Field (Oe)"]

    # Create an interpolation function
    interp_func = interp1d(x, y, kind='linear', fill_value='extrapolate')

    # Call the interpolation function with magnetic_field_values to get the interpolated values
    true_field_interpolated = interp_func(magnetic_field_values)
    return true_field_interpolated

def interpolateMvsH(separated_MvsH, correction_tables):
    """
    Adds the true field values column to the existing measurement dataframes.

    Parameters
    ----------
    separated_MvsH : LIST OF LIST OF DATAFRAMES
        Nested list with MvsH measurements dataframes.
    correction_tables : DICT OF {int : Dataframe}
        Dict with the correction tables, field_value : values_at_that_field format.

    Returns
    -------
    interpolated_dict : LIST OF LIST OF DATAFRAMES
        "separated_MvsH" with the added true field column

    """

    interpolated_dict = {}
    for pair in separated_MvsH:
        
        # max_val = max(val_pair[0]["Magnetic Field (Oe)"]) # [0] ei pruugi alati õige max anda
        # print("max val",max_val)
        
        first_max = max(pair[0]["Magnetic Field (Oe)"])
        second_max = max(pair[1]["Magnetic Field (Oe)"])
        max_range = None
        
        # print(f"{first_max=}")
        # print(f"{second_max=}")
        if first_max > second_max:
            max_range = first_max
        elif first_max < second_max:
            max_range = second_max
            
        for key in correction_tables:
            if key - 200 <= max_range <= key + 200:

                #print(f"{max_range} kukub {key} vahemikku")
                
                for val in pair:
                    
                    magnetic_field_values = val["Magnetic Field (Oe)"]
                    # print(len(magnetic_field_values))
                    true_field_interpolated = interpolateTrueField(magnetic_field_values, correction_tables[key])
                    #print(type(true_field_interpolated))
                    # print(len(true_field_interpolated))
                    val["True Field (Oe)"] = true_field_interpolated
                                
    return interpolated_dict #!!1 siin vaata üle func returnib aga ei omista muutujale

def plotMvsH(dict_MvsH):
    """
    Plots the MvsH measurement pictures wtih a legend where it specifies the ascending and descending part with different shades.
    Plots the original field and true field values separately with different colors.
    
    It also saves the figures in the same folder the measurement file is in.

    Parameters
    ----------
    dict_MvsH : DICT OF {float : list of list of dataframes}
        Dictionary with MvsH measurements, temp_values : measurement_dataframes format.

    Returns
    -------
    Direct return None but saves the plots in a dedicated folder.

    """
    #Plots the MvsH measurement pair with different colors
    
    for key1 in dict_MvsH:
        
        i_pair = 1
        
        for df in dict_MvsH[key1]:
            
            colorIdx = df[0].iloc[1].name
            Color = ORIGINAL_DATAFRAME["color"].loc[colorIdx]
            
            fig, ax = plt.subplots()
            H1 = df[0]["Magnetic Field (Oe)"]
            M1 = df[0]["Moment (emu)"]
            H2 = df[1]["Magnetic Field (Oe)"]
            M2 = df[1]["Moment (emu)"] 
            
            ax.plot(H1 ,M1, color = "grey", label = "Descending")#!!! mis siin värvidega jääb
            ax.plot(H2, M2, color = "grey", label = "Ascending")
            
            if "True Field (Oe)" in df[0]:
                
                H1_true = df[0]["True Field (Oe)"]
                H2_true = df[1]["True Field (Oe)"]   
                ax.plot(H1_true, M1, color = Color, label = "True Field Descending", alpha = 0.5)
                ax.plot(H2_true, M2, color = Color, label = "True Field Ascending")
            
            const_val = int(key1)
            
            plot_title = f"M vs H at {const_val} K"
            ax.set_title(plot_title)
            ax.set_xlabel("Magnetic field (Oe)")
            ax.set_ylabel("Moment (emu)")
            ax.legend() #Hetkel legend nimetab selle järgi et esimene tsükkel on kasvav ja teine kahanev ehk eeldus et mõõtmisel temp algas kasvamisest
            ax.grid(True)
            #fig.savefig(f"C:/Users/kevin/OneDrive/Desktop/Andmetöötlus/Projekt_andmed1/MvsH_graph_at_{val}K.png",bbox_inches = "tight", dpi = 200) #laptop
            fig.savefig(os.path.join(folder_name,f'MvsH_graph_at_{const_val}K_{i_pair}.png'),bbox_inches = "tight", dpi = 200) #PC
            i_pair = i_pair + 1
            plt.show()
        
    return None
#-----------------------------MvsT specific functions--------------------------------

def getMeasurementMvsT(const):
    """
    Saves the initial row indices of the MvsT measurement points that are equal to the const field value.

    Parameters
    ----------
    const : LIST OF FLOAT
        Measurement temperatures as floats.

    Returns
    -------
    row_indices : LIST OF LIST OF INT
        Nested list with row indices at each const field

    """
    #Saves all the indices of the points that are equal to the predetermined H value
    row_indices = {}
    table = ORIGINAL_DATAFRAME['Magnetic Field (Oe)']

    indices = table.index[table == const].tolist()
    row_indices[const] = indices

    return row_indices

def plotMvsT(dict_MvsT):
    """
    Plots the MvsT measurement pictures wtih a legend where it specifies the ascending and descending part with different shades.
    
    If there are multiple measurements on the same field, plots them on the same graph.
    It also saves the figures in the same folder the measurement file is in.
    
    Parameters
    ----------
    dict_MvsT : DICT OF {float : list of list of dataframes}
        Dictionary with MvsT measurements, field_value : measurement_dataframes format.

    Returns
    -------
    Direct return None but saves the plots in a dedicated folder.

    """
   
    for key in dict_MvsT:
        fig, ax = plt.subplots()
        i_pair = 1
        for df in dict_MvsT[key]:

            
            T1 = df[0]["Temperature (K)"]
            M1 = df[0]["Moment (emu)"]
            T2 = df[1]["Temperature (K)"] if len(df) > 1 else None
            M2 = df[1]["Moment (emu)"] if len(df) > 1 else None
            
            colorIdx = df[0].iloc[0].name
            Color = ORIGINAL_DATAFRAME["color"].loc[colorIdx]
            
            ax.plot(T1,M1,color = Color, label = f"Ascending {i_pair}", alpha = 0.5) # peaks tegelt kontrollima kas kasvab või kahaneb
            ax.plot(T2,M2,color = Color, label = f"Descending {i_pair}") if len(df) > 1 else None #, marker = "o") #descending ei pea paika kui on alt üle > alt üles mõõtmine
            ax.set_title(f"M vs T at {key} Oe")
            ax.set_xlabel("Temperature (K)")
            ax.set_ylabel("Moment (emu)")
            ax.legend() #Hetkel legend nimetab selle järgi et esimene tsükkel on kasvav ja teine kahanev ehk eeldus et mõõtmisel temp algas kasvamisest
            ax.grid(True)
            i_pair = i_pair + 1
        fig.savefig(os.path.join(folder_name,f'MvsT_graph_at_{key}K.png'),bbox_inches = "tight", dpi = 200)
        plt.show()
        
    return None

#---------------------------------------Universal functions that both paths use----------------------------------------------

#Separates all the series in the indices that grow by one and returns the longest one, that one being the measurement indices
def filterMeasurementIndices(unfiltered_indices):
    """
    Filters the measurement indices so that the longest consecutive list is the measurement.

    Parameters
    ----------
    unfiltered_indices : LIST OF LIST OF INT
        Nested list with unfiltered measurement data indices.

    Returns
    -------
    filtered : LIST OF LIST OF INT
        Nested list with filtered measurement data indices.

    """
    filtered = []
    
    for unfiltered in unfiltered_indices:
        
        consecutive_sequences = []
        current_sequence = [unfiltered[0]]
    
        for i in range(1, len(unfiltered)):
            if unfiltered[i] - unfiltered[i - 1] == 1:
                current_sequence.append(unfiltered[i])
            else:
                if len(current_sequence) > 1:
                    consecutive_sequences.append(current_sequence)
                current_sequence = [unfiltered[i]]
    
        # Check if the last sequence is consecutive and has more than one element
        if len(current_sequence) > 1:
            consecutive_sequences.append(current_sequence)
        
        longest_sequence = max(consecutive_sequences, key=len, default=[])
        filtered.append(longest_sequence)
    
    return filtered


#Returns the separation indices for ascending and descending points based on the extrema
def separationIndexForSingleSeries(data, column_name, x): #!!! https://stackoverflow.com/questions/48023982/pandas-finding-local-max-and-min
    """
    Returns the indices of the series local peaks (maxima and minima) in a list

    Parameters
    ----------
    data : SERIES OF FLOAT.
        The series to analyze.
    column_name : STRING
        Name of the column to analyze based on the measurement.
    x : FLOAT
        Percentage parameter, determines the percentage of n (number of points to compare around the extrema) to use. 

    Returns
    -------
    min_indices : PANDAS.INDEX
        Series minima index (int) values 
    max_indices : PANDAS.INDEX
        Series maxima index (int) values
    """
    
    if isinstance(data, pd.Series):
        # Convert a Series to a DataFrame with a specified column name
        data = pd.DataFrame({column_name: data})
        
    index = data.index
    
    n = int(x*len(data))
    #print(f"n for this run: {n}")
    
    # Find local peaks
    relative_min_indices = argrelextrema(data[column_name].values, np.less_equal, order=n)[0]
    relative_max_indices = argrelextrema(data[column_name].values, np.greater_equal, order=n)[0]
    
    min_indices = index[relative_min_indices]
    max_indices = index[relative_max_indices]
    
    title = ""
    
    def removeSpecialCaseIndex():
        """
        Helper function that removes the first min/max index if it's value is near 0 tesla
        because it is not needed from the data perspective. All the points from the first removed
        index till the next index will be removed from the data plot due to this.
        """
        nonlocal data, min_indices, max_indices, title

        if column_name == "Magnetic Field (Oe)":
            
            if not len(max_indices) == 1:
                title = "(Esimene min indeks eemaldatud)"
                first_index = min_indices[0]
                first_val = data.loc[first_index, column_name]
                # print(f"{first_index=}")
                # print(f"{first_val=}")
                if -1 < first_val < 1:
                    min_indices = min_indices[1:]
                
    removeSpecialCaseIndex()
    
    # Create a DataFrame to store results
    local_peaks = pd.DataFrame(index=index)
    local_peaks['min'] = np.nan
    local_peaks['max'] = np.nan
    
    # Fill in the 'min' and 'max' columns with peak values 
    local_peaks.loc[min_indices, 'min'] = data.loc[min_indices, column_name]
    local_peaks.loc[max_indices, 'max'] = data.loc[max_indices, column_name]
    
    # Extract max min indices
    mask = local_peaks.notna()
    max_indices = local_peaks["max"].index[mask["max"]]
    min_indices = local_peaks["min"].index[mask["min"]]
    
    # Plot results, tegelikult pole vaja, aga hea kontrollida kas tegi õigesti
    plt.scatter(index, local_peaks['min'], c='r', label='Minima')
    plt.scatter(index, local_peaks['max'], c='g', label='Maxima')
    plt.plot(index, data[column_name], label=column_name)
    plt.title(f"Extrema graph \n {title}")
    plt.legend()
    plt.show()
        
    return min_indices, max_indices

# Iterates the functionality of separationIndexForSingleSeries 
def separationIndexForMultipleSeries(indices, column_name, x=0.1):
    """
    Iterates the separationIndexForSingleSeries function over all the separate indices and appends it's results.
    
    This is used for separating single measurements made on the same const value into two ascending/descending pair.

    Parameters
    ----------
    indices : LIST OF LIST OF INT
        Nested list with measurement indices.
    column_name : STRING
        Name of the column to analyze based on the measurement.

    Returns
    -------
    indices_for_separation : LIST OF TUPLE OF (pandas.Index,pandas.Index)
        List with the indices (int) for the separation in a tuple.

    """
    indices_for_separation = []
    
    for measurement_index in indices:
        
        measurement = ORIGINAL_DATAFRAME.loc[measurement_index, column_name]
        separation_indices = separationIndexForSingleSeries(measurement, column_name,x)
        indices_for_separation.append(separation_indices)
        
    return indices_for_separation

def separateMeasurementWithColorIdx(indices_for_separation, measurement_indices, column_name):
    """
    Separates the points based on the separation indices and returns the separated series in pairs

    Parameters
    ----------
    indices_for_separation : LIST OF TUPLE OF (pandas.Index,pandas.Index)
        List with the indices (int) for the separation in a tuple..
    measurement_indices : LIST OF LIST OF INT
        Nested list with measurement indices.
    column_name : STRING
        Name of the column to use based on the measurement.

    Returns
    -------
    separated_pair : LIST OF LIST OF DATAFRAME
        DESCRIPTION.
    pair_indices : TYPE
        DESCRIPTION.

    """
    
    separated_pair_all = []
    pair_indices = []
    
    if not isinstance(indices_for_separation, list):
        indices_for_separation = [indices_for_separation]
    
    color_index = 0
    
    for min_max_index in indices_for_separation: 
        
        if min_max_index[0][0] < min_max_index[1][0]: #First assigns the correct indices, assumes data is shaped as + - + or - + -
            min_index_list = min_max_index[0].tolist()
            max_index_list = min_max_index[1].tolist()
        else:
            min_index_list = min_max_index[1].tolist()
            max_index_list = min_max_index[0].tolist()
    
        iteration_index = 0
        
        for indices in measurement_indices: #Iterates through the separate indices and separates them
    
            measurement = ORIGINAL_DATAFRAME[[column_name,"Moment (emu)"]].loc[indices]
        
            for max_index in max_index_list: #Slices them from min to max and from max to min pairs
                
                separated_pair = []
                indices_pair1 = []
                indices_pair2 = []
                
                if max_index in indices: #Checks if 
                    
                    sliced1 = measurement.loc[min_index_list[iteration_index]:max_index] #paaride data
                    separated_pair.append(sliced1)
                    
                    ORIGINAL_DATAFRAME.loc[min_index_list[iteration_index]:max_index, "color"] = COLORS[color_index] #värvid paaridele
                    
                    indices_pair1 = ORIGINAL_DATAFRAME.loc[min_index_list[iteration_index]:max_index].index.tolist() #paaride indeksid
                    
                    if iteration_index == len(min_index_list) - 1:
                        sliced2 = measurement.loc[max_index+1:min_index_list[iteration_index]]
                        
                        ORIGINAL_DATAFRAME.loc[max_index+1:min_index_list[iteration_index], "color"] = COLORS[color_index]
                        
                        indices_pair2 = ORIGINAL_DATAFRAME.loc[max_index+1:min_index_list[iteration_index]].index.tolist()
                        
                    else:
                        sliced2 = measurement.loc[max_index+1:min_index_list[iteration_index+1]]
                        
                        ORIGINAL_DATAFRAME.loc[max_index+1:min_index_list[iteration_index+1], "color"] = COLORS[color_index]
                        
                        indices_pair2 = ORIGINAL_DATAFRAME.loc[max_index+1:min_index_list[iteration_index+1]].index.tolist()
                else:
                    max_index_list = max_index_list[iteration_index:]

                    break
                
                pair_indices.append(indices_pair1 + indices_pair2)
                separated_pair.append(sliced2)
                separated_pair_all.append(separated_pair)
                
                color_index += 1
                iteration_index += 1
                if color_index == len(COLORS): #if the colors are used up it goes back to the beginning
                    color_index = 0
                    
    return separated_pair_all, pair_indices

def separateIntoDictValuePair(separated_pairs, const_val, column_name, token):
    # Creates a dict{const value the measurement was made at: measurement} 
    raamat = {}
    
    for const in const_val:
        
        for val in separated_pairs:
            
            index_to_check = val[0].index[0]
            val_to_check = ORIGINAL_DATAFRAME[column_name].iloc[index_to_check]
            
            #Rounds the values if the measurement needs to check for MvsH, for MvsT can use direct == becasue they are precise
            val_to_check, const = (round(val_to_check), round(const)) if token == "MvsH" else (val_to_check, const)
            
            if val_to_check == const: #or round(val_to_check) == round(const):
                
                # print("Check", val_to_check, "Rounded", round(val_to_check))
                # print("const", const, "rounded", round(const))
                key = const
                
                if key not in raamat:
                    raamat[key] = []  # Create an empty list for this key if it doesn't exist
                    
                raamat[key].append(val)  # Append the value to the list associated with the key

    return raamat

def plotMeasurementTimeseries(): 

    # Create subplots with shared x-axis
    fig, axes = plt.subplots(nrows=3, sharex=True)
    
    # Plot data on each subplot
    time_axis = "Time Stamp (sec)"
    ORIGINAL_DATAFRAME.plot(x=time_axis, y="Temperature (K)", kind="scatter", c=ORIGINAL_DATAFRAME["color"].values, ax=axes[0])
    ORIGINAL_DATAFRAME.plot(x=time_axis, y="Magnetic Field (Oe)", kind="scatter", c=ORIGINAL_DATAFRAME["color"].values, ax=axes[1])
    ORIGINAL_DATAFRAME.plot(x=time_axis, y="Moment (emu)", kind="scatter", c=ORIGINAL_DATAFRAME["color"].values, ax=axes[2])
    
    # Customize axes labels and other properties
    
    axes[0].set_ylabel("Temperature (K)", fontsize = 8)
    axes[1].set_ylabel("Magnetic Field (Oe)", fontsize = 8)
    axes[1].tick_params(axis = "y", labelsize = 8)
    axes[2].set_ylabel("Moment (emu)", fontsize = 8)
    axes[-1].set_xlabel("Time Stamp (sec)")
    fig.suptitle("Timeseries \n (black dots are discarded from individual plots)", fontsize = 10)
    
    fig.savefig(os.path.join(folder_name,'timeseries.png'),bbox_inches = "tight", dpi = 200)
    
    # Show the plot
    plt.show()
    
    return None


def setColumnForType(indices, type_string):
    
    for idx in indices:
        
        ORIGINAL_DATAFRAME.loc[idx, "Type"] = type_string
        
    return None


def addParameterColumns(separated, type_string):
    
    temp = "Temperature (K)"
    temp_unit = "K"
    
    field = "Magnetic Field (Oe)"
    field_unit = "Oe"
    
    moment = "Moment (emu)"
    moment_unit = "emu"
    error = "M. Std. Err. (emu)"
    
    momentDivMass = "Moment (emu)/mass(g)"
    momentDivMass_unit = "emu/g"
    
    momentDivArea = "Moment (emu)/area(cm2)"
    momentDivArea_unit = "emu/cm2"
    
    momentDivVolume = "Moment (emu)/volume(cm3)"
    momentDivVolume_unit = "emu/cm3"
    
    susceptibility = "Susceptibility (emu/g Oe)"
    susceptibility_unit = "emu/g*Oe"
    
    oneOverSusceptibility = "1/Susceptibility"
    oneOverSusceptibility_unit = "g*Oe/emu"
    
    volume = SAMPLE_AREA_CM2*THICKNESS if SAMPLE_AREA_CM2 and THICKNESS else None

    unit_row = None
    
    for i ,pair in enumerate(separated):
        
        for j, series in enumerate(pair):
            
            if type_string == "MvsH":
                
                indices = series.index
                series[temp] = ORIGINAL_DATAFRAME.loc[indices, temp]
                unit_row = pd.DataFrame({ field: [field_unit], moment: [moment_unit], "True Field (Oe)": [field_unit], temp: [temp_unit],
                                        error: [moment_unit], momentDivMass: [momentDivMass_unit], momentDivArea: [momentDivArea_unit],
                                        momentDivVolume: [momentDivVolume_unit],susceptibility: [susceptibility_unit], oneOverSusceptibility: [oneOverSusceptibility_unit] }, index=['unit'])
                
            elif type_string == "MvsT":
                
                indices = series.index
                series[field] = ORIGINAL_DATAFRAME.loc[indices, field]
                unit_row = pd.DataFrame({ temp: [temp_unit], moment: [moment_unit], field: [field_unit],
                                        error: [moment_unit], momentDivMass: [momentDivMass_unit], momentDivArea: [momentDivArea_unit],
                                        momentDivVolume: [momentDivVolume_unit],susceptibility: [susceptibility_unit], oneOverSusceptibility: [oneOverSusceptibility_unit] }, index=['unit'])
            
            series[ error] = ORIGINAL_DATAFRAME.loc[indices, error]
            series[momentDivMass] = series[moment]/SAMPLE_MASS_g if SAMPLE_MASS_g else None
            series[ momentDivArea] = series[moment]/SAMPLE_AREA_CM2 if SAMPLE_AREA_CM2 else None
            series[ momentDivVolume] = series[moment]/volume if volume else None
            series[susceptibility] = series[moment]/(SAMPLE_MASS_g*series[field]) if SAMPLE_MASS_g else None
            series[oneOverSusceptibility] = 1/(series[moment]/(SAMPLE_MASS_g*series[field])) if SAMPLE_MASS_g else None
            
            # Concatenate the new row DataFrame and the original DataFrame
            pair[j] = pd.concat([unit_row, series])
            
        separated[i] = pair
        
    return None

# Siin on ka üks huvitav error kui file excelis avatud sama aeg siis ei luba uut üle salvestada
def appendAndSave(dictionary, dType):
    i_key = 1
    
    for key in dictionary:

        i_key = i_key + 1
        i_pair = 1
        for pair in dictionary[key]:
            
            result = pd.concat([pair[0], pair[1].tail(pair[1].shape[0]-1)])
            
            file_name = f'{dType}_data_at_{key}_{i_pair}.csv'
            
            full_path = os.path.join(folder_name, file_name)
            
            result.to_csv(full_path, index = False)
            
            i_pair = i_pair + 1
            
    return None

def showExtremaError(message):
    
    root = tk.Tk()
    root.withdraw()  # Hide the main window

    messagebox.showerror("Error", message)
    
#---------------Measurement errors----------------------------

def BTypeMeasurementError(measurement_deviation):
    StudentFactor = 2 #Inf, reliability = 0.95
    formula = (StudentFactor/3)*measurement_deviation
    return formula

def momentDivDimensionUncertaintyError(separated, dimension, measurement_deviation):
    
    for pair in separated:
        
        for df in pair:
            
            momentStdError = df["M. Std. Err. (emu)"].iloc[1:]
            moment = df["Moment (emu)"].iloc[1:]
            # momentDivMassError = math.sqrt(((1/SAMPLE_MASS_g)*momentStdError)**2+((-moment/SAMPLE_MASS_g**2)*momentStdError)**2)
            momentDivMassError = (((1/dimension)*momentStdError)**2 + 
                                  ((-moment/dimension**2)*measurement_deviation)**2)**0.5
            df["Error test"] = momentDivMassError
    return None

#-------------- Actually Run the program here -------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

# Ask user input to get a path to datafile
DATAFILE_PATH = askNewDatafilePath()
save_to_path = os.path.dirname(DATAFILE_PATH)

# Read the datafile
HEADER, ORIGINAL_DATAFRAME = readDatafile(DATAFILE_PATH)
ORIGINAL_COPY = copy.deepcopy(ORIGINAL_DATAFRAME)
# VSM or ACMSII? or maybe HC or ETO in the future
OPTION_TYPE = determineDatafileType(HEADER)
#print(OPTION_TYPE) kas teha nii et funk ise ei prindi ja prindid kui tahad või funktsioon prindib anyway?


#Selle lisasin juurde kuna moment tulbas võib olla nan values ja enne pead kõik õiged tulbad võtma, et need eraldada, muidu eemaldab kõik read,
# sest igas reas on mingi tulp nan value'ga
if OPTION_TYPE == "VSM":
    ORIGINAL_DATAFRAME = ORIGINAL_DATAFRAME[["Time Stamp (sec)", "Temperature (K)", "Magnetic Field (Oe)", "Moment (emu)", "M. Std. Err. (emu)"]].dropna()
    ORIGINAL_DATAFRAME = ORIGINAL_DATAFRAME.reset_index(drop = True)
    
elif OPTION_TYPE == "ACMS": #!!! siia veel eraldi check et kas AC või DC
    ORIGINAL_DATAFRAME = ORIGINAL_DATAFRAME[["Time Stamp (sec)", "Temperature (K)", "Magnetic Field (Oe)", "DC Moment (emu)", "DC Std. Err. (emu)"]].dropna()
    ORIGINAL_DATAFRAME = ORIGINAL_DATAFRAME.rename(columns={"DC Moment (emu)": "Moment (emu)", "DC Std. Err. (emu)": "M. Std. Err. (emu)"})
    ORIGINAL_DATAFRAME = ORIGINAL_DATAFRAME.reset_index(drop = True)

#parse HEADER
SAMPLE_MASS_g = getMassInGrams(HEADER)
SAMPLE_AREA_CM2 = getAreaCM2(HEADER)
#siin peaks olema try ja catch exepction ümber
THICKNESS = getThickness(HEADER)


#Color list for color idx and initializing the starting color idx with just black
COLORS = ["red", "green", "blue", "yellow", "brown", "purple", "orange", "pink", "olive", "magenta"]
ORIGINAL_DATAFRAME["color"] = "black"

print("_________chechMeasurementType2-----------")  #mis peaks olema selle funktsiooni ebaõnnestumise/veateade? return None? või lihtsalt kaks tühja Seriet?
# Tagastab kaks Pandas.Seriest: temperatures_of_interest, magnetic_fields_of_interest
# edasi peaks veel kontrollima välja filtreerima üksikud punktid mis ei kuulu MvsH ja MvsT andmete hulka
TEMPERATURES_OF_INTEREST, MAGNETIC_FIELDS_OF_INTEREST= checkMeasurementType2(ORIGINAL_DATAFRAME)
#print(TEMPERATURES_OF_INTEREST, MAGNETIC_FIELDS_OF_INTEREST)
print("_________end-----------")


print('--------<<<<<<<<<>>>>>>>>>>-----------')
print('--------<<<<<<<<<>>>>>>>>>>-----------')

#creates a column "Type" for each data point type
ORIGINAL_DATAFRAME["Type"] = ""

#Creates a folder for the current data file to save related files
folder_name = os.path.splitext(DATAFILE_PATH)[0] + ""
os.makedirs(folder_name, exist_ok = True)

if MAGNETIC_FIELDS_OF_INTEREST.size <= 0:
    print('no MvsT detected')
    
else:
    print(' MvsT data detected')
    print(MAGNETIC_FIELDS_OF_INTEREST)
    
    unfiltered_MvsT_indices = getMeasurementMvsT(MAGNETIC_FIELDS_OF_INTEREST) # SELLE SAAB ÜHE INPUTIGA

    def cycle(const):#TEE KÕIK FUNKTSIOONID NIIMOODI ÜMBER ET TERVE TSÜKKEL TEEB KÕIK MÖÖTMISED ÜHEL CONST VÄÄRTUSEL ALGUSEST LÕPUNI
        unfiltered_MvsT_indices = getMeasurementMvsT(const)
        MvsT_INDICES = filterMeasurementIndices(unfiltered_MvsT_indices)
        separation_index_MvsT = separationIndexForMultipleSeries(MvsT_INDICES, "Temperature (K)")
        return None

    for const in MAGNETIC_FIELDS_OF_INTEREST:
        try:
            cycle(const)
        except:
            #mingi indikaator näiteks timeseries et need punktid feilisid
            pass
    MvsT_INDICES = filterMeasurementIndices(unfiltered_MvsT_indices)
    test_MvsT = copy.deepcopy(MvsT_INDICES)
    #loop list Xidest (0.01, 0.1, 0.5)
    separation_index_MvsT = separationIndexForMultipleSeries(MvsT_INDICES, "Temperature (K)")# the indices where the separation is going to be done
    
    try:#siin võib juhtuda et liiga väike n siis tulevad valesti ekstreemumid, siis custom error
        SEPARATED_MvsT, MvsT_pair_indices = separateMeasurementWithColorIdx(separation_index_MvsT, MvsT_INDICES, "Temperature (K)")
        DICT_MvsT = separateIntoDictValuePair(SEPARATED_MvsT, MAGNETIC_FIELDS_OF_INTEREST, "Magnetic Field (Oe)", "MvsT")
        
    except IndexError:
        print('temp sõltuvuse error')
        # raise ValueError("separationIndexForSingleSeries funktsiooni n argumenti peab muutma, ekstreemumid tulevad valesti sellise n puhul") from ie
        separation_index_MvsT = separationIndexForMultipleSeries(MvsT_INDICES, "Temperature (K)", x=0.5)# the indices where the separation is going to be done
        SEPARATED_MvsT, MvsT_pair_indices = separateMeasurementWithColorIdx(separation_index_MvsT, MvsT_INDICES, "Temperature (K)")
        DICT_MvsT = separateIntoDictValuePair(SEPARATED_MvsT, MAGNETIC_FIELDS_OF_INTEREST, "Magnetic Field (Oe)", "MvsT")
        

    plotMvsT(DICT_MvsT)
    
    setColumnForType(MvsT_INDICES, "MvsT")
    addParameterColumns(SEPARATED_MvsT, "MvsT")#this function modifies SEPARATED_MvsT which inturn modifies DICT_MvsT since it's a global mutable variable
    appendAndSave(DICT_MvsT, "MvsT")
    
print('--------<<<<<<<<<>>>>>>>>>>-----------')
print('--------<<<<<<<<<>>>>>>>>>>-----------')

#!!!
if TEMPERATURES_OF_INTEREST.size <= 0:
    print('no MvsH detected')
     
    
else:
    print(' MvsH data detected')
    print(TEMPERATURES_OF_INTEREST)
    
    unfiltered_MvsH_INDICES = getMeasurementMvsH(TEMPERATURES_OF_INTEREST)
    MvsH_INDICES = filterMeasurementIndices(unfiltered_MvsH_INDICES)
    
    separation_indices_MvsH = separationIndexForMultipleSeries(MvsH_INDICES, "Magnetic Field (Oe)")
    # test_MvsH1 = copy.deepcopy(MvsH_INDICES)
    #try:#siin võib juhtuda et liiga väike n siis tulevad valesti ekstreemumid, siis custom error
    SEPARATED_MvsH, MvsH_pair_indices = separateMeasurementWithColorIdx(separation_indices_MvsH, MvsH_INDICES, "Magnetic Field (Oe)")
    # except IndexError as ie:
    #     raise ValueError("separationIndexForSingleSeries funktsiooni n argumenti peab muutma, ekstreemumid tulevad valesti sellise n puhul") from ie
    
    #test_MvsH2 = copy.deepcopy(SEPARATED_MvsH)
    
#    SEPARATED_MvsH = removeBleedingElement(SEPARATED_MvsH)
    try:
        correction_field_value = roundFieldForCorrection(SEPARATED_MvsH)
        CORRECTION_TABLES = CorrectionTableToDict(correction_field_value)
        
        interpolateMvsH(SEPARATED_MvsH, CORRECTION_TABLES)
        
        DICT_MvsH = separateIntoDictValuePair(SEPARATED_MvsH, TEMPERATURES_OF_INTEREST, "Temperature (K)", "MvsH")
        
        plotMvsH(DICT_MvsH)
        
        setColumnForType(MvsH_INDICES, "MvsH")
        addParameterColumns(SEPARATED_MvsH, "MvsH")#this function modifies SEPARATED_MvsH which inturn modifies DICT_MvsH since it's a global mutable variable
        appendAndSave(DICT_MvsH, "MvsH")
    except IndexError:
       print("välja sõltuvuse error")

print('--------<<<<<<<<<>>>>>>>>>>-----------')
print('--------<<<<<<<<<>>>>>>>>>>-----------')


if MAGNETIC_FIELDS_OF_INTEREST.size <= 0 and TEMPERATURES_OF_INTEREST.size <= 0:
    print('Error, ei suutnud eraldada MvsH ja MvsT mõõtmisi')


#Plots temp, field and moment against time
plotMeasurementTimeseries()


#Error       
# momentDivDimensionUncertaintyError(SEPARATED_MvsH, SAMPLE_MASS_g, 0.0001) #for moment/mass uncertainty



