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
from tkinter import filedialog
from scipy.interpolate import interp1d
from scipy.signal import argrelextrema


USER_PATH = os.getcwd()


#-------------- OPENING THE FILE AND INDEXING IT -------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------
root = tk.Tk()
root.wm_attributes('-topmost', 1)
root.withdraw()


#save_to_path = os.path.dirname(file_path)

#ask for user input for a datafile        
def askNewDatafilePath():
    
    file_path = filedialog.askopenfilename()
    
    return file_path


#Read datafile form file_path, return header and data as pd df
def readDatafile(file_path):
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



# determine if its VSM or ACMS datafile, return "token"
# Headers of VSM and ACMS files are similiar. DATA columns of those files have a difference in the Moment column. IN VSM the columns is named Moment (emu), while in ACMS its named DC Moment (emu) 
def determineDatafileType(header):
    #data_type checks from the first line under the header whether the file is VSM or ACMS and returns a token for further use
    token = "error - unknown datafile format"
    
    option_specific_line = header['1'][1]
    
    
    if "VSM" in option_specific_line:
        print("\nThis is a VSM data file \n")
        token = "VSM"
        
    elif "ACMS" in option_specific_line:
        print("\nThis is an ACMS data file \n")
        token = "ACMS"
        
    return token

        
# text parsing function to help with sample parameters
def extractFloatWithUnit(string):
    #extract_float_with_unit function uses regex to index the input string
    #into a (float, unit) format if it has units, if no units (float, None) format, if no value (None)
    regex = r'^([\d.]+)\s*([a-zA-Z]+(\d)?)?$'# regular expression to match float value and unit
    match = re.search(regex, string)
    
    if match:
        float_str = match.group(1)
        unit = match.group(2)
        float_val = float(float_str)
        print(f"SAMPLE MASS :{float_val}, {unit}")
        return (float_val, unit)
    else:
        return None
    
    
# text parsing function before   extractFloatWithUnit function, returns Header property value and unit
def headerValueCheck(header, sample_property):
    #header_value_check checks the value of interest, is the value nan or is it just a string in 
    #the wrong format, otherwise outputs the float value
    header_column = header['2']

    print(f"Checking:{sample_property}, {header_column[sample_property]}")#Hetkel jääb
    
    if not isinstance(header_column[sample_property], str) and np.isnan(header_column[sample_property]): #Checks whether the value is not a string and if it is a nan value
        print(f"NO VALID VALUE FOR {sample_property}, value is NAN \n")
        return None
    
    match = extractFloatWithUnit(header_column[sample_property])
    
    if match is None: #condition for extract_float_with_unit when it didn't find a match for a desired format therefore being a string
        print(f"{sample_property} had an input put it is not in the valid format: {header_column[sample_property]} \n")
        return None
    
    float_val = match[0]
    #print("float value:", float_val)#Hetkel jääb
    unit = match[1]
    #print("units:", unit)#Hetkel jääb
    
    return float_val, unit #Ühikutega peaks veel midagi tegema


#parse HEADER, return mass, area, volume
def parseHeader(header):
    return "NaN"
    

#returned parsed MASS from Header
def getMass(header):
    parameter = "SAMPLE_MASS"
    
    mass = headerValueCheck(header, parameter)

    return mass


def getMassInGrams(header):
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

def getAreaCM2(header):
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
def getThickness(data):
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
        return None





#---------------------------------------------------------------------------------------------------------------------------------------------------

# Tuvastab kas sisaldab MvsH ja MvsT mõõtmisi. Tagastab kaks Pandas.Seriest: temperatures_of_interest, magnetic_fields_of_interest
# Mida peaks tagastama errori korral?

def checkMeasurementType2(measurement_table, discrete_detection_ration = 0.02):
    #Checks what measurements the file contains
    
    #global uniques_T,uniques_H,codes_T
    
    type_token = {"Temperature": None, "Field": None}
    
    rounded_dataset_T = measurement_table["Temperature (K)"].round(decimals=0)
    rounded_dataset_H = measurement_table["Magnetic Field (Oe)"]#.round(decimals=1)
    
    #returnitavad Seried
    magnetic_fields_of_interest = pd.Series(dtype=float)
    temperatures_of_interest = pd.Series(dtype=float)
    
    
    tempCount = rounded_dataset_T.value_counts()
    fieldCount = rounded_dataset_H.value_counts()
    
    
    codes_T, uniques_T = pd.factorize(rounded_dataset_T)
    codes_H, uniques_H = pd.factorize(rounded_dataset_H)
    
    ratio_T = uniques_T.size/rounded_dataset_T.size
    ratio_H = uniques_H.size/rounded_dataset_H.size
    print(f"T : {ratio_T}, H : {ratio_H }")

    
    if ratio_T < discrete_detection_ration: #discrete
    
        if ratio_H < discrete_detection_ration: #discrete
        
            type_token["Temperature"] = "discrete"
            type_token["Field"] = "discrete"
            print("T discrete, H discrete = error \n")
        
        else: #continous
        
            type_token["Temperature"] = "discrete"
            type_token["Field"] = "continous"
            print("T discrete, H continous = MvsH \n")
            temperatures_of_interest = pd.concat([temperatures_of_interest,pd.Series(tempCount.index.values)], ignore_index = True)
            
    else: #continous
    
        if ratio_H < discrete_detection_ration: #discrete
        
            type_token["Temperature"] = "continous"
            type_token["Field"] = "discrete"
            print("T continous, H discrete = MvsT \n")
            magnetic_fields_of_interest = pd.concat([magnetic_fields_of_interest,pd.Series(fieldCount.index.values)], ignore_index = True)
        
        else: #continous
        
            type_token["Temperature"] = "continous"
            type_token["Field"] = "continous"
            print("T continous, H continous = both \n")
            
            meanTempCount=tempCount.mean()                     #see osa siin annab segafailide puhul mõistlikud numbrid. Aga mõne erilisema faili korral ei saa ta aru et sega fail on, aga need read annavad ka siis mõistlikud numbrid. Äkki saame kuidagi hoopis neid prinditavaid liste kasutada faili ära tundmiseks.
            muutuja = rounded_dataset_T.value_counts() > meanTempCount*10
            temperatures_of_interest = pd.Series(rounded_dataset_T.value_counts()[muutuja].index.values)
            
            meanFieldCount=fieldCount.mean()
            muutuja2 = rounded_dataset_H.value_counts() > meanFieldCount*10
            magnetic_fields_of_interest = pd.Series(rounded_dataset_H.value_counts()[muutuja2].index.values)
            
    return type_token, temperatures_of_interest, magnetic_fields_of_interest

#----------------------------------------MvsH specific functions-------------------------

#const väärtuste põhjal mõõtmise punktid võtab
def getMeasurementMvsH(const_T_values, bound = 0.15):
    #Saves all the indices of the points that fall between the bound
    table = ORIGINAL_DATAFRAME[['Temperature (K)', "color"]]
    filtered_dfs = []
    all_indices = []
    
    for value in const_T_values:
        lower = value - bound
        upper = value + bound
        filtered_df = table[(table['Temperature (K)'] >= lower) & (table['Temperature (K)'] <= upper) & (table['color'] == "black") ]
        indices = filtered_df.index.tolist()
        filtered_dfs.append(filtered_df)
        all_indices.append(indices)
        
    return all_indices #filtered_dfs

#ümardab iga mõõtmise max välja ilusaks numbriks, et saaks võtta õige korrektsiooni tabeli
def roundFieldForCorrection(indices):
    #rounds the magnetic field to a nice number (100041.221 -> 100000)
    values = []
    for indices in indices:
        MvsH = ORIGINAL_DATAFRAME["Magnetic Field (Oe)"].loc[indices]
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

#otsib vastava max välja järgi õiged korrektsiooni tabelid
def searchCorrectionTable(folder_path, number):
    #searches for the right table based on the value in the title and returns the filepath
    
    for filename in os.listdir(folder_path):
        # Split the filename on a specific character, e.g., underscore
        parts = filename.split('_') # Modify this based on your file naming convention
        # Check if the number is an exact match in any part of the filename
        
        if str(number) in parts:
            file_path = os.path.join(folder_path, filename)
            return file_path
        
    return None

correction_folder_path = os.path.join(USER_PATH,'PdCorrection tables')

def CorrectionTableToDict(numbers_to_search):
    #returns the corresponding amount of error tables for each unique value in a dictionary form with the key being the value the table is for
    error_tables = {}
    
    for nr in numbers_to_search:
        error_table_path = searchCorrectionTable(correction_folder_path, nr)
        
        if error_table_path is None:
            print(f"{nr} Oe jaoks ei leia parandustabelit ")
            continue
        error_table_measurements = pd.read_csv(error_table_path, index_col=0, encoding="ANSI")
        error_tables[nr] = error_table_measurements
        
    return error_tables

def interpolateTrueField(magnetic_field_values, error_table_measurements):
    x = error_table_measurements["Magnetic Field (Oe)"]
    y = error_table_measurements["True Field (Oe)"]

    # Create an interpolation function
    interp_func = interp1d(x, y, kind='linear', fill_value='extrapolate')

    # Call the interpolation function with magnetic_field_values to get the interpolated values
    true_field_interpolated = interp_func(magnetic_field_values)
    return true_field_interpolated

def interpolateMvsH(separated_MvsH, error_tables):
    #replaces the old field values with the correct interpolated value
    global true_field_interpolated1
    
    interpolated_dict = {}
    for val_pair in separated_MvsH:
        max_val = max(val_pair[0]["Magnetic Field (Oe)"]) # [0] ei pruugi alati õige max anda
        
        for key in error_tables:
            if max_val - 100 <= key <= max_val + 100:
                #print(f"{key} kukub {max_val} vahemikku")
                
                interpolated = []
                for val in val_pair:
                    
                    magnetic_field_values = val["Magnetic Field (Oe)"]
                    true_field_interpolated1 = interpolateTrueField(magnetic_field_values, error_tables[key])

                    true_field_pd = pd.DataFrame({"True Field (Oe)": true_field_interpolated1})
                    interpolated.append(true_field_pd)
                    
                if key not in interpolated_dict:
                    interpolated_dict[key] = []
                    
                interpolated_dict[key].append(interpolated) #algul oli üks tab edasi
            #else:
                #print(f"{key} ei kuku {max_val} vahemikku")
                
    return interpolated_dict

def plotMvsH(raamat, const_T_values, interpolated_MvsH, MvsH_indices):
    #Plots the MvsH measurement pair with different colors
    
    for key1 in raamat:
        
        for df in raamat[key1]:
            
            fig, ax = plt.subplots()
            H1 = df[0]["Magnetic Field (Oe)"]
            M1 = df[0]["Moment (emu)"]
            H2 = df[1]["Magnetic Field (Oe)"]
            M2 = df[1]["Moment (emu)"] 
            
            colorIdx = df[0].iloc[0].name
            Color = ORIGINAL_DATAFRAME["color"].loc[colorIdx]
            
            for key2 in interpolated_MvsH:
                
                for df_interpolated in interpolated_MvsH[key2]:
                    
                    H3 = df_interpolated[0]["True Field (Oe)"]
                    H4 = df_interpolated[1]["True Field (Oe)"]
                    
                    if len(H3) != len(M1) or len(H4) != len(M2):
                        #Siin ongi see, et ei ole interpoleeritud punkte, siis x/y erineva pikkusega ja ei saa joonistada
                        #print(f"Warning: Length mismatch between True Field and Moment data for {key1} K.") #error selles,et pole correction tabelit iga välja tugevuse jaoks
                        continue  # Continue to the next iteration
                        
                    ax.plot(H3, M1, color = Color, label = "True Field Descending", alpha = 0.5)
                    ax.plot(H4, M2, color = Color, label = "True Field Ascending")
                    ax.legend()
            
            
            ax.plot(H1 ,M1, color = "grey", label = "Descending") 
            ax.plot(H2, M2, color = "grey", label = "Ascending")
            
            val = int(key1)
            
            plot_title = f"M vs H at {val} K"
            ax.set_title(plot_title)
            ax.set_xlabel("Magnetic field (Oe)")
            ax.set_ylabel("Moment (emu)")
            ax.legend() #Hetkel legend nimetab selle järgi et esimene tsükkel on kasvav ja teine kahanev ehk eeldus et mõõtmisel temp algas kasvamisest
            ax.grid(True)
            #fig.savefig(f"C:/Users/kevin/OneDrive/Desktop/Andmetöötlus/Projekt_andmed1/MvsH_graph_at_{val}K.png",bbox_inches = "tight", dpi = 200) #laptop
            fig.savefig(os.path.join(save_to_path,f'MvsH_graph_at_{val}K.png'),bbox_inches = "tight", dpi = 200) #PC
            plt.show()
        
    return None
#-----------------------------MvsT specific functions--------------------------------

#const väärtuste põhjal mõõtmise punktid võtab
def getMeasurementMvsT(const_H_values):
    #Saves all the indices of the points that are equal to the predetermined H value
    row_indices = []
    table = ORIGINAL_DATAFRAME['Magnetic Field (Oe)']
    
    for value in const_H_values:

        indices = table.index[table == value].tolist()
        row_indices.append(indices)

    return row_indices

def plotMvsT(raamat, const_H_values, MvsT_indices):
    #Plots the MvsT measurement pair with different colors

    # ascending_color = ""
    # descending_color = ""
    
    for key in raamat:
        
        for df in raamat[key]:

            fig, ax = plt.subplots()
            T1 = df[0]["Temperature (K)"]
            M1 = df[0]["Moment (emu)"]
            T2 = df[1]["Temperature (K)"] if len(df) > 1 else None
            M2 = df[1]["Moment (emu)"] if len(df) > 1 else None
            
            colorIdx = df[0].iloc[0].name
            Color = ORIGINAL_DATAFRAME["color"].loc[colorIdx]
            
            ax.plot(T1,M1,color = Color, label = "Ascending", alpha = 0.5) # peaks tegelt kontrollima kas kasvab või kahaneb
            ax.plot(T2,M2,color = Color, label = "Descending") if len(df) > 1 else None #, marker = "o") #descending ei pea paika kui on alt üle > alt üles mõõtmine
            ax.set_title(f"M vs T at {key} Oe")
            ax.set_xlabel("Temperature (K)")
            ax.set_ylabel("Moment (emu)")
            ax.legend() #Hetkel legend nimetab selle järgi et esimene tsükkel on kasvav ja teine kahanev ehk eeldus et mõõtmisel temp algas kasvamisest
            ax.grid(True)
            fig.savefig(os.path.join(save_to_path,f'MvsT_graph_at_{key}K.png'),bbox_inches = "tight", dpi = 200)
            plt.show()
        
    return None

#---------------------------------------Universal functions that both paths use----------------------------------------------

#eraldab kõik järjestikku kasvavad indeksid ja siis väljastab kõige pikema vahemiku,
#kõige pikem vahemik on mõõtmine ise
def filterMeasurementIndices(unfiltered_indices):
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

#otsib mingi series'e lokaalsed ekstreemumid hilisemaks eraldamiseks
def separationIndexForSingleSeries(data, column_name, n = 100): # https://stackoverflow.com/questions/48023982/pandas-finding-local-max-and-min
    """
    Find local peaks indices (maxima and minima) in a DataFrame or Series.

    Parameters:
    - data: DataFrame or Series.
    - column_name: Name of the column to analyze.
    - n: Number of points to be checked before and after.

    Returns:
    - DataFrame with 'min' and 'max' columns indicating local minima and maxima.
    """
    
    if isinstance(data, pd.Series):
        # Convert a Series to a DataFrame with a specified column name
        data = pd.DataFrame({column_name: data})
    
    index = data.index
    
    # Find local peaks
    relative_min_indices = argrelextrema(data[column_name].values, np.less_equal, order=n)[0]
    relative_max_indices = argrelextrema(data[column_name].values, np.greater_equal, order=n)[0]
    
    min_indices = index[relative_min_indices]
    max_indices = index[relative_max_indices]
    
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
    plt.scatter(data.index, local_peaks['min'], c='r', label='Minima')
    plt.scatter(data.index, local_peaks['max'], c='g', label='Maxima')
    plt.plot(data.index, data[column_name], label=column_name)
    plt.legend()
    plt.show()
        
    return min_indices, max_indices

#Võtab originaalsest tabelist vastavale mõõtmisele vastavad punktid ja edastab need separationIndexForSingleSeries
#funktsioonile eraldamiseks
def separationIndexForMultipleSeries(indices, column_name):
    series = []
    
    for index in indices:
        
        filtered = ORIGINAL_DATAFRAME[column_name].loc[index]
        separation_index = separationIndexForSingleSeries(filtered, column_name)
        series.append(separation_index)
        
    return series

def separateMeasurementWithColorIdx(separation_index, measurement_indices, column):
    #Separates the points based on the index and returns the separated series in pairs with their color
    # global separated_pair1
    separated_pair1 = []
    pair_indices = []
    
    if not isinstance(separation_index, list):
        separation_index = [separation_index]
    
    for min_max_index in separation_index:
        
        if min_max_index[0][0] < min_max_index[1][0]: #kontroll selleks et kas andmed on + - + või - + -
            min_index_list = min_max_index[0].tolist()
            max_index_list = min_max_index[1].tolist()
        else:
            min_index_list = min_max_index[1].tolist()
            max_index_list = min_max_index[0].tolist()
    
        j = 0
        
        for indices in measurement_indices:
    
            tabel = ORIGINAL_DATAFRAME[[column,"Moment (emu)"]].loc[indices] #VB SIIA VÄRV PANNA KUI EI OLE VAHET
        
            for max_index in max_index_list:
                
                separated = []
                index1 = []
                index2 = []
                if max_index in indices:
                    
                    sliced1 = tabel.loc[min_index_list[j]:max_index] #paaride data
                    separated.append(sliced1)
                    
                    ORIGINAL_DATAFRAME.loc[min_index_list[j]:max_index, "color"] = COLORS[0] #värvid paaridele
                    
                    index1 = ORIGINAL_DATAFRAME.loc[min_index_list[j]:max_index].index.tolist() #paaride indeksid
                    
                    if j == len(min_index_list) - 1:
                        sliced2 = tabel.loc[max_index+1:min_index_list[j]]
                        
                        ORIGINAL_DATAFRAME.loc[max_index+1:min_index_list[j], "color"] = COLORS[0]
                        
                        index2 = ORIGINAL_DATAFRAME.loc[max_index+1:min_index_list[j]].index.tolist()
                        
                    else:
                        sliced2 = tabel.loc[max_index+1:min_index_list[j+1]]
                        
                        ORIGINAL_DATAFRAME.loc[max_index+1:min_index_list[j+1], "color"] = COLORS[0]
                        
                        index2 = ORIGINAL_DATAFRAME.loc[max_index+1:min_index_list[j+1]].index.tolist()
                else:
                    max_index_list = max_index_list[j:]

                    break
                
                pair_indices.append(index1 + index2)
                separated.append(sliced2)
                separated_pair1.append(separated)
                
                COLORS.pop(0)
                j += 1
        
    return separated_pair1, pair_indices

def separateIntoDictValuePair(separated_pairs, const_val, column):
    # Creates a dict{const value the measurement was made at: measurement} 
    raamat = {}
    
    for const in const_val:
        
        for val in separated_pairs:
            
            index_to_check = val[0].index[0]
            val_to_check = ORIGINAL_DATAFRAME[column].iloc[index_to_check]
            
            if val_to_check == const or round(val_to_check) == round(const):
                
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
    ORIGINAL_DATAFRAME.plot(x=time_axis, y="Moment (emu)", kind="scatter", c=ORIGINAL_DATAFRAME["color"].values, ax=axes[1])
    ORIGINAL_DATAFRAME.plot(x=time_axis, y="Magnetic Field (Oe)", kind="scatter", c=ORIGINAL_DATAFRAME["color"].values, ax=axes[2])
    
    # Customize axes labels and other properties
    
    axes[0].set_ylabel("Temperature (K)")
    axes[1].set_ylabel("Moment (emu)")
    axes[2].set_ylabel("Magnetic Field (Oe)")
    axes[-1].set_xlabel("Time Stamp (sec)")
    
    # Show the plot
    plt.show()
    
    return None


#-------------- Actually Run the program here -------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

# Ask user input to get a path to datafile
DATAFILE_PATH = askNewDatafilePath()
save_to_path = os.path.dirname(DATAFILE_PATH)

# Read the datafile
HEADER, ORIGINAL_DATAFRAME = readDatafile(DATAFILE_PATH)

# VSM or ACMSII? or maybe HC or ETO in the future
OPTION_TYPE = determineDatafileType(HEADER)
#print(OPTION_TYPE) kas teha nii et funk ise ei prindi ja prindid kui tahad või funktsioon prindib anyway?

#Selle lisasin juurde kuna moment tulbas võib olla nan values ja enne pead kõik õiged tulbad võtma, et need eraldada, muidu eemaldab kõik read,
# sest igas reas on mingi tulp nan value'ga
if OPTION_TYPE == "VSM":
    ORIGINAL_DATAFRAME = ORIGINAL_DATAFRAME[["Time Stamp (sec)", "Temperature (K)", "Magnetic Field (Oe)", "Moment (emu)", "M. Std. Err. (emu)"]].dropna()
    ORIGINAL_DATAFRAME = ORIGINAL_DATAFRAME.reset_index(drop = True)
    
elif OPTION_TYPE == "ACMS":
    ORIGINAL_DATAFRAME = ORIGINAL_DATAFRAME[["Time Stamp (sec)", "Temperature (K)", "Magnetic Field (Oe)", "DC Moment (emu)", "M. Std. Err. (emu)"]].dropna()
    ORIGINAL_DATAFRAME = ORIGINAL_DATAFRAME.reset_index(drop = True)

#parse HEADER
SAMPLE_MASS_g = getMassInGrams(HEADER)
SAMPLE_AREA_CM2 = getAreaCM2(HEADER)
THICKNESS = getThickness(HEADER)

#Color list for color idx and initializing the starting color idx with just black
COLORS = ["red", "green", "blue", "yellow", "brown", "purple", "orange"]
ORIGINAL_DATAFRAME["color"] = "black"

print("_________chechMeasurementType2-----------")  #mis peaks olema selle funktsiooni ebaõnnestumise/veateade? return None? või lihtsalt kaks tühja Seriet?
# Tagastab kaks Pandas.Seriest: temperatures_of_interest, magnetic_fields_of_interest
# edasi peaks veel kontrollima välja filtreerima üksikud punktid mis ei kuulu MvsH ja MvsT andmete hulka
TYPE_TOKEN, TEMPERATURES_OF_INTEREST, MAGNETIC_FIELDS_OF_INTEREST= checkMeasurementType2(ORIGINAL_DATAFRAME)
#print(TEMPERATURES_OF_INTEREST, MAGNETIC_FIELDS_OF_INTEREST)
print("_________end-----------")


print('--------<<<<<<<<<>>>>>>>>>>-----------')
print('--------<<<<<<<<<>>>>>>>>>>-----------')

if MAGNETIC_FIELDS_OF_INTEREST.size <= 0:
    print('no MvsT detected')
    
else:
    print(' MvsT data detected')
    print(MAGNETIC_FIELDS_OF_INTEREST)
    
    unfiltered_MvsT_indices = getMeasurementMvsT(MAGNETIC_FIELDS_OF_INTEREST)
    MvsT_INDICES = filterMeasurementIndices(unfiltered_MvsT_indices)
    
    separation_index_MvsT = separationIndexForMultipleSeries(MvsT_INDICES, "Temperature (K)")# the indices where the separation is going to be done
    
    try:#siin võib juhtuda et liiga väike n siis tulevad valesti ekstreemumid, siis custom error
        separated_MvsT, MvsT_pair_indices = separateMeasurementWithColorIdx(separation_index_MvsT, MvsT_INDICES, "Temperature (K)")
    except IndexError as ie:
        raise ValueError("separationIndexForSingleSeries funktsiooni n argumenti peab muutma, ekstreemumid tulevad valesti sellise n puhul") from ie
        
    DICT_MvsT = separateIntoDictValuePair(separated_MvsT, MAGNETIC_FIELDS_OF_INTEREST, "Magnetic Field (Oe)")
    
    plotMvsT(DICT_MvsT, MAGNETIC_FIELDS_OF_INTEREST, MvsT_INDICES)
print('--------<<<<<<<<<>>>>>>>>>>-----------')
print('--------<<<<<<<<<>>>>>>>>>>-----------')

if TEMPERATURES_OF_INTEREST.size <= 0:
    print('no MvsH detected')
    
else:
    print(' MvsH data detected')
    print(TEMPERATURES_OF_INTEREST)
    
    unfiltered_MvsH_INDICES = getMeasurementMvsH(TEMPERATURES_OF_INTEREST)
    MvsH_INDICES = filterMeasurementIndices(unfiltered_MvsH_INDICES)
    
    separation_indices_MvsH = separationIndexForMultipleSeries(MvsH_INDICES, "Magnetic Field (Oe)")
    
    try:#siin võib juhtuda et liiga väike n siis tulevad valesti ekstreemumid, siis custom error
        SEPARATED_MvsH, MvsH_pair_indices = separateMeasurementWithColorIdx(separation_indices_MvsH, MvsH_INDICES, "Magnetic Field (Oe)")
    except IndexError as ie:
        raise ValueError("separationIndexForSingleSeries funktsiooni n argumenti peab muutma, ekstreemumid tulevad valesti sellise n puhul") from ie
    
    DICT_MvsH = separateIntoDictValuePair(SEPARATED_MvsH, TEMPERATURES_OF_INTEREST, "Temperature (K)")
    
    correction_field_value = roundFieldForCorrection(MvsH_INDICES)
    CORRECTION_TABLES = CorrectionTableToDict(correction_field_value)
    INTERPOLATED_MvsH = interpolateMvsH(SEPARATED_MvsH, CORRECTION_TABLES)
    
    plotMvsH(DICT_MvsH, TEMPERATURES_OF_INTEREST, INTERPOLATED_MvsH, MvsH_INDICES)
    
print('--------<<<<<<<<<>>>>>>>>>>-----------')
print('--------<<<<<<<<<>>>>>>>>>>-----------')

if MAGNETIC_FIELDS_OF_INTEREST.size <= 0 and TEMPERATURES_OF_INTEREST.size <= 0:
    print('Error, ei suutnud eraldada MvsH ja MvsT mõõtmisi')

#Plots temp, field and moment against time
plotMeasurementTimeseries()



#Creating a pandas dataframe for different parameters
# parameters = {"Temperature (K)": None, "Magnetic Field (Oe)": None, "True Field (Oe)": None, "Moment (emu)": None, "M. Std. Err. (emu)": None,\
#         "Moment (emu)/mass(g)": None, "Moment (emu)/area(cm2)": None, "Moment (emu)/volume(cm3)": None, "Susceptibility (emu/g Oe)": None, "1/Susceptibility": None}

# columns = ["Temperature (K)", "Magnetic Field (Oe)", "Moment (emu)", "M. Std. Err. (emu)"]
    
# MvsT_parameters = ORIGINAL_DATAFRAME.loc[MvsT_pair_indices[0], columns]

# MvsT_parameters["Test column"] = range(491)

# parameter_template = pd.DataFrame(data = parameters, index = [MvsT_pair_indices[0]])

# MvsT_test = parameter_template
# MvsH_test = parameter_template

# values = ORIGINAL_DATAFRAME["Temperature (K)"].loc[MvsT_pair_indices[0]]
# MvsT_test["Temperature (K)"] = values.copy






































