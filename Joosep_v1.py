# -*- coding: utf-8 -*-
"""
Created on Thu Aug 31 11:41:28 2023

@author: joose
"""

import os
import statistics as stat
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
    data = data.reset_index(drop=True)
    
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

# MEASUREMENT_TABLE["Temperature (K)"].round(decimals=1).value_counts().median()
# meanTemp=MEASUREMENT_TABLE["Temperature (K)"].round(decimals=1).value_counts().mean()
# MEASUREMENT_TABLE["Temperature (K)"].round(decimals=1).value_counts()>meanTemp
# MEASUREMENT_TABLE["Temperature (K)"].round(decimals=1).value_counts().mode()
# MEASUREMENT_TABLE["Temperature (K)"].round(decimals=1).value_counts()

def checkMeasurementType(measurement_table, discrete_detection_ration = 0.07):
    #Checks what measurements the file contains
    
    #global uniques_T,uniques_H,codes_T
    
    type_token = {"Temperature": None, "Field": None}
    
    dataset_T = measurement_table["Temperature (K)"]
    dataset_H = measurement_table["Magnetic Field (Oe)"]
      
    smoothed_dataset_T = dataset_T.round(decimals=1)
    smoothed_dataset_H = dataset_H.round(decimals=1)
    
    
    codes_T, uniques_T = pd.factorize(smoothed_dataset_T)
    codes_H, uniques_H = pd.factorize(smoothed_dataset_H)
    
    ratio_T = uniques_T.size/smoothed_dataset_T.size
    ratio_H = uniques_H.size/smoothed_dataset_H.size
    print(f"T : {ratio_T}, H : {ratio_H }")

    
    if ratio_T < discrete_detection_ration: #discrete
    
        if ratio_H < discrete_detection_ration: #discrete
        
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
    
        if ratio_H < discrete_detection_ration: #discrete
        
            type_token["Temperature"] = "continous"
            type_token["Field"] = "discrete"
            print("T continous, H discrete = MvsT \n")
            
            return type_token
        
        else: #continous
        
            type_token["Temperature"] = "continous"
            type_token["Field"] = "continous"
            print("T continous, H continous = both \n")
            
            return type_token


#-------------- Actually Run the program here -------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------------------

# Ask user input to get a path to datafile
DATAFILE_PATH = askNewDatafilePath()

# Read the datafile
HEADER, ORIGINAL_DATAFRAME = readDatafile(DATAFILE_PATH)

# VSM or ACMSII? or maybe HC or ETO in the future
OPTION_TYPE = determineDatafileType(HEADER)
print(OPTION_TYPE)


#parse HEADER
SAMPLE_MASS_g = getMassInGrams(HEADER)
SAMPLE_AREA_CM2 = getAreaCM2(HEADER)
THIKNESS = getThickness(HEADER)

type_token = checkMeasurementType(ORIGINAL_DATAFRAME)

if type_token["Temperature"] == "discrete" and type_token["Field"] == "continous":
    print("\n Kas jõudis MvsH algusesse?")
    
    print(ORIGINAL_DATAFRAME["Temperature (K)"].round(decimals=1).value_counts())   #annab mõistlikud numbrid
    
    # ORIGINAL_DATAFRAME["Temperature (K)"].round(decimals=1).value_counts().median()
    # meanTemp=ORIGINAL_DATAFRAME["Temperature (K)"].round(decimals=1).value_counts().mean()
    # ORIGINAL_DATAFRAME["Temperature (K)"].round(decimals=1).value_counts()>meanTemp
    # ORIGINAL_DATAFRAME["Temperature (K)"].round(decimals=1).value_counts().mode()
    # ORIGINAL_DATAFRAME["Temperature (K)"].round(decimals=1).value_counts()
    
    
    
elif type_token["Temperature"] == "continous" and type_token["Field"] == "discrete":
    print(ORIGINAL_DATAFRAME["Magnetic Field (Oe)"].round(decimals=1).value_counts())   #annab mõistlikud numbrid
    
elif type_token["Temperature"] == "continous" and type_token["Field"] == "continous": #!!!  
    print("segafail")
  
meanTempCount=ORIGINAL_DATAFRAME["Temperature (K)"].round(decimals=1).value_counts().mean()                     #see osa siin annab segafailide puhul mõistlikud numbrid. Aga mõne erilisema faili korral ei saa ta aru et sega fail on, aga need read annavad ka siis mõistlikud numbrid. Äkki saame kuidagi hoopis neid prinditavaid liste kasutada faili ära tundmiseks.
muutuja = ORIGINAL_DATAFRAME["Temperature (K)"].round(decimals=1).value_counts() > meanTempCount*10
print(ORIGINAL_DATAFRAME["Temperature (K)"].round(decimals=1).value_counts()[muutuja])
type(ORIGINAL_DATAFRAME["Temperature (K)"].round(decimals=1).value_counts()[muutuja].index)
    #ORIGINAL_DATAFRAME["Temperature (K)"].round(decimals=1).value_counts()[muutuja].index.toSeries()
    
meanTempCount=ORIGINAL_DATAFRAME["Magnetic Field (Oe)"].round(decimals=1).value_counts().mean()
muutuja = ORIGINAL_DATAFRAME["Magnetic Field (Oe)"].round(decimals=1).value_counts() > meanTempCount*10
print(ORIGINAL_DATAFRAME["Magnetic Field (Oe)"].round(decimals=1).value_counts()[muutuja])

#ORIGINAL_DATAFRAME["Temperature (K)"].round(decimals=1).value_counts()[muutuja].index[ORIGINAL_DATAFRAME["Temperature (K)"].round(decimals=1).value_counts()[muutuja].index.size-1]
