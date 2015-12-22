# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 14:08:19 2015

@author: michaelcrawley
"""

###Import modules
import os
import re
import numpy as np
from glob import glob
from scipy.interpolate import InterpolatedUnivariateSpline
from datetime import datetime

###Class and function definitions
class model_params():
    def __init__(self,P_inf,Temp,Velocity,area,chord,span,loadcell,Xmrp=0,Ymrp=0,Zmrp=0,AoAB=0,AoAM=0,AoSB=0,AoSM=0,AoRB=0,AoRM=0):
        self.area = area
        self.chord = chord
        self.span = span
        self.velocity = Velocity
        self.loadcell = loadcell
        self.P_inf = P_inf
        self.Temp = Temp
        self.Xmrp = Xmrp
        self.Ymrp = Ymrp
        self.Zmrp = Zmrp
        self.AoAB = AoAB
        self.AoAM = AoAM
        self.AoSB = AoSB
        self.AoSM = AoSM
        self.AoRB = AoRB
        self.AoRM = AoRM
        
def get_dim(nparray,dim):
    size = nparray.shape
    return size[dim-1]
    
def get_AoX(fname):
    tmp = re.search('_AoA[0-9.-]+',fname)
    tmp = tmp.group(0)
    AoA = float(tmp[4:])
    tmp = re.search('_AoS[0-9.-]+',fname)
    tmp = tmp.group(0)
    AoS = float(tmp[4:])
    return AoA, AoS
    
def Calibration(loadcell):
    loadcell = loadcell.lower()
    if loadcell == 'jr3':
        cal = np.loadtxt('JR3_CalibrationMatrix.txt')
        forces = np.array((-1,1,-1))
        moments = np.array((-1,1,1))
        f_indx = lambda x: x[:,(2,1,0)]*np.tile(forces,(get_dim(x,1),1))
        m_indx = lambda x: x[:,(5,4,3)]*np.tile(moments,(get_dim(x,1),1))
    elif loadcell == 'jr3_fz_up':
        cal = np.loadtxt('JR3_CalibrationMatrix.txt')
        forces = np.array((1,-1,1))
        moments = np.array((1,-1,1))
        f_indx = lambda x: x[:,(0,1,2)]*np.tile(forces,(get_dim(x,1),1))
        m_indx = lambda x: x[:,(3,4,5)]*np.tile(moments,(get_dim(x,1),1))
    elif loadcell == 'ati_n25_fz_up':
        cal = np.loadtxt('ATI_N25_FT14574.txt')
        f_indx = lambda x: x[:,(0,1,2)]
        m_indx = lambda x: x[:,(3,4,5)]
    elif loadcell == 'ati_n25_fz_down':
        cal = np.loadtxt('ATI_N25_FT14574.txt')
        forces = np.array((1,-1,-1))
        moments = np.array((1,-1,-1))
        f_indx = lambda x: x[:,(0,1,2)]*np.tile(forces,(get_dim(x,1),1))
        m_indx = lambda x: x[:,(3,4,5)]*np.tile(moments,(get_dim(x,1),1))
    else:
        raise NameError('Incorrect Load Cell Definition')
    return cal, f_indx, m_indx
    
def process_tare(fname,cal):
    voltage = np.loadtxt(fname,skiprows=1)
    loads = np.dot(voltage,cal)
    return np.mean(loads,axis=0)

def process_data(fname,cal,tare_means,tare_AoA):
    voltage = np.loadtxt(fname,skiprows=1)
    loads = np.dot(voltage,cal)
    raw_mean = np.mean(loads,axis=0)
    raw_std = np.std(loads,axis=0)
    
    int_tare_fun = InterpolatedUnivariateSpline(tare_AoA,tare_means)
    
def aero_processing(src,params):
    start=datetime.now()
    #Constant parameters
    rho = params.P_inf/(params.Temp*287.05)
    Q = 0.5 * rho * params.velocity**2 #dynamic pressure
    AoSB = params.AoSB*np.pi/180 #convert to radians
    AoSM = (params.AoSM + params.AoSM)*np.pi/180
    RollB = params.AoRB*np.pi/180
    RollM = (params.AoRB + params.AoRM)*np.pi/180
    delAoS = AoSM-AoSB
    delRoll = RollM-RollB
    
    #Define coordinate transformation matrices
    Rot1 = lambda angle: np.array((np.cos(angle),np.sin(angle),0),(-np.sin(angle),np.cos(angle),0),(0,0,1))
    Rot2 = lambda angle: np.array((np.cos(angle),0,np.sin(angle)),(0,1,0),(-np.sin(angle),0,np.cos(angle)))
    Rot3 = lambda angle: np.array((1,0,0),(0,np.cos(angle),np.sin(angle)),(0,-np.sin(angle),np.cos(angle)))
    
    #Get calibration data & functions
    cal, f_indx, m_indx = Calibration(params.loadcell)
    
    #Find all tare and data files
    tare_files = glob(os.path.join(src,'Tare','*Raw.wtd'))
    data_files = glob(os.path.join(src,'*Raw.wtd'))
    
    #Sort and Process Tare data   
    num_tare = len(tare_files)      #number of tare files
    #num_cpu_tare = min((num_tare,mp.cpu_count()))         #number of cores to use when processing tare files
    tare_AoA = np.zeros((num_tare))
    for k in range(num_tare):
        tare_AoA[k] = get_AoX(tare_files[k])[0]
        
    tare_indx = np.argsort(tare_AoA,axis=0)
    tare_AoA = tare_AoA[tare_indx]
    tare_files = [tare_files[k] for k in tare_indx]     #because python is total bullshit when it comes to lists...
    tare_means = np.zeros((num_tare,6)) #output of the load cells is always of DIM 6
    for k in range(len(tare_files)):
        tare_means[k,:] = process_tare(tare_files[k],cal)
        
    print(datetime.now()-start)
    
    return tare_means