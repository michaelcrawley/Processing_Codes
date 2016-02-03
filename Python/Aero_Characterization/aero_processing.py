# -*- coding: utf-8 -*-
"""
@author: michaelcrawley
"""

###Import modules
import os
import re
import time
import numpy as np
import pandas as pd
from glob import glob

###Class and function definitions
class model_params:
    def __init__(self,velocity,area,chord,span,loadcell,Xmrp=0,Ymrp=0,Zmrp=0,AoAB=0,AoAM=0,AoSB=0,AoSM=0,AoRB=0,AoRM=0):
        self.area = area
        self.chord = chord
        self.span = span
        self.velocity = velocity
        self.loadcell = loadcell
        self.P_inf = None  #we'll get this info later
        self.Temp = None   #we'll get this info later
        self.Xmrp = Xmrp
        self.Ymrp = Ymrp
        self.Zmrp = Zmrp
        self.AoAB = AoAB
        self.AoAM = AoAM
        self.AoSB = AoSB * np.pi/180
        self.AoSM = (AoSB + AoSM)*np.pi/180;
        self.AoRB = AoRB
        self.AoRM = AoRM
        self.AoSB = AoSB
        self.RollB = AoRB*np.pi/180
        self.RollM = (AoRB + AoRM)*np.pi/180
        self.delAoS = (AoSM - AoSB)*np.pi/180
        self.delRoll = self.RollB + self.RollM
        self.delAoA = AoAM * np.pi / 180

class data_struct:
    def __init__(self,CL,CD,CS,CRM,CPM,CYM):
        self.CL = np.mean(CL,axis=0)
        self.CD = np.mean(CD,axis=0)
        self.CS = np.mean(CS,axis=0)
        self.CRM = np.mean(CRM,axis=0)
        self.CPM = np.mean(CPM,axis=0)
        self.CYM = np.mean(CYM,axis=0)
        self.CLstd = np.std(CL,axis=0)
        self.CDstd = np.std(CD,axis=0)
        self.CSstd = np.std(CS,axis=0)
        self.CRMstd = np.std(CRM,axis=0)
        self.CPMstd = np.std(CPM,axis=0)
        self.CYMstd = np.std(CYM,axis=0)

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

def interp_extrap(AoA_list,int_AoA,loads):
    for_sorting = abs(AoA_list - int_AoA)
    indx = np.argsort(for_sorting)
    sorted_AoA = AoA_list[indx]
    sorted_loads = loads[(indx,)]
    int_loads = (int_AoA-sorted_AoA[0])/(sorted_AoA[0]-sorted_AoA[1])*(sorted_loads[0,]-sorted_loads[1,]) + sorted_loads[0,]
    return int_loads

def Calibration(loadcell):
    loadcell = loadcell.lower()
    if loadcell == 'jr3_fz_down':
        fname = 'JR3_CalibrationMatrix.txt'
        forces = np.array((-1,1,-1))
        moments = np.array((-1,1,1))
        f_indx = lambda x: x[:,(2,1,0)]*np.tile(forces,(get_dim(x,1),1))
        m_indx = lambda x: x[:,(5,4,3)]*np.tile(moments,(get_dim(x,1),1))
    elif loadcell == 'jr3_fz_up':
        fname = 'JR3_CalibrationMatrix.txt'
        forces = np.array((1,-1,1))
        moments = np.array((1,-1,1))
        f_indx = lambda x: x[:,(0,1,2)]*np.tile(forces,(get_dim(x,1),1))
        m_indx = lambda x: x[:,(3,4,5)]*np.tile(moments,(get_dim(x,1),1))
    elif loadcell == 'ati_n25_fz_up':
        fname = 'ATI_N25_FT14574.txt'
        f_indx = lambda x: x[:,(0,1,2)]
        m_indx = lambda x: x[:,(3,4,5)]
    elif loadcell == 'ati_n25_fz_down':
        fname = 'ATI_N25_FT14574.txt'
        forces = np.array((1,-1,-1))
        moments = np.array((1,-1,-1))
        f_indx = lambda x: x[:,(0,1,2)]*np.tile(forces,(get_dim(x,1),1))
        m_indx = lambda x: x[:,(3,4,5)]*np.tile(moments,(get_dim(x,1),1))
    else:
        raise NameError('Incorrect Load Cell Definition')

    cal = np.loadtxt(fname)
    return cal, f_indx, m_indx

def Save_TXT(src,params,aero):
    fname = os.path.join(src,'Results.txt')
    with open(fname,'w') as fid:
        fid.write('Ref_Area[m2]\tRef_Chord_Length[m]\tRef_Span_b[m]\tXc[m]\tYc[m]\tZc[m]\tBalance_AoA_Offset[deg]\tModel_AoA_Offset[deg]\tBalance_Side_Slip_Offset[deg]\tModel_Side_Slip_Slip_Offset[deg]\tBalance_Roll_Offset[deg]\tModel_Roll_Offset[deg]\n')
        params_str = str(params.Area)
        params_str += '\t' + str(params.Chord)
        params_str += '\t' + str(params.span)
        params_str += '\t' + str(params.Xmrp)
        params_str += '\t' + str(params.Yrmp)
        params_str += '\t' + str(params.Zmrp)
        params_str += '\t' + str(params.AoAB)
        params_str += '\t' + str(params.AoAM)
        params_str += '\t' + str(params.AoSB)
        params_str += '\t' + str(params.AoSM)
        params_str += '\t' + str(params.AoRB)
        params_str += '\t' + str(params.AoRM)
        fid.write(params_str)
        fid.write('AoA\tCL\tCD\tCS\tCRM\tCPM\tCYM\tSD(CL)\tSD(CD)\tSD(CS)\tSD(CRM)\tSD(CPM)\tSD(CYM)')

        for k in range(len(aero)):
            data_str = str(aero(k).AoA)
            data_str += '\t' + str(aero(k).CL)
            data_str += '\t' + str(aero(k).CD)
            data_str += '\t' + str(aero(k).CS)
            data_str += '\t' + str(aero(k).CRM)
            data_str += '\t' + str(aero(k).CPM)
            data_str += '\t' + str(aero(k).CYM)
            data_str += '\t' + str(aero(k).CLstd)
            data_str += '\t' + str(aero(k).CDstd)
            data_str += '\t' + str(aero(k).CSstd)
            data_str += '\t' + str(aero(k).CRMstd)
            data_str += '\t' + str(aero(k).CPMstd)
            data_str += '\t' + str(aero(k).CYMstd)
            fid.write(data_str + '\n')

def aero_processing(src,params):

    #Constant parameters
    rho = params.P_inf/(params.Temp*287.05)
    Q = 0.5 * rho * params.velocity**2 #dynamic pressure

    #Define coordinate transformation matrices
    Rot1 = lambda angle: np.array((np.cos(angle),np.sin(angle),0),(-np.sin(angle),np.cos(angle),0),(0,0,1))
    Rot2 = lambda angle: np.array((np.cos(angle),0,np.sin(angle)),(0,1,0),(-np.sin(angle),0,np.cos(angle)))
    Rot3 = lambda angle: np.array((1,0,0),(0,np.cos(angle),np.sin(angle)),(0,-np.sin(angle),np.cos(angle)))

    #Get calibration data & functions
    cal, f_indx, m_indx = Calibration(params.loadcell)

    #Find all tare and data files
    tare_files = glob(os.path.join(src,'Tare','*Raw.wtd'))
    data_files = glob(os.path.join(src,'*Raw.wtd'))
    version = 8;
    if len(data_files) == 0: # either we have no data, or it is in binary files
        tare_files = glob(os.path.join(src,'*.tare'))
        data_files = glob(os.path.join(src,'*.raw'))
        version = 9

    ####Sort and Process Tare data
    num_tare = len(tare_files)      #number of tare files
    tare_AoA = np.zeros((num_tare))
    for k in range(num_tare):
        tare_AoA[k] = get_AoX(tare_files[k])[0]

    tare_indx = np.argsort(tare_AoA,axis=0)
    tare_AoA = tare_AoA[tare_indx]
    tare_files = [tare_files[k] for k in tare_indx]     #because python is total bullshit when it comes to lists...
    tare_means = np.zeros((num_tare,6)) #output of the load cells is always of DIM 6
    for k in range(num_tare):
        if version == 8:
            voltage = pd.read_csv(tare_files[k],sep='\t',comment='R')
        elif version == 9:
            fid = open(tare_files[k],'rb')
            tmp = np.fromfile(fid,dtype='<f')
            voltage = np.reshape(tmp,(len(tmp)/6,6),'C')
            fid.close()
        loads = np.dot(cal,voltage.T).T
        tare_means[k,:] = np.mean(loads,axis=0)


    ####Sort and Process Aero data
    num_data = len(data_files)
    data_AoA = np.zeros((num_data))
    for k in range(num_data):
        data_AoA[k] = get_AoX(data_files[k])[0]

    data_indx = np.argsort(data_AoA,axis=0)
    data_AoA = data_AoA[data_indx]
    data_files = [data_files[k] for k in data_indx]
    aero = [] #initialize data container
    for k in range(num_data):
        #Calculate rotation matrices
        AoAB = (AoA[k] + params.AoAB)*np.pi/180
        AoAM = (AoA[k] + params.AoAB + params.AoAM)*np.pi/180
        delAoA = AoAM-AoAB;    #Angle difference from y axis
        aircraft_body_rot = np.dot(np.dot(Rot3(delRoll),Rot2(delAoA)),Rot1(delAoS))
        stability_rot = np.dot(Rot1(AoSM),Rot2(AoAM))
        rotation = np.dot(stability_rot,aircraft_body_rot)

        #Load data
        if version == 8:
            voltage = pd.read_csv(data_files[k],sep='\t',comment='R')
        elif version == 9:
            fid = open(data_files[k],'rb')
            tmp = np.fromfile(fid,dtype='<f')
            voltage = np.reshape(tmp,(len(tmp)/6,6),'C')
            fid.close()
        loads = np.dot(cal,voltage.T).T
        raw_mean = np.mean(loads,axis=0)
        raw_std = np.std(loads,axis=0)

        #Correct for wind-off tares
        int_tare = interp_extrap(tare_AoA,AoA,tare_means)
        corrected = loads - np.tile(int_tare,(loads.shape[0],1))
        corrected_mean = np.mean(corrected,axis=0)

        #convert to aero axis
        forces = np.dot(rotation,f_indx(corrected).T)
        moment = np.dot(rotation,m_indx(correct).T)
        MRP = np.dot(rotation,np.array([params.Xmrp,params.Ymrp,params.Zmrp]).T)

        D = np.transpose(forces[0,])
        S = np.transpose(forces[1,])
        L = np.transpose(forces[2,])
        RMT = np.transpose(moments[0,])
        PMT = np.transpose(moments[1,])
        YMT = np.transpose(moments[2,])
        XP = MRP[0]
        YP = MRP[1]
        ZP = MRP[2]

        #Calculate Aerodynamic coefficients
        RM = RMT+L*YP-S*ZP
        PM = PMT-L*XP-D*ZP
        YM = YMT-D*YP-S*XP

        CL = L/Q/params.Area
        CD = D/Q/params.Area
        CS = S/Q/params.Area

        CRM = RM/Q/params.Area/params.span
        CPM = PM/Q/params.Area/params.Chord
        CYM = YM/Q/params.Area/params.span

        #append data container
        aero.append(data_struct(CL,CD,CS,CRM,CPM,CYM))

    #save formatted data
    processing_file = os.path.basename(__file__)
    processing_date = time.strftime("%d/%m/%Y")
    Save_TXT(src,params,aero)
    np.savez(os.path.join(src,'Results.npz'),params,src,aero,tare,processing_date,processing_file)
    return aero
