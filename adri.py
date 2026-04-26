import numpy as np
from scipy.interpolate import interp1d
# Interpolation function
def interp1(x, y, value):
    # Perform cubic interpolation using interp1d
    f = interp1d(x, y, kind='cubic')
    return f(value)

def currentMirror(moslk, spec, param):
        #Arg :  
        #spec. : input_current
        #param. : l_mos, gm_id1
    #return 
        #perf : Rout, vin_min, vout_min
        #MOS. : Name, W, L, VGS, VDS, gm_id
    #Init output struc
    param_length = len(param)
    perf = {'rout' : np.full(param_length, np.nan), 'vin_min' : np.full(param_length, np.nan), 'vout_max' : np.full(param_length, np.nan)}
    mos_unit = {'NAME' : np.full(param_length, np.nan),
           'W' : np.full(param_length, np.nan),
           'L' : np.full(param_length, np.nan),
           'VGS' : np.full(param_length, np.nan),
           'VDS' : np.full(param_length, np.nan),
           'GM_ID' : np.full(param_length, np.nan)}
    mos = {'MOS1' : mos_unit,'MOS2' : mos_unit}


    #Get parameter from args
    gm_id1 = param['gm_id1']
    l_mos = param['l_mos']
    input_current = spec['input_current']

    #Two Step evaluation of VGS, with 1er eval to get a 1st VDS value
    vgs1 = moslk.lookupVGS(GM_ID=gm_id1, L=l_mos)
    vds1 = vgs1
    vgs1 = moslk.lookupVGS(GM_ID=gm_id1, VDS=vds1, L=l_mos)
    jd1  = moslk.lookup('ID_W', VGS=vgs1, VDS=vds1, L=l_mos)
    w1 = input_current / jd1
    rout = 1 / (w1 * moslk.look_up('GDS_W',VGS=vgs1, VDS=vds1, L=l_mos))
    vin_min = vgs1
    vout_min = vds1


    perf = {'Rout' : rout, 'vin_min' : vin_min, 'vout_min' : vout_min}
    mos['MOS1'] = {'NAME' : "MOS1", 'W' : w1, 'L' : l_mos, 'VGS' : vgs1, 'VDS' : vds1, 'GM_ID' : gm_id1}
    mos['MOS2'] = {'NAME' : "MOS2", 'W' : w1, 'L' : l_mos, 'VGS' : vgs1, 'VDS' : vds1, 'GM_ID' : gm_id1}

    return mos, perf

def currentMirrorCascode(moslk, spec, param):
        #Arg :  
        #spec. : input_current
        #param. : l_mos, gm_id1
    #return 
        #perf : Rout, vin_min, vout_min
        #MOS. : Name, W, L, VGS, VDS, gm_id
    #Init output struc
    param_length = len(param)
    perf = {'rout' : np.full(param_length, np.nan), 'vin_min' : np.full(param_length, np.nan), 'vout_max' : np.full(param_length, np.nan)}
    mos_unit = {'NAME' : np.full(param_length, np.nan),
           'W' : np.full(param_length, np.nan),
           'L' : np.full(param_length, np.nan),
           'VGS' : np.full(param_length, np.nan),
           'VDS' : np.full(param_length, np.nan),
           'GM_ID' : np.full(param_length, np.nan)}
    mos = {'MOS1' : mos_unit,'MOS2' : mos_unit,'MOS3' : mos_unit,'MOS4' : mos_unit,'MOS6' : mos_unit,'MOS7' : mos_unit}
    s = np.linspace(0,400) * 0.001

    #Get parameter from args
    gm_id1 = param['gm_id1']
    l_mos = param['l_mos']
    input_current = spec['input_current']

    vds1 = 2.0 / gm_id1 + 0.05

    #Start compute MOS
    vgs1 = moslk.lookupVGS(GM_ID=gm_id1, VDS=vds1, L=l_mos)
    jd1  = moslk.lookup('ID_W', VGS=vgs1, VDS=vds1, L=l_mos)
    jd2  = moslk.lookup('ID_W', VGS=vgs1+s, VDS=vgs1-vds1, L=l_mos, VSB=-vds1)
    vgs2 = interp1(jd2/jd1, vgs1+s, 1)
    vbias = vds1+vgs2
    w1 = input_current / jd1

    vea1 = moslk.look_up('ID_GDS', VGS=vgs1, VDS=vds1, L=l_mos,VSB=0)
    gm4_id = moslk.look_up('GM_ID', VGS=vgs2, VDS=vgs1-vds1, VSB=-vds1, L=l_mos)
    gmb4_id = moslk.look_up('GMB_ID', VGS=vgs2, VDS=vgs1-vds1, VSB=-vds1, L=l_mos)
    gds4_id = moslk.look_up('GDS_ID', VGS=vgs2, VDS=vgs1-vds1, VSB=-vds1, L=l_mos)

    A4 = (gm4_id+gmb4_id+gds4_id+1/vea1) / gds4_id

    VEA_cascode = vea1 * A4
    rout = VEA_cascode / input_current

    #Bias circuit
    jd6 = moslk.look_up('ID_W', VGS=vgs2, VDS=vgs2, VSB=-vds1, L=l_mos)
    gm6_id = moslk.look_up('GM_ID', VGS=vgs2, VDS=vgs2, VSB=-vds1, L=l_mos)
    id6 = w1 * jd6
    jd7 = moslk.look_up('ID_W', VGS=vbias, VDS = vds1, L=l_mos)

    w7 = input_current / jd7

    vin_min = vds1 + vgs2 #VT+2Vov
    vout_min = vds1 + 2.0/gm4_id #2Vov
    #vin_min = 2.0 * 2.0/gm_id1 + 0.6
    #vout_min = 2.0 * 2.0/gm_id1

    gm2_id = moslk.look_up('GM_ID', VGS=vgs2, VDS=vgs1-vds1, VSB=-vds1, L=l_mos)
    gm3_id = moslk.look_up('GM_ID', VGS=vgs2, VDS=vbias-vgs2, VSB=0, L=l_mos)
    gm7_id = moslk.look_up('GM_ID', VGS=vbias, VDS = vds1, L=l_mos)



    perf = {'Rout' : rout, 'vin_min' : vin_min, 'vout_min' : vout_min}
    mos['MOS1'] = {'NAME' : "MOS1", 'W' : w1, 'L' : l_mos, 'VGS' : vgs1, 'VDS' : vds1, 'GM_ID' : gm_id1}
    mos['MOS2'] = {'NAME' : "MOS2", 'W' : w1, 'L' : l_mos, 'VGS' : vgs2, 'VDS' : vgs1-vds1, 'GM_ID' : gm2_id}
    mos['MOS3'] = {'NAME' : "MOS3", 'W' : w1, 'L' : l_mos, 'VGS' : vgs1, 'VDS' : vbias-vgs2, 'GM_ID' : gm3_id}
    mos['MOS4'] = {'NAME' : "MOS4", 'W' : w1, 'L' : l_mos, 'VGS' : vgs2, 'VDS' : vgs1-vds1, 'GM_ID' : gm4_id}
    mos['MOS6'] = {'NAME' : "MOS6", 'W' : w1, 'L' : l_mos, 'VGS' : vgs2, 'VDS' : vgs2, 'GM_ID' : gm6_id}
    mos['MOS7'] = {'NAME' : "MOS7", 'W' : w7, 'L' : l_mos, 'VGS' : vbias, 'VDS' : vds1, 'GM_ID' : gm7_id}

    return mos, perf