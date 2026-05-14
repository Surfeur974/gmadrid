import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import inspect, re

def coolPlot(xvar, yvar, xname, yname, title):
    fig, ax1 = plt.subplots(figsize=(6, 4))
    ax1.plot(xvar, yvar, color='black', linewidth=1.5)
    ax1.set_xlabel(xname, fontsize=12)
    ax1.set_ylabel(yname, fontsize=12)
    ax1.grid(True, which='both', linestyle='-', linewidth=0.5, color='gray', alpha=0.7)
    #ax1.ticklabel_format(axis='y', style='sci', scilimits=(-4, -4), useMathText=True)
    ax1.legend(loc='center right', bbox_to_anchor=(1, 0.8), labels='label')
    ax1.set_title(title, loc='left', fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.show()

def interp1(x, y, value):
    # Interpolation function
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

def ota_5T_L_gmid(nmoslk,pmoslk, spec, param):
        #Arg :  
        #spec. : FU, CL, VDD, V_IC
        #param. : l_mos, gm_id1, gm_id2
    #return 
        #perf : AV0, ibias, Vs
        #MOS. : Name, W, L, VGS, VDS, gm_id
    #Init output struc
    param_length = len(param)
    perf = {'av0' : np.full(param_length, np.nan), 'ibias' : np.full(param_length, np.nan)}
    mos_unit = {'NAME' : np.full(param_length, np.nan),
           'W' : np.full(param_length, np.nan),
           'L' : np.full(param_length, np.nan),
           'VGS' : np.full(param_length, np.nan),
           'VDS' : np.full(param_length, np.nan),
           'GM_ID' : np.full(param_length, np.nan)}
    mos = {'MOS1' : mos_unit,'MOS2' : mos_unit}

    #Get spec from args
    fu = spec['fu']
    cload = spec['cload']
    vdd = spec['vdd']
    vic = spec['vic']
    #Get parameter from args
    l_mos = param['l_mos']
    gm_id1 = param['gm_id1']
    gm_id2 = param['gm_id2']
    ####################################
    ####################################
    ####################################
    #gm and ID, from Spec and Hypothesis
    gm1 = 2*3.14159*fu*cload
    ibias = gm1 / gm_id1
    #VGS from lookupVGS, L and GM_ID
    # Find VDS1 and 
    vgs1 = nmoslk.lookupVGS(GM_ID=gm_id1, L=l_mos, VDS=0.6, VSB=0) #Fist estimation with vds default
    vgs2 = pmoslk.lookupVGS(GM_ID=gm_id2, L=l_mos, VDS=0.6, VSB=0) #Fist estimation with vds default
    vgs2 = pmoslk.lookupVGS(GM_ID=gm_id2, L=l_mos, VDS=vgs2, VSB=0) #Second with VDS = VGS (current mirror)
    vd = vdd - vgs2
    vs = vic - vgs1
    vds1 = vd - vs
    vgs1 = nmoslk.lookupVGS(GM_ID=gm_id2, L=l_mos, VDS=vds1, VSB=0) #Second with VDS more accurate
    vds2 = vgs2

    #Calcul Gain
    gds_id1 = nmoslk.look_up('GDS_ID', GM_ID=gm_id1, VDS=vds1, L=l_mos, VSB=0)
    gds_id2 = pmoslk.look_up('GDS_ID', GM_ID=gm_id2, VDS=vds2, L=l_mos, VSB=0)
    av0 = gm_id1 / (gds_id1 + gds_id2)

    #print('AV0 = %.2F' % AV0)
    #print('gm1 = %.2F uS' % (gm1 * 1e6))
    #print('ID = %.2F uA before self loading' % (ID*1e6))

    cselfloading = 0
    for m in range(1,10,1):
        #self loading
        gm1 = 2*3.14159*fu*(cload+cselfloading)
        ibias = gm1 / gm_id1
        #Denormalization
        jd1 = nmoslk.lookup('ID_W', GM_ID=gm_id1, VDS=vds1, L=l_mos, VSB=0)
        jd2 = pmoslk.lookup('ID_W', GM_ID=gm_id2, VDS=vds2, L=l_mos, VSB=0)
        w1 = ibias/jd1
        w2 = ibias/jd2
        #print('ID = %.2F uA' % (id* 1e6))
        #print('W1 = %.2F um' % (W1))
        #print('W2 = %.2F um' % (W2))
        cdd1 = w1 * nmoslk.look_up('CDD_W', GM_ID=gm_id1, L=l_mos, VSB=0)
        cdd2 = w2 * pmoslk.look_up('CDD_W', GM_ID=gm_id2, L=l_mos, VSB=0)
        cgg1 = w1 * nmoslk.look_up('CGG_W', GM_ID=gm_id1, L=l_mos, VSB=0)
        cselfloading = cdd1 + cdd2 
        #+cgg1 if follower
    #print('ID = %.2F uA before self loading' % (ID*1e6))

    ###################################
    ###################################
    ###################################
    ###################################
    perf = {'AV0' : av0, 'ibias' : ibias, 'vs' : vs}
    mos['MOS1'] = {'NAME' : "MOS1", 'W' : w1, 'L' : l_mos, 'VGS' : vgs1, 'VDS' : vds1, 'GM_ID' : gm_id1}
    mos['MOS2'] = {'NAME' : "MOS2", 'W' : w2, 'L' : l_mos, 'VGS' : vgs2, 'VDS' : vds2, 'GM_ID' : gm_id2}
    #print('AV0 = %.2F dB' % (20*np.log10(AV0)))
    #print('VS = %.2F V ' % vs)
    return mos, perf

def cool_print(*vars):
    frame = inspect.currentframe().f_back
    call_line = inspect.getframeinfo(frame).code_context[0].strip()
    
    match = re.search(r'cool_print\((.+)\)', call_line)
    if match:
        var_names = [v.strip() for v in match.group(1).split(',')]
    else:
        var_names = ["?"] * len(vars)
    
    for name, val in zip(var_names, vars):
        print(f"{name} = {val!r}")
def get_all_mos(moslk, GM_ID, L, GM=None, VDS=None, VSB=None, ID=None):
    mos = {}
    if VSB is None:
        VSB = 0
    if VDS is None:
        VDS = moslk['VDS'][-1]/2
    if ID is None:
        ID = GM / GM_ID
    if GM is None:
        GM = GM_ID * ID
    
    mos['ID'] = ID
    mos['GM_ID'] = GM_ID
    mos['ID_W'] = moslk.look_up('ID_W', GM_ID=GM_ID, L=L, VDS=VDS, VSB=VSB)
    mos['GDS_ID'] = moslk.look_up('GDS_ID', GM_ID=GM_ID, L=L, VDS=VDS, VSB=VSB)
    mos['W'] = ID / mos['ID_W']
    mos['L'] = L
    return mos
    