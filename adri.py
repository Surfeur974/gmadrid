import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import inspect, re
from IPython.display import display, HTML

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
def print_design_summary(*args):
    BG, FG, COL_KEY, COL_NUM = "#062c50", "#cf7e48", "120px", "85px"

    def fmt(val):
        if not isinstance(val, (int, float, np.floating)):
            return str(val), ""
        a = abs(val)
        for threshold, divisor, suffix in [
            (1e9, 1e9, "G"), (1e6, 1e6, "M"), (1, 1, ""), (1e-3, 1e-3, "m"),
            (1e-6, 1e-6, "u"), (1e-9, 1e-9, "n"), (1e-12, 1e-12, "p"), (1e-15, 1e-15, "f"),
        ]:
            if a >= threshold:
                return f"{val/divisor:.4g}", suffix
        return f"{val:.4e}", ""

    def make_table(title, blk):
        SEP, total = f"border-top:1px solid {FG}", int(COL_KEY[:-2]) + int(COL_NUM[:-2])
        header = (f"<tr><td colspan='2' style='width:{total}px; padding:8px 0 3px; font-size:25px; "
                  f"letter-spacing:.1em; opacity:.5; text-transform:uppercase; font-weight:bold; {SEP}'>{title}</td></tr>")
        footer = f"<tr><td colspan='2' style='{SEP}; padding-top:2px'></td></tr>"
        rows = "".join(
            f"<tr>"
            f"<td style='width:{COL_KEY}; padding:1px 4px 1px 0; opacity:.6; white-space:nowrap; overflow:hidden'>{k}</td>"
            f"<td style='width:{COL_NUM}; padding:1px 0; text-align:right; font-variant-numeric:tabular-nums; font-weight:500; white-space:nowrap'>"
            f"{num}{'&thinsp;' + unit if unit else '&thinsp;&nbsp;'}</td>"
            f"</tr>"
            for k, v in blk.items()
            for num, unit in [fmt(np.asarray(v).flat[0] if not np.isscalar(v) else v)]
        )
        return (f"<table style='font-family:monospace; font-size:18px; border-collapse:collapse; "
                f"color:{FG}; background:{BG}; table-layout:fixed; width:{total}px'>{header}{rows}{footer}</table>")

    frame = inspect.currentframe().f_back
    src = inspect.getframeinfo(frame).code_context[0].strip()
    m = re.search(r'print_design_summary\((.*)\)', src)
    names = [a.strip() for a in m.group(1).split(',')] if m else [f"Block {i}" for i in range(len(args))]

    chunks = [list(zip(names, args))[i:i+4] for i in range(0, len(args), 4)]
    tr_html = "".join(
        f"<tr>{''.join(f'<td style=vertical-align:top;padding-right:16px>{make_table(n,b)}</td>' for n,b in chunk)}</tr>"
        for chunk in chunks
    )
    display(HTML(f"<table style='border-collapse:collapse; background:{BG}'>{tr_html}</table>"))

def folded_cascode(nmoslk, pmoslk, s, d): ##Version Matriced
    """
    Size a folded-cascode amplifier using the gm/ID methodology.

    Sweeps β and gm/ID simultaneously to find the design points that meet
    the noise and speed specifications, then computes the self-loading
    capacitance at each point.

    nmoslk / pmoslk are gm/ID lookup tables for the NMOS and PMOS devices.

    s must contain:
        fu1       — unity-gain frequency (Hz)
        vod_noise — output-referred noise voltage (V)
        FO        — fan-out
        G         — closed-loop gain
        vod_final — voltage ouput step max
        ed — erreur settling %
        VDD
        VCM
        V_SWING_P2P_DIFF

    d must contain:
        gm_ID1   — gm/ID sweep vector for the input transistor (S/A)
        beta     — feedback factor sweep vector
        L1       — input transistor length (µm)
        Lcas     — cascode length (µm)
        gm_IDcas — cascode gm/ID operating point (S/A)
        gamma    — excess noise coefficient
        cself    — self-loading estimate, set to 0 on first pass (F)

    Returns m1 (ID, gm_ID) and p (cltot, cself) at each matched design point.
    """
    VDD = s['VDD']
    VCM = s['VCM']
    V_SWING_P2P_DIFF = s['V_SWING_P2P_DIFF']
    VOUT_CM_HIGH = VCM+V_SWING_P2P_DIFF/2
    VOUT_CM_LOW = VCM-V_SWING_P2P_DIFF/2
    VDS_MOS = VOUT_CM_LOW/2
    beta_max = 1 / (1 + s['G'])
    kB = 1.3806488e-23
    temp = 300
    # ── Device lookups over the gm/ID sweep vector ────────────────────────────
    gm_gds1  = pmoslk.look_up('GM_GDS',  GM_ID=d['gm_ID1'],   L=d['L1'])
    wt1      = pmoslk.look_up('GM_CGG',  GM_ID=d['gm_ID1'],   L=d['L1'])
    cgd_cgg1 = pmoslk.look_up('CGD_CGG', GM_ID=d['gm_ID1'],   L=d['L1'])
    gm_gds2  = nmoslk.look_up('GM_GDS',  GM_ID=d['gm_IDcas'], L=d['Lcas'])
    cdd_gm3  = nmoslk.look_up('CDD_GM',  GM_ID=d['gm_IDcas'], L=d['Lcas'])
    cdd_gm4  = pmoslk.look_up('CDD_GM',  GM_ID=d['gm_IDcas'], L=d['Lcas'])

    # ── Step 2 : sweep β — broadcast to 2D (N_beta, N_gm_ID1) ────────────────
    beta   = d['beta'][:, None]    # (N, 1)
    gm_ID1 = d['gm_ID1'][None, :] # (1, M)

    # Step 3a. Excess noise factor α (eq. 6.54)
    alpha = 2 * d['gamma'] * (1 + 3 * d['gm_IDcas'] / gm_ID1)

    # Step 3b. C_Ltot from noise spec (eq. 6.53) — also fixes C_S, C_F, C_L
    cltot = alpha / beta * kB * temp / s['vod_noise']**2
    CF    = (cltot - d['cself']) / (s['FO'] * s['G'] + 1 - beta)
    CS    = s['G'] * CF
    CL   = s['FO'] * CS

    # Step 3c. Current division factor κ (eq. 6.43)
    kappa = 1 / (1 + gm_ID1 / (gm_gds1 * d['gm_IDcas']) + 2 / gm_gds2)

    #Slewing compute X, then more accurate wu1 
    X = s['vod_final'] * beta / 2 * gm_ID1
    wu1 = 1 / s['ts']*(X - 1 + np.log(1/(s['ed']*X)))

    # Step 3d. g_m1 from settling formula above 
    #gm1 = 2 * np.pi * s['fu1'] * cltot / beta / kappa
    gm1 = wu1 * cltot / beta / kappa

    # Step 3e. I_D1 and f_Ti from the gm/ID vector
    ID1 = gm1 / gm_ID1

    # Step 3f. C_gg1 from g_m and f_Ti → compute C_gd1 and C_in (eq. 6.47)
    cgs1 = gm1 / wt1
    cin  = cgs1 * (1 + cgd_cgg1 * gm_ID1 / d['gm_IDcas'])

    # Step 3g. Actual β values along the gm/ID vector (eq. 6.46)
    beta_actual = CF / (CF + CS + cin)

    # ── Step 4 : find closest match between β_actual and β_k ─────────────────
    idx = np.argmin(np.abs(beta_actual - beta), axis=1)  # (N,)
    j   = np.arange(len(d['beta']))
    


    m1_gm1 = gm1[j, idx]
    m1_gm_id1 = d['gm_ID1'][idx]
    ID1 = gm1[j, idx] / d['gm_ID1'][idx]
    
    #Fill output with MOS charac
    m1 = get_all_mos(pmoslk, GM_ID=m1_gm_id1,L=d['L1'], GM=m1_gm1, ID=ID1)
    m2 = get_all_mos(nmoslk, GM_ID=d['gm_IDcas'], L=d['Lcas'], ID=2*ID1, VDS=VDS_MOS)
    m3 = get_all_mos(nmoslk, GM_ID=d['gm_IDcas'], L=d['Lcas'], ID=ID1, VDS=VDS_MOS)
    m4 = get_all_mos(pmoslk, GM_ID=d['gm_IDcas'], L=d['Lcas'], ID=ID1, VDS=VDS_MOS)
    m5 = get_all_mos(pmoslk, GM_ID=d['gm_IDcas'], L=d['Lcas'], ID=ID1, VDS=VDS_MOS)

    #######Calc GAIN
    gds5    = pmoslk.look_up( 'GDS_ID', GM_ID=d['gm_IDcas'], L=d['Lcas'], VDS=VDS_MOS) * ID1
    gds4    = pmoslk.look_up( 'GDS_ID', GM_ID=d['gm_IDcas'], L=d['Lcas'], VDS=VCM-VDS_MOS, VSB = VDS_MOS) * ID1
    gm4     = pmoslk.look_up( 'GM_ID', GM_ID=d['gm_IDcas'], L=d['Lcas'], VDS=VCM-VDS_MOS, VSB = VDS_MOS) * ID1
    gds3    = nmoslk.look_up( 'GDS_ID', GM_ID=d['gm_IDcas'], L=d['Lcas'], VDS=VCM-VDS_MOS, VSB = VDS_MOS) * ID1
    gm3     = nmoslk.look_up( 'GM_ID', GM_ID=d['gm_IDcas'], L=d['Lcas'], VDS=VCM-VDS_MOS) * ID1
    gds2    = nmoslk.look_up( 'GDS_ID', GM_ID=d['gm_IDcas'], L=d['Lcas'], VDS=VDS_MOS) * 2 * ID1
    gds1    = pmoslk.look_up( 'GDS_ID', GM_ID=m1_gm_id1, L=d['L1']) * 2 * ID1
    gds_gm1    = pmoslk.look_up( 'GDS_GM', GM_ID=m1_gm_id1, L=d['L1'])
    gm_gds2    = pmoslk.look_up( 'GM_GDS', GM_ID=d['gm_IDcas'], L=d['Lcas'], VDS=VDS_MOS)
    gm1 = m1_gm_id1 * ID1

    rout = 1 / ( gds4 / (1+gm4/gds5) + gds3 / (1 + gm3 / (gds1+gds2) ))
    kappa = 1 / ( 1 + gds1 / gm3 + 2 / gm_gds2)
    L0 = 20*np.log10(beta * kappa * gm1 * rout)

    p = {
        'cltot': cltot[j, idx],
        'cself': m1['ID'] * d['gm_IDcas'] * (cdd_gm3 + cdd_gm4),
        #'beta' : beta,
        'L0' : L0,
        'CF' : CF,
        'CS' : CS,
        'CL' : CL,
        'fu1' : wu1 * 1/(2*np.pi),
        'beta_betamax' : beta_actual / d['beta'],
        'FP2 / FU1' : d['fp2'] / (wu1 * 1/(2*np.pi))
    }
    return m1,m2,m3,m4,m5,p

def two_stage_sr_opti(nmos, pmos, s, d):
    kBT = 1.3806488e-23 * 300
    beta_arr     = np.atleast_1d(d['beta'])
    cltot_cc_arr = np.atleast_1d(d['cltot_cc'])
    N = len(beta_arr)   # axe 0
    M = len(cltot_cc_arr)  # axe 1
    m1 = {k: np.full((N, M), np.nan) for k in ['W','gm','gmid','cgg','cgd','cdd','ft','id']}
    m2 = {k: np.full((N, M), np.nan) for k in ['W','gm','gmid','cgg','cgd','cdd','ft','id','cgs','fts']}
    m3 = {k: np.full((N, M), np.nan) for k in ['W','gmid','cdd']}
    m4 = {k: np.full((N, M), np.nan) for k in ['W','gmid','cdd']}
    p  = {k: np.full((N, M), np.nan) for k in ['cc','cltot','cf','cs','cl','c1','rz','cn','rself1','rself2','cc_add']}
    m1['L'] = d['L1']
    m2['L'] = d['L2']
    m3['L'] = d['L3']
    m4['L'] = d['L4']

    for j, beta in enumerate(beta_arr):  # j = indice beta → ligne [j,:]
        cc    = (2/beta * kBT * d['gam1'] * (1 + d['gam3']/d['gam1'] * d['gm3_gm1'])
               + 1/cltot_cc_arr * kBT * (1 + d['gam2'] * (1 + d['gam4']/d['gam2'] * d['gm4_gm2']))
               ) / s['vod_noise']**2
        cltot = cc * cltot_cc_arr
        cf    = cltot / (1 + d['rself_2']) / (1 - beta + s['FO'] * s['G'])
        cs    = cf * s['G']
        cl    = cs * s['FO']
        cgs2  = cc * d['cgs2_cc']
        c1    = cgs2 * (1 + d['rself_1'])
        gmR    = np.sqrt(s['L0'] / beta)
        
        gm_id_test = d['gm_id_test'] 
        N_iter_gm_id_test = len(d['gm_id_test'])

        X = np.zeros((N_iter_gm_id_test,M))
    
        wu1 = np.zeros((N_iter_gm_id_test,M))
        gm1 = np.zeros((N_iter_gm_id_test,M))
        cgg1 = np.zeros((N_iter_gm_id_test,M))
        ft1 = np.zeros((N_iter_gm_id_test,M))
        #gmid1 = np.zeros((N_iter_gm_id_test,M))
        id1 = np.zeros((N_iter_gm_id_test,M))
        cgs1 = np.zeros((N_iter_gm_id_test,M))
        cin = np.zeros((N_iter_gm_id_test,M))
        beta_actual = np.zeros((N_iter_gm_id_test,M))

        for k in range(len(gm_id_test)):
            
            X[k,:] = s['vod_final'] * beta / 2 * gm_id_test[k]
            X[k, X[k,:] < 1] = 1
            #wu1[k,:] = 1 / s['ts']*(X[k,:] - 1 + np.log(X[k,:] * 1/s['ed']))
            wu1[k,:] = 1 / s['ts']*(X[k,:] - 1 + np.log(1/(s['ed']*X[k,:])))
            gm1[k,:]    = (wu1[k,:] * cc / beta * (1 + (1 + c1/cc) / gmR + (1 + cltot/cc) / gmR))
            cgg1[k,:]  =  pmos.lookup('CGG_GM', GM_ID=gm_id_test[k], L=d['L1'], WARNING=False) * gm1[k,:]
            cgs1[k,:]  =  pmos.lookup('CGS_GM', GM_ID=gm_id_test[k], L=d['L1'], WARNING=False) * gm1[k,:]
            #cgg1[k,:]   = cf * (1/beta - 1 - s['G'])
            ft1[k,:]    = gm1[k,:] / (2*np.pi * cgg1[k,:])
            #gmid1[k,:]  = pmos.lookup('GM_ID', GM_CGG=2*np.pi*ft1[k,:], L=d['L1'], WARNING=False)
            id1[k,:]    = gm1[k,:] / gm_id_test[k]
            ##Compute input cap for real beta
            #cgs1[k,:] = cgg1[k,:]#gm1[k,:] / (gm1[k] / cgg1[k])
            cin[k,:] = cgs1[k,:]
            beta_actual[k,:] = cf / (cf + cs + cin[k,:])

        #Take closest point to beta, thorw away others
        index_min_gmid = np.argmin(np.abs(beta_actual - beta), axis=0)

        cols = np.arange(M)

        X       = X[index_min_gmid, cols]
        wu1     = wu1[index_min_gmid, cols]
        gm1     = gm1[index_min_gmid, cols]
        cgg1    = cgg1[index_min_gmid, cols]
        ft1     = ft1[index_min_gmid, cols]
        gmid1   = gm_id_test[index_min_gmid]
        id1     = id1[index_min_gmid, cols]
        cgs1    = cgs1[index_min_gmid, cols]
        cin     = cin[index_min_gmid, cols]

        slew_pct = (X - 1) / (X - 1 + np.log(X / s['ed']))
        #print(slew_pct)
        #######Compute actaul beta, find closest to beta, store only this value



        gm2    = 2*np.pi * s['fp2'] * c1 * (1 + cltot/cc + cltot/c1)
        fts2   = gm2 / (2*np.pi * cgs2)
        gmid2  = nmos.lookup('GM_ID', GM_CGS=2*np.pi*fts2, L=d['L2'], WARNING=False)
        id2    = gm2 / gmid2

        # Storage [j,:] → ligne j = beta j, colonnes = tous les cltot_cc
        p['cc'][j,:]     = cc
        p['cltot'][j,:]  = cltot
        p['cf'][j,:]     = cf
        p['cs'][j,:]     = cs
        p['cl'][j,:]     = cl
        p['c1'][j,:]     = c1
        m2['cgs'][j,:]   = cgs2
        m1['gm'][j,:]    = gm1
        m1['cgg'][j,:]   = cgg1
        m1['ft'][j,:]    = ft1
        m1['gmid'][j,:]  = gmid1
        m1['id'][j,:]    = id1
        m2['gm'][j,:]    = gm2
        m2['fts'][j,:]   = fts2
        m2['gmid'][j,:]  = gmid2
        m2['id'][j,:]    = id2
        p['rz'][j,:]     = 1 / gm2

        m3['gmid'][j,:]  = gmid1 * d['gm3_gm1']
        m4['gmid'][j,:]  = gmid2 * d['gm4_gm2']
        m1['W'][j,:]     = id1 / pmos.lookup('ID_W', GM_ID=gmid1, L=d['L1'], WARNING=False)
        m2['W'][j,:]     = id2 / nmos.lookup('ID_W', GM_ID=gmid2, L=d['L2'], WARNING=False)
        m3['W'][j,:]     = id1 / nmos.lookup('ID_W', GM_ID=m3['gmid'][j,:], L=d['L3'], WARNING=False)
        m4['W'][j,:]     = id2 / pmos.lookup('ID_W', GM_ID=m4['gmid'][j,:], L=d['L4'], WARNING=False)

        m1['cgd'][j,:]   = m1['W'][j,:] * pmos.lookup('CGD_W', GM_ID=gmid1, L=d['L1'], WARNING=False)
        p['cn'][j,:]     = m1['cgd'][j,:]
        m2['cgd'][j,:]   = m2['W'][j,:] * nmos.lookup('CGD_W', GM_ID=gmid2, L=d['L2'], WARNING=False)
        p['cc_add'][j,:] = cc - m2['cgd'][j,:]

        m1['cdd'][j,:]   = m1['W'][j,:] * pmos.lookup('CDD_W', GM_ID=gmid1, L=d['L1'], WARNING=False)
        m2['cdd'][j,:]   = m2['W'][j,:] * nmos.lookup('CDD_W', GM_ID=gmid2, L=d['L2'], WARNING=False)
        m3['cdd'][j,:]   = m3['W'][j,:] * nmos.lookup('CDD_W', GM_ID=m3['gmid'][j,:], L=d['L3'], WARNING=False)
        m4['cdd'][j,:]   = m4['W'][j,:] * pmos.lookup('CDD_W', GM_ID=m4['gmid'][j,:], L=d['L4'], WARNING=False)

        p['rself1'][j,:] = (m1['cdd'][j,:] + m3['cdd'][j,:]) / m2['cgs'][j,:]
        p['rself2'][j,:] = ((m2['cdd'][j,:] - m2['cgd'][j,:]) + m4['cdd'][j,:]) / (cl + (1 - beta) * cf)
        #p['rself2'][j,:] = 0

        gds1 = m1['W'][j,:] * pmos.lookup('GDS_W', GM_ID = gmid1, L = d['L1'], WARNING=False)
        gds3 = m3['W'][j,:] * nmos.lookup('GDS_W', GM_ID = m3['gmid'][j,:], L = d['L3'], WARNING=False)
        gds2 = m2['W'][j,:] * nmos.lookup('GDS_W', GM_ID = gmid2, L = d['L2'], WARNING=False)
        gds4 = m4['W'][j,:] * pmos.lookup('GDS_W', GM_ID = m4['gmid'][j,:], L = d['L4'], WARNING=False)

        R1 = 1/(gds1 + gds3)
        R2 = 1/(gds2 + gds4)
        L0 = beta * gm1 * R1 * gm2 * R2
        L0 = 20*np.log10(L0)
        p['L0'] = L0

    return m1, m2, m3, m4, p