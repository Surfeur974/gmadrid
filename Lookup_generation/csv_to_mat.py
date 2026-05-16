import pandas as pd
import numpy as np
#from numpy.core.records import fromarrays
from numpy.rec import fromarrays  # Updated to use the public API
from scipy.io import savemat


def convert_csv_to_mat(device, w, nfing, txt_path, infomos):
    ##HOW TO USE ##
    ##Get csv convert to .MAT -> pray it work -> done
    # device = 'nmos'
    # w = 1.0
    # nfing = 1
    # txt_path = '../lookup_table/nmos_xh018.txt'
    # infomos = 'XH018 NMOS NEI 180nm
    # pathout = convert_txt_to_mat(device, w, nfing, txt_path, infomos)

    df_raw = pd.read_csv(txt_path)
    df = df_raw
    df = df_raw.drop(df.columns[:7], axis=1)
    df.columns = df.columns.str.replace('.1', '')
    df.columns = df.columns.str.replace('sfl', 'n1overf')
    df = df.apply(pd.to_numeric)
    df=abs(df)

    #print(df)
    # sweep variable vectors
    l =   np.unique(df['l'])*1e6
    vgs = np.unique(df['vg'])
    vds = np.unique(df['vd'])
    vsb = np.unique(df['vb'])

    # --- détection simple ordre sweep ---
    def get_sweep_order(df):
        cols = ['l','vb','vd','vg']
        runs = {}
        for c in cols:
            v = df[c].values
            run = 1
            max_run = 1
            for i in range(1, len(v)):
                if v[i] == v[i-1]:
                    run += 1
                    if run > max_run:
                        max_run = run
                else:
                    run = 1
            runs[c] = max_run
        return sorted(runs, key=runs.get)  # fastest → slowest

    order = get_sweep_order(df)

    # tailles
    sizes = {
        'l': len(l),
        'vb': len(vsb),
        'vd': len(vds),
        'vg': len(vgs)
    }
    # reshape dans l'ordre réel (slowest → fastest)
    shape = tuple(sizes[c] for c in reversed(order))

    def reshape_param(name):
        arr = np.reshape(df[name].values, shape)

        # ordre actuel
        current = list(reversed(order))

        # ordre pygmid
        target = ['l','vg','vd','vb']

        perm = [current.index(t) for t in target]

        return np.transpose(arr, perm)
    print('l','vb','vd','vg')
    print('Shape is : ', shape)

    id  = reshape_param('id')
    vt  = reshape_param('vth')
    gm  = reshape_param('gm')
    gmb = reshape_param('gmbs')
    gds = reshape_param('gds')

    cgg = reshape_param('cgg') + reshape_param('cgdo') + reshape_param('cgso')
    cgb = reshape_param('cgb')
    cgd = reshape_param('cgd') + reshape_param('cgdo')
    cgs = reshape_param('cgs') + reshape_param('cgso')
    cdd = reshape_param('cdd') + reshape_param('capbd') + reshape_param('cgdo')
    css = reshape_param('css') + reshape_param('capbs') + reshape_param('cgso')

    ####A verifier que c'est **2 ou pas
    sth = reshape_param('nid')
    sfl = reshape_param('n1overf')

    dic = {
    "INFO": infomos,
    "CORNER": "NOM",
    "TEMP": 300.0,
    "VGS": vgs,
    "VDS": vds,
    "VSB": vsb,
    "L": l,
    "W": w,
    "NFING": nfing,
    "ID": id,
    "VT": vt,
    "GM": gm,
    "GMB": gmb,
    "GDS": gds,
    "CGG": cgg,
    "CGB": cgb,
    "CGD": cgd,
    "CGS": cgs,
    "CDD": cdd,
    "CSS": css,
    "STH": sth,
    "SFL": sfl
    }
    pathout = '../lookup_table/'+device+'.mat'
    savemat(pathout, {device: dic})
    #savemat('../lookup_table/nmos_xh018.mat', {'nmos_xh018': dic})
    input_table = df
    print('TXT to MAT DONE')
    return pathout, input_table