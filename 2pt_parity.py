# %%
import numpy as np
import h5py as h5
import gvar as gv
import lsqfit as lsf

from funcs import *
from plot import *

#!# If you simultaneously fit the positive and negative parity correlators, does the spectrum significantly change based on if you fit only the positive parity correlation functions?

#!# Model 1 -- only fit positive parity to original correlation functions
    #* Cpp = Sum [ A_n^pp * A_n^pp * Exp[ -E_n^pp * t ] ]
    #* Cnp = Sum [ A_n^np * A_n^np * Exp[ -E_n^np * t ] ]
#!# Model 2 -- fit positive parity and negative parity correlation functions together
    #* Cpp = Sum [ A_n^pp * A_n^pp * Exp[ -E_n^pp * t ] ] + Sum [ B_m^np * B_m^np * Exp[ -E_m^np * t ] ]
    #* Cnp = Sum [ A_n^np * A_n^np * Exp[ -E_n^np * t ] ] + Sum [ B_m^pp * B_m^pp * Exp[ -E_m^pp * t ] ]

#* hadron: 'proton', 'proton_np'





p_sq_cut = 35  #* cut off of p^2
mom_list = [1, 2, 4, 5, 8, 9, 10, 13, 16, 17, 18, 20, 25, 32] # no zero-mom
pt2_n = 3


def find_data_2pt(hadron):
    fname = 'data/spec_4D_a12m130_a_avg_cfgs_300-5295_srcs_16-23_fft_n6.h5' #* full 2pt data

    # fname = 'data/spec_4D_a12m130_a_tslice_avg_cfgs_300-5295_srcs_0-31_fft_n6.h5' #* previous

    basekey = 'gf1p0_w3p0_n30_M51p2_L520_a3p0/spec_4D/ml0p00195'
    basekey_2 = '4D_correlator/spin_avg'

    myfile = h5.File(fname, 'r')[basekey]
    data = myfile[hadron][basekey_2]

    return data

def rotation_avg(data_ls):
    #!# need to separate correlators by their pz, not just the magnitude of the momentum
    #!# just take pz = 0 here, so it is (0,:,:)

    #* for (1000, xx, 13, 13, 13)

    hash_dic = {}

    for px in range(-6, 7):
        for py in range(-6, 7):
            pz = 0
            p_sq = px ** 2 + py ** 2 + pz ** 2
            if p_sq <= p_sq_cut: #* cut off
                hash_key = 'p_sq_{}_pz_{}'.format(p_sq, pz)
                if hash_key not in hash_dic:
                    hash_dic[hash_key] = []
                hash_dic[hash_key].append( data_ls[:,:,pz,py,px] )

    data_avg_ls = []
    hash_key_ls = []
    for key in hash_dic:
        hash_key_ls.append(key)
        data_avg_ls.append( np.average(hash_dic[key], axis=0) ) # shape = (len(hash), 1000, xx)

    data_avg_ls = np.swapaxes(data_avg_ls, 0, 1)
    data_avg_ls = np.swapaxes(data_avg_ls, 1, 2) # shape = (1000, xx, len(hash))
    
    return hash_key_ls, data_avg_ls

def data_check_meff(data_set_tidy, p_sq):
    hash_key = 'p_sq_{}_pz_0'.format(p_sq)


    pt2_ls = data_set_tidy[hash_key]['pp_avg'] #! positive parity
    # pt2_ls = pt2_ls[::-1] #! negative parity
    t_ls = np.arange(len(pt2_ls))
    meff_ls = pt2_to_meff(pt2_ls)

    errorbar_plot(t_ls[:-1][:20], [v.mean for v in meff_ls][:20], [v.sdev for v in meff_ls][:20], 'meff_pp_avg_{}'.format(hash_key))


    pt2_ls = data_set_tidy[hash_key]['np_avg']
    # pt2_ls = -pt2_ls[::-1] #! positive parity
    # pt2_ls = pt2_ls[::-1] #! negative parity
    t_ls = np.arange(len(pt2_ls))
    meff_ls = pt2_to_meff(pt2_ls)

    errorbar_plot(t_ls[:-1][:20], [v.mean for v in meff_ls][:20], [v.sdev for v in meff_ls][:20], 'meff_np_avg_{}'.format(hash_key))


    return 




def set_priors(mom_ls=None):
    #* XXpp{m}_{n} m means energy states, n means momentum

    prior = gv.BufferDict()

    #!# zero mom
    for i in range(pt2_n):
        prior['App{}_0'.format(i)] = gv.gvar(0, 0.01)
        prior['Anp{}_0'.format(i)] = gv.gvar(0, 0.01)
        prior['Bpp{}_0'.format(i)] = gv.gvar(0, 0.01)
        prior['Bnp{}_0'.format(i)] = gv.gvar(0, 0.01)




    prior['Epp0_0'] = gv.gvar(0.6, 0.06)
    prior['log(dEpp1_0)'] = gv.gvar(-1.07, 0.3) # ln(0.94-0.6)
    prior['log(dEpp2_0)'] = gv.gvar(-0.69, 0.5) # ln(1.44-0.94)

    prior['log(dEppmax_0)'] = gv.gvar(-1.25, 0.5*5)


    prior['Enp0_0'] = gv.gvar(0.6, 0.06)
    prior['log(dEnp1_0)'] = gv.gvar(-1.07, 0.3) # ln(0.94-0.6)
    prior['log(dEnp2_0)'] = gv.gvar(-0.69, 0.5) # ln(1.44-0.94)

    prior['log(dEnpmax_0)'] = gv.gvar(-1.25, 0.5*5)




    #!# non zero mom
    if mom_ls == None:
        mom_ls = mom_list

    for mom in mom_ls:
        mo = '_'+str(mom)

        for i in range(pt2_n):
            prior['App{}'.format(i)+mo] = gv.gvar(0, 0.008)
            prior['Anp{}'.format(i)+mo] = gv.gvar(0, 0.008)
            prior['Bpp{}'.format(i)+mo] = gv.gvar(0, 0.008)
            prior['Bnp{}'.format(i)+mo] = gv.gvar(0, 0.008)



        Emean = np.sqrt(0.6 ** 2 + (2*np.pi*np.sqrt(mom)/48) ** 2)
        prior['Epp0'+mo] = gv.gvar(Emean, Emean/5)
        prior['Enp0'+mo] = gv.gvar(Emean, Emean/5)

        dE1mean =  np.sqrt((0.94) ** 2 + (2*np.pi*np.sqrt(mom)/48) ** 2) - Emean
        prior['log(dEpp1'+mo+')'] = gv.gvar(np.log(dE1mean), 0.5)
        prior['log(dEnp1'+mo+')'] = gv.gvar(np.log(dE1mean), 0.5)

        dE2mean =  np.sqrt((1.44) ** 2 + (2*np.pi*np.sqrt(mom)/48) ** 2) - Emean - dE1mean
        prior['log(dEpp2'+mo+')'] = gv.gvar(np.log(dE2mean), 0.8)
        prior['log(dEnp2'+mo+')'] = gv.gvar(np.log(dE2mean), 0.8)

        prior['log(dEppmax'+mo+')'] = gv.gvar(-1.25, 0.5*10)
        prior['log(dEnpmax'+mo+')'] = gv.gvar(-1.25, 0.5*10)

    return prior


def fit_mode_1(pt2_t, p, mom):
    mo = '_'+str(mom) # to indicate the value of momentum

    E_list = {}
    for i in range(pt2_n): # initialize       
        E_list['Epp'+str(i)] = p['Epp0'+mo]
        E_list['Enp'+str(i)] = p['Enp0'+mo]

    for i in range(1, pt2_n): # define Ei  
        for j in range(1, i+1):
            if i == pt2_n-1 and j == i:
                E_list['Epp'+str(i)] += p['dEppmax'+mo]
                E_list['Enp'+str(i)] += p['dEnpmax'+mo]
            else:
                E_list['Epp'+str(i)] += p['dEpp'+str(j)+mo]
                E_list['Enp'+str(i)] += p['dEnp'+str(j)+mo]

    val = {}
    val['pp'] = p['App0'+mo] * p['App0'+mo] * np.exp( - E_list['Epp0'] * pt2_t['pp'] )
    val['np'] = p['Anp0'+mo] * p['Anp0'+mo] * np.exp( - E_list['Enp0'] * pt2_t['np'] )

    for i in range(1, pt2_n):
        val['pp'] += p['App'+str(i)+mo] * p['App'+str(i)+mo] * np.exp( - E_list['Epp'+str(i)] * pt2_t['pp'] )
        val['np'] += p['Anp'+str(i)+mo] * p['Anp'+str(i)+mo] * np.exp( - E_list['Enp'+str(i)] * pt2_t['np'] )

    return val

def fit_mode_2(pt2_t, p, mom):
    mo = '_'+str(mom) # to indicate the value of momentum

    E_list = {}
    for i in range(pt2_n): # initialize       
        E_list['Epp'+str(i)] = p['Epp0'+mo]
        E_list['Enp'+str(i)] = p['Enp0'+mo]

    for i in range(1, pt2_n): # define Ei  
        for j in range(1, i+1):
            if i == pt2_n-1 and j == i:
                E_list['Epp'+str(i)] += p['dEppmax'+mo]
                E_list['Enp'+str(i)] += p['dEnpmax'+mo]
            else:
                E_list['Epp'+str(i)] += p['dEpp'+str(j)+mo]
                E_list['Enp'+str(i)] += p['dEnp'+str(j)+mo]

    return

def fit_func(mom_ls):
    #!# 0 is not in mom_ls
    def fcn(x, p):
        val = {}

        pt2_t = {}
        pt2_t['pp'] = x['pt2_0_pp']
        pt2_t['np'] = x['pt2_0_np']

        val['pt2_0_pp'] = fit_mode_1(pt2_t, p, 0)['pp']
        val['pt2_0_np'] = fit_mode_1(pt2_t, p, 0)['np']

        for mom in mom_ls:
            mo = '_'+str(mom)
            pt2_t = {}
            pt2_t['pp'] = x['pt2'+mo+'_pp']
            pt2_t['np'] = x['pt2'+mo+'_np']

            val['pt2'+mo+'_pp'] = fit_mode_1(pt2_t, p, mom)['pp']
            val['pt2'+mo+'_np'] = fit_mode_1(pt2_t, p, mom)['np']

        return val
    return fcn

def fit(data_dic, t_dic, mom_ls, priors=None):
    if priors == None:
        priors = set_priors(mom_ls)

    x = {}
    amp = {}

    #* mom = 0
    x['pt2_0_pp'] = t_dic['p_sq_0_pz_0']['pp'] 
    x['pt2_0_np'] = t_dic['p_sq_0_pz_0']['np'] 

    ti = x['pt2_0_pp'][0]
    tf = x['pt2_0_pp'][-1] + 1
    amp['pt2_0_pp'] = data_dic['p_sq_0_pz_0']['pp_avg'][ti:tf]

    ti = x['pt2_0_np'][0]
    tf = x['pt2_0_np'][-1] + 1
    amp['pt2_0_np'] = data_dic['p_sq_0_pz_0']['np_avg'][ti:tf]


    #* mom != 0 
    for mom in mom_ls:
        mo = '_'+str(mom)
        x['pt2'+mo+'_pp'] = t_dic[hash_key]['pp']
        x['pt2'+mo+'_np'] = t_dic[hash_key]['np']

        ti = x['pt2'+mo+'_pp'][0]
        tf = x['pt2'+mo+'_pp'][-1] + 1
        amp['pt2'+mo+'_pp'] = data_dic[hash_key]['pp_avg'][ti:tf]

        ti = x['pt2'+mo+'_np'][0]
        tf = x['pt2'+mo+'_np'][-1] + 1
        amp['pt2'+mo+'_np'] = data_dic[hash_key]['np_avg'][ti:tf]

    fit_res = lsf.nonlinear_fit(data=(x, amp), prior=priors, fcn=fit_func(mom_ls), maxit=10000, fitter='scipy_least_squares') # scipy_least_squares   # gsl_multifit

    return fit_res


# %%
p_data = np.array( find_data_2pt('proton') )
p_np_data = np.array( find_data_2pt('proton_np') )

hash_key_ls, p_ro_avg = rotation_avg(p_data)
hash_key_ls, p_np_ro_avg = rotation_avg(p_np_data)


data_set = {}
data_set['pp'] = p_ro_avg
data_set['np'] = p_np_ro_avg

data_set_avg = gv.dataset.avg_data(data_set)

#* make the dict as [hash_key][pp/np][:]
data_set_tidy = {}
for hash_key in hash_key_ls:
    data_set_tidy[hash_key] = gv.BufferDict()


for key in data_set_avg:
    for idx in range(len(hash_key_ls)):
        hash_key = hash_key_ls[idx]
        data_set_tidy[hash_key][key] = np.swapaxes(data_set_avg[key], 0, 1)[idx]


#!# combine _pp and _np to get positive / negative parity
for hash_key in data_set_tidy:
    data_set_tidy[hash_key]['pp_avg'] = ( data_set_tidy[hash_key]['pp'] - data_set_tidy[hash_key]['np'][::-1] ) / 2
    data_set_tidy[hash_key]['np_avg'] = ( data_set_tidy[hash_key]['pp'][::-1] - data_set_tidy[hash_key]['np'] ) / 2

data_check_meff(data_set_tidy, 0)




# %%
#!# fit
t_dic = {}
for hash_key in hash_key_ls:
    t_dic[hash_key] = {}
    t_dic[hash_key]['pp'] = np.arange(3, 10)
    t_dic[hash_key]['np'] = np.arange(3, 9)

mom_ls = []
fit_res = fit(data_set_tidy, t_dic, mom_ls)

print(fit_res.format(100))




# %%
hash_key = 'p_sq_0_pz_0'
pp_np = 'pp'
mom_plot = 0
pt2_ls = data_set_tidy[hash_key][pp_np+'_avg'][:20]

meff_plot(pt2_ls, ti=3, tf=9, fit_res=fit_res, mom_ls=mom_ls, mom_plot=mom_plot, title='meff_{}_fit_on_data_{}'.format(pp_np, hash_key), pp_np=pp_np)





# %%
#!# teeeeest

test = data_set_tidy['p_sq_2_pz_0']['np']

errorbar_plot(np.arange(64)[-15:], [v.mean for v in test][-15:], [v.sdev for v in test][-15:], 'test')

test = data_set_tidy['p_sq_2_pz_0']['pp'][:20]

meff_ls = pt2_to_meff(test)

errorbar_plot(np.arange(64)[:15], [v.mean for v in meff_ls][:15], [v.sdev for v in meff_ls][:15], 'test')


# %%
test = {}
for mom in [0, 1, 2, 4]:
    test = data_set_tidy['p_sq_{}_pz_0'.format(mom)]['pp']

    meff_ls = pt2_to_meff(test)

    errorbar_plot(np.arange(64)[:15], [v.mean for v in test][:15], [v.sdev for v in test][:15], 'test_C2_pp_Q2_{}'.format(mom))
# %%
