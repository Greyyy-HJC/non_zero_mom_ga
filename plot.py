# %%
import matplotlib.pyplot as plt
import numpy as np
from funcs import *

grey = "#808080" 
red = "#FF6F6F" 
peach = "#FF9E6F" 
orange = "#FFBC6F" 
sunkist = "#FFDF6F"
yellow = "#FFEE6F"
lime = "#CBF169"
green = "#5CD25C" 
turquoise = "#4AAB89"
blue = "#508EAD" 
grape = "#635BB1"
violet = "#7C5AB8" 
fuschia = "#C3559F"

color_ls = ['orange','dodgerblue','blueviolet','deeppink','indigo','rosybrown','greenyellow','cyan','fuchsia','royalblue', 'red','green','orange','dodgerblue','blueviolet','deeppink','indigo','rosybrown','greenyellow','cyan','fuchsia','royalblue', 'red','green']

def fit_on_data_R(data_set_tidy, mom, current, title, ylim=None):
    mo = '_' + str(mom)
    hash_key = 'p_sq_{}_pz_0'.format(mom)
    pt2_0_ls = data_set_tidy['p_sq_0_pz_0']['2pt']
    pt2_mom_ls = data_set_tidy[hash_key]['2pt']





    #!# fill between
    # t_tsep_tau = {}
    # t_tsep_tau['pt2_0'] = np.arange(33)
    # t_tsep_tau['pt2_{}'.format(mom)] = np.arange(33)

    # tsep_ls = []
    # tau_ls = []
    # for tsep in range(3, 10): #* keep same 
    #     for tau in range(1, tsep): #* keep same 
    #         tsep_ls.append(tsep)
    #         tau_ls.append(tau)

    # t_tsep_tau['pt3_A3_0'] = [np.array(tsep_ls), np.array(tau_ls)]
    # t_tsep_tau['pt3_V4_0'] = [np.array(tsep_ls), np.array(tau_ls)]
    # t_tsep_tau['pt3_A3_{}'.format(mom)] = [np.array(tsep_ls), np.array(tau_ls)]
    # t_tsep_tau['pt3_V4_{}'.format(mom)] = [np.array(tsep_ls), np.array(tau_ls)]

    # fit_pt2_0_ls = fcn(t_tsep_tau, post)['pt2_0']
    # fit_pt2_mom_ls = fcn(t_tsep_tau, post)['pt2'+mo]
    # fit_pt3_collect = list(fcn(t_tsep_tau, post)['pt3_'+current+mo])


    # fit_x_dic = {}
    # fit_y_dic = {}
    # for tsep in range(3, 10): #* keep same
    #     tau_ls = []
    #     fit_pt3_ls = []

    #     for tau in range(1, tsep): #* keep same 
    #         tau_ls.append(tau)
    #         fit_pt3_ls.append(fit_pt3_collect.pop())

    #     fit_x_dic['tsep_{}'.format(tsep)] = np.array(tau_ls) - tsep/2
    #     fit_y_dic['tsep_{}'.format(tsep)] = pt2_pt3_to_R(tsep, tau_ls, fit_pt2_0_ls, fit_pt2_mom_ls, np.array(fit_pt3_ls))







    fig = plt.figure(figsize=fig_size)
    ax = plt.axes(plt_axes)
    for tsep in range(3, 10):
        pt3_ls = data_set_tidy[hash_key][current+'_tsep_{}'.format(tsep)][1:-1]
        tau_ls = np.arange(tsep+1)[1:-1]
        R_tsep = pt2_pt3_to_R(tsep, tau_ls, pt2_0_ls, pt2_mom_ls, pt3_ls)

        ax.errorbar(tau_ls - tsep/2, [v.mean for v in R_tsep], [v.sdev for v in R_tsep], color=color_ls[tsep], label='tsep {}'.format(tsep), **errorb)




        #!# fill between
        # y_gv = fit_y_dic['tsep_{}'.format(tsep)]

        # print(fit_x_dic['tsep_{}'.format(tsep)])
        # print(y_gv)

        # ax.fill_between(fit_x_dic['tsep_{}'.format(tsep)], [v.mean + v.sdev for v in y_gv], [v.mean - v.sdev for v in y_gv], alpha=0.4, color=color_ls[tsep])



    ax.tick_params(direction='in', **ls_p)
    ax.grid(linestyle=':')
    ax.set_ylim(ylim)
    plt.title(title, fs_p)
    plt.legend()
    plt.show()

    return
