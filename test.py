# %%
import gvar as gv
import numpy as np
import matplotlib.pyplot as plt

data_set_tidy = gv.load('dump/data_set_tidy')

print(data_set_tidy['p_sq_0_pz_0']['2pt'])

pt2 = data_set_tidy['p_sq_0_pz_0']['2pt']

plt.errorbar(np.arange(len(pt2)), [v.mean for v in pt2], [v.sdev for v in pt2])

# %%
print([key for key in data_set_avg])
print([np.shape(data_set_avg[key]) for key in data_set_avg])

# %%
print(len(hash_key_ls))


# %%
print( data_set_tidy['p_sq_34_pz_0']['A3_tsep_3'] )
print( np.swapaxes(data_set_avg['A3_tsep_3'],0,1)[0] )


# %%
import gvar as gv

n3 = gv.load('dump/A3_00_0_ls_till_p_sq_32_n3')
n4 = gv.load('dump/A3_00_0_ls_till_p_sq_32_n4')

print(n3)
print(n4)
# %%
