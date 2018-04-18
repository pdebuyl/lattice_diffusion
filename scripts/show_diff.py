import sys
import numpy as np
import matplotlib.pyplot as plt

conf = {}
with open(sys.argv[1], 'r') as f:
    for l in f:
        if '=' in l:
            idx = l.find('=')
            k = l[:idx]
            v = l[idx+1:]
            if len(k)>0 and len(v)>0:
                conf[k.strip()] = v.strip()

p_flip = float(conf['p_flip'])
p_move = float(conf['p_move'])
n_inner = int(conf['n_inner'])
n_particles = int(conf['n_particles'])

do_diffusion = conf['do_diffusion'] == 'T'
do_drive = conf['do_drive'] == 'T'

D = 0
if do_diffusion:
    D = D + 0.5*p_move
if do_drive:
    D = D + 0.5 + 1/(np.exp(2*p_flip)-1)

print(do_diffusion, do_drive, D)

data = np.loadtxt('fort.31')

n_bins = data.shape[1]
xr = np.arange(n_bins)

def diffusion(x, x0, D, t):
    norm_factor = 1/np.sqrt(4*np.pi*D*t)
    res = np.exp(-(x-x0)**2/(4*D*t))
    return res/np.sum(res)

for i, step in enumerate(data):
    l, = plt.plot(xr, step)
    plt.plot(xr, n_particles*diffusion(xr, n_bins//2, D, n_inner*(i+1)), ls='--', color=l.get_color())

plt.show()
