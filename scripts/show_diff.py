#!/usr/bin/env python3
import argparse
import numpy as np
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser()
parser.add_argument('config_file')
parser.add_argument('data_file', nargs='+')
parser.add_argument('--show-n', action='store_true')
parser.add_argument('--show-d', action='store_true')
parser.add_argument('--show-rho0', action='store_true')
parser.add_argument('--stride', default=1, type=int)
parser.add_argument('--no-plot', action='store_false')
args = parser.parse_args()

conf = {}
with open(args.config_file, 'r') as f:
    for l in f:
        if '=' in l:
            idx = l.find('=')
            k = l[:idx]
            v = l[idx+1:]
            if len(k)>0 and len(v)>0:
                conf[k.strip()] = v.strip()

p_flip = float(conf['p_flip'])
p_drive = float(conf['p_drive'])
dt = float(conf['dt'])
n_inner = int(conf['n_inner'])
n_particles = int(conf['n_particles'])
rho_0 = int(conf['rho_0'])

do_diffusion = conf['do_diffusion'] == 'T'
do_drive = conf['do_drive'] == 'T'
cst_bc = conf['cst_bc'] == 'T'

D = 0
if do_diffusion:
    D = D + 1
if do_drive:
    D = D + (0.5 + 1/(np.exp(2*p_flip)-1))*p_drive

print("Diffusive move", do_diffusion)
print("Driven motion", do_drive)
print("D", D)

m = np.loadtxt(args.data_file[0])

for f in args.data_file[1:]:
    data = np.loadtxt(f)
    m += data

m /= len(args.data_file)

n_bins = m.shape[1]
xr = np.arange(n_bins)

def diffusion(x, x0, D, t):
    norm_factor = 1/(np.sqrt(4*np.pi*D*t) * (x[1]-x[0]))
    res = norm_factor * np.exp(-(x-x0)**2/(4*D*t))
    return res

if not args.show_d:
    ax1 = plt.subplot(211)
    ax2 = plt.subplot(212)
else:
    ax1 = plt.subplot(111)

for i, step in enumerate(m[::args.stride]):
    t = n_inner * dt * (1+i*args.stride)

    l, = ax1.plot(xr, step)

    if args.show_d:
        ax1.plot(xr, n_particles*diffusion(xr, n_bins//2-1, D, t), ls='--', color=l.get_color())

if args.show_rho0:
    plt.axhline(rho_0)

if args.show_n:
    n = m.sum(axis=1)
    t = (np.arange(len(n)))*n_inner*dt
    plt.figure()
    plt.plot(t, n)
    fit = np.polyfit(t[len(n)//2:], n[len(n)//2:], 1)
    plt.plot(t, np.poly1d(fit)(t))
    vel = fit[0]/(rho_0*D)
    if not cst_bc:
        vel /= 2
    print("Front velocity          :", vel)

if not args.show_d:
    for i, step in enumerate(m[::args.stride]):
        t = n_inner * dt * (1+i*args.stride)
        ax2.plot(xr - vel*t, step)

if args.no_plot:
    plt.show()
