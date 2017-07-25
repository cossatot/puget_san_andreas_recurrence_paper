import json
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib as mpl

import sys; sys.path.append('/Users/itchy/research/culpable/')
import culpable as cp

from culpable.offset_marker import OffsetMarker
from culpable.recurrence import RecKDE


time_data_file = '../data/puget_sound_eq_probs.json'
fault_zones_file = '../data/puget_sound_fault_zones.json'

with open(time_data_file, 'r') as f:
    eq_times = json.load(f)

with open(fault_zones_file, 'r') as f:
    fault_zones = json.load(f)

fault_zones['puget_lowland'] = []
for fz, eqs in fault_zones.items():
    if fz != 'puget_lowland':
        for eq in eqs:
            fault_zones['puget_lowland'].append(eq)

def year_to_cal_year(eq_year_list):
    cal_year = [(-1 * (x - 1950)) for x in eq_year_list]
    return cal_year


def make_om(eq_key, eq):
    om = OffsetMarker(name=eq_key,
                      age=np.array(year_to_cal_year(eq['year'])),
                      age_err=np.array(eq['probs']),
                      age_dist_type='arbitrary'
                      )
    om.init()
    om.trim_ages(min=-66)
    return om

eqs = {k: make_om(k, v) for k, v in eq_times.items()}


rec_fault_zones = {k : v for k, v in fault_zones.items() if len(v) > 1}





def get_rec_ints(fz, n_quakes=int(1e4)):
    
    eq_list = [eq for k, eq in eqs.items() if k in fz]

    eq_times = np.array([eq.sample_ages(n_quakes) for eq in eq_list]).T
    rec_ints = np.diff(np.sort(eq_times, axis=1), axis=1).ravel()
    
    return rec_ints

    
def get_rec_pdf(fz, n_quakes=int(1e4)):
    rec_ints = get_rec_ints(fz, n_quakes)
    rec_int_pdf = RecKDE(rec_ints)
    rec_int_pdf.fit()

    return rec_int_pdf


recs = {k : get_rec_pdf(v) for k, v in rec_fault_zones.items()}

rec_means = {k: cp.stats.pdf_mean(v.x, v.px) for k, v in recs.items()}


def disp(fz, n_quakes=int(1e4)):
    rec_ints = get_rec_ints(fz, n_quakes)
    return np.var(rec_ints) / np.mean(rec_ints)**2


def burstiness(fz, n_quakes=int(1e4)):
    rec_ints = get_rec_ints(fz, n_quakes)
    return ((np.std(rec_ints) - np.mean(rec_ints)) 
            / (np.std(rec_ints) + np.mean(rec_ints)))


def memory(fz, n_quakes=int(1e4)):
    rec_ints = get_rec_ints(fz, n_quakes)

    n = len(eqs)
    m = rec_ints.mean()
    v = rec_ints.var()

    return (1 / (n-1)) * np.sum(((rec_ints[i]-m) * (rec_ints[i+1] - m)
                                 for i in range(n-1))) / v



disps = {k: disp(v) for k, v in rec_fault_zones.items()}
Bs = {k: burstiness(v) for k, v in rec_fault_zones.items()}
Ms = {k: memory(v) for k, v in rec_fault_zones.items()}

f, axs = plt.subplots(nrows=9, sharex=True, figsize=(8,12))

for i, (eq_name, rec_pdf) in enumerate(recs.items()):
    axs[i].plot(rec_pdf.x, rec_pdf.px, label=eq_name)

    axs[i].legend()
plt.show()



#def get_cluster(fz


#plt.plot(sfz_rec.x, sfz_rec.px)
#plt.show()
