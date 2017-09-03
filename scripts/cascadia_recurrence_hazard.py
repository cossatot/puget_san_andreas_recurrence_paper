import json
from collections import OrderedDict
import numpy as np
import matplotlib.pyplot as plt

import sys; sys.path.append('/Users/itchy/research/culpable/')
import culpable as cp

from culpable.offset_marker import OffsetMarker
from culpable.recurrence import RecKDE

time_data_file = '../data/puget_sound_eq_probs.json'
fault_zones_file = '../data/puget_sound_fault_zones.json'

# Load data
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

# Process data
eqs = {k: make_om(k, v) for k, v in eq_times.items()}

rec_fault_zones = {k : v for k, v in fault_zones.items() if len(v) > 1}

sfz_eqs = {k: v for k, v in eqs.items() if k in rec_fault_zones['Seattle']}


# Calculate mean ages
eq_mean_ages = {k: cp.stats.Pdf(eq.age, eq.age_err).mean() 
                for k, eq in eqs.items()}

eq_mean_ages = OrderedDict(sorted(eq_mean_ages.items(), key=lambda x: x[1]))

sfz_mean_ages = {k: cp.stats.Pdf(eq.age, eq.age_err).mean() 
                 for k, eq in sfz_eqs.items()}
sfz_mean_ages = OrderedDict(sorted(sfz_mean_ages.items(), key=lambda x: x[1]))


# Do recurrence interval calculations
def get_rec_ints(fz, n_quakes=int(1e4), ravel=True):
    
    eq_list = [eq for k, eq in eqs.items() if k in fz]

    eq_times = np.array([eq.sample_ages(n_quakes) for eq in eq_list]).T
    rec_ints = np.diff(np.sort(eq_times, axis=1), axis=1)

    if ravel == True:
        rec_ints = rec_ints.ravel()
    else:
        rec_ints = rec_ints.T
    
    return rec_ints

    
def get_rec_pdf(fz, n_quakes=int(1e4)):
    rec_ints = get_rec_ints(fz, n_quakes)
    rec_int_pdf = RecKDE(rec_ints)

    return rec_int_pdf


sfz_rec_ints = get_rec_ints(rec_fault_zones['Seattle'], ravel=False)
sfz_rec_pdfs = [RecKDE(row) for row in sfz_rec_ints]

sfz_tot_rec_pdf = get_rec_pdf(rec_fault_zones['Seattle'])



pug_rec_ints = get_rec_ints(rec_fault_zones['puget_lowland'], ravel=False)
pug_rec_pdfs = [RecKDE(row) for row in pug_rec_ints]

pug_tot_rec_pdf = get_rec_pdf(rec_fault_zones['puget_lowland'])


print('Puget mode: ', pug_tot_rec_pdf.mode(), '\n',
      'Puget median: ', pug_tot_rec_pdf.median(), '\n',
      'Puget mean: ', pug_tot_rec_pdf.mean(), '\n',
      )

print('Seattle mode: ', sfz_tot_rec_pdf.mode(), '\n',
      'Seattle median: ', sfz_tot_rec_pdf.median(), '\n',
      'Seattle mean: ', sfz_tot_rec_pdf.mean(), '\n',
      )


## Plotting

#####
# SFZ
#####

f0, (ax00, ax01) = plt.subplots(2, 1, figsize=(7,4))
f0.subplots_adjust(top=0.95, bottom=0.1, hspace=0.4)

for k, eq in sfz_eqs.items():
    ax00.plot(eq.age, eq.age_err,
              lw=0.5)

ax00.invert_xaxis()
ax00.set_xlim([16000, 0])
ax00.set_ylim([0., 0.02])
ax00.set_xlabel('calendar year')
ax00.set_ylabel('eq time probability')
ax00.set_title('a', loc='left', weight='bold')

ax02 = ax01.twinx()
ax02.plot(sfz_tot_rec_pdf.x, sfz_tot_rec_pdf.y,
          color='grey', lw=1)

ax02.plot(sfz_tot_rec_pdf.cdf.x, 
          sfz_tot_rec_pdf.cdf.y * sfz_tot_rec_pdf.mode()[1],
          color='grey', linestyle='--', lw=1)

for rec_pdf in sfz_rec_pdfs:
    ax01.plot(rec_pdf.x, rec_pdf.y,
              lw=0.5)
#ax00.set_title('Seattle Fault Zone recurrence')

ax01.set_title('b', loc='left', weight='bold')
ax01.set_xlabel('recurrence interval')
ax01.set_ylabel('probability')
ax02.set_ylabel('probability')

f0.savefig('../manuscript/figures/sfz_recurrence.pdf')

#####
# PUG
#####

f1, (ax10, ax11) = plt.subplots(2, 1, figsize=(7,4))
f1.subplots_adjust(top=0.95, bottom=0.1, hspace=0.4)

for k, eq in eqs.items():
    ax10.plot(eq.age, eq.age_err,
              lw=0.5)

ax10.invert_xaxis()
ax10.set_xlim([16000, 0])
ax10.set_ylim([0., 0.02])
ax10.set_xlabel('calendar year')
ax10.set_ylabel('eq time probability')
ax10.set_title('a', loc='left', weight='bold')

ax12 = ax11.twinx()
ax12.plot(pug_tot_rec_pdf.x, pug_tot_rec_pdf.y,
          color='grey', lw=1)

ax12.plot(pug_tot_rec_pdf.cdf.x, 
          pug_tot_rec_pdf.cdf.y * pug_tot_rec_pdf.mode()[1],
          color='grey', linestyle='--', lw=1)

for rec_pdf in pug_rec_pdfs:
    ax11.plot(rec_pdf.x, rec_pdf.y,
              lw=0.5)
#ax10.set_title('Puget Lowland recurrence')

ax11.set_xlabel('recurrence interval')
ax11.set_ylabel('probability')
ax12.set_ylabel('probability')
ax11.set_title('b', loc='left', weight='bold')


f1.savefig('../manuscript/figures/pug_recurrence.pdf')



f2, (ax20, ax21, ax22) = plt.subplots(3, 1, figsize=(7,6))
f2.subplots_adjust(top=0.95, bottom=0.1, hspace=0.4)

ax20.plot(pug_tot_rec_pdf.x, cp.recurrence.hazard(pug_tot_rec_pdf.x,
                                                  pug_tot_rec_pdf),
          label='Puget Lowland')

ax20.plot(sfz_tot_rec_pdf.x, cp.recurrence.hazard(sfz_tot_rec_pdf.x,
                                                  sfz_tot_rec_pdf),
          label='Seattle Fault Zone')

ax20.axhline(1 / pug_tot_rec_pdf.mean(), color='C0', linestyle='--', lw=0.5)
ax20.axhline(1 / sfz_tot_rec_pdf.mean(), color='C1', linestyle='--', lw=0.5)

ax20.legend(loc='best')
ax20.set_xlabel('years since last event')
ax20.set_ylabel('earthquake hazard λ(t)')
ax20.set_title('a', loc='left', weight='bold')


ax21.plot(pug_tot_rec_pdf.x + (1950 - list(eq_mean_ages.values())[0]), 
          cp.recurrence.hazard(pug_tot_rec_pdf.x, pug_tot_rec_pdf),
          )

ax21.axvline(1950 - list(eq_mean_ages.values())[0], 
             color='grey', lw=0.5, label='last earthquake')

ax21.axvline(2017, linestyle='--', color='grey', lw=0.5)
ax21.axhline(1 / pug_tot_rec_pdf.mean(), 
             color='C0', linestyle='--', lw=0.5, label='mean (Poisson) hazard')

ax21.set_xlabel('calendar year')
ax21.set_ylabel('earthquake hazard λ(t)')
ax21.set_title('b', loc='left', weight='bold')

ax21.legend(loc='upper right')

ax22.plot(sfz_tot_rec_pdf.x + (1950 - list(sfz_mean_ages.values())[0]), 
          cp.recurrence.hazard(sfz_tot_rec_pdf.x, sfz_tot_rec_pdf),
          color='C1')

ax22.axvline(1950 - list(sfz_mean_ages.values())[0], 
             color='grey', lw=0.5, label='last earthquake')
ax22.axvline(2017, linestyle='--', color='grey', lw=0.5)
ax22.axhline(1 / sfz_tot_rec_pdf.mean(), 
             color='C1', linestyle='--', lw=0.5, label='mean (Poisson) hazard')

ax22.set_xlabel('calendar year')
ax22.set_ylabel('earthquake hazard λ(t)')
ax22.set_title('c', loc='left', weight='bold')
ax22.legend(loc='upper left')

f2.savefig('../manuscript/figures/pug_hazard.pdf')




#f3, (ax30, ax31, ax32) = plt.subplots(3, 1, figsize=(7,6))
#f3.subplots_adjust(top=0.95, bottom=0.1, hspace=0.4)

#ax30.plot(pug_tot_rec_pdf.x, cp.recurrence.S(pug_tot_rec_pdf.x,
#                                             pug_tot_rec_pdf),
#          label='Puget Lowland')

#ax30.plot(sfz_tot_rec_pdf.x, cp.recurrence.S(sfz_tot_rec_pdf.x,
#                                             sfz_tot_rec_pdf),
#          label='Seattle Fault Zone')


plt.show()
