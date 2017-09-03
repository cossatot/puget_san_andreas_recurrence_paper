from collections import OrderedDict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

import sys; sys.path.append('/Users/itchy/research/culpable/')
import culpable as cp

from culpable.offset_marker import OffsetMarker
from culpable.recurrence import RecKDE


wrightwood_data = pd.read_csv('../data/WW_ages.csv')
pallet_creek_data = pd.read_csv('../data/PC_Scharer2011_ages.csv')

# remove zero probability rows
#wrightwood_data = wrightwood_data[wrightwood_data.probs > 0.]
#pallet_creek_data = pallet_creek_data[pallet_creek_data.probs > 0.]

def om_from_dataframe(eq_name, dataframe):
    event_dataframe = dataframe.loc[dataframe.event == eq_name]

    if event_dataframe.shape[0] > 1:
        om = OffsetMarker(name=eq_name,
                          age=event_dataframe.age.values,
                          age_err=event_dataframe.probs.values,
                          age_dist_type='arbitrary')
    
    elif event_dataframe.shape[0] == 1:
        om = OffsetMarker(name=eq_name,
                          age=event_dataframe.age.values,
                          age_err=event_dataframe.probs.values,
                          age_dist_type='scalar')
    om.init()
    return om


ww_remove = []
#ww_remove = ['PCT', 'W520', 'W380', 'W390']
ww_events = []
for ev in wrightwood_data.event.values:
    if ev not in ww_events:
        if ev not in ww_remove:
            ww_events.append(ev) 

pc_events = []
for ev in pallet_creek_data.event.values:
    if ev not in pc_events:
        pc_events.append(ev) 



wrightwood_eqs = OrderedDict((eq, om_from_dataframe(eq, wrightwood_data))
                  for eq in ww_events)

pallet_creek_eqs = OrderedDict((eq, om_from_dataframe(eq, pallet_creek_data))
                    for eq in pc_events)


# get mean ages
ww_mean_ages = {k: cp.stats.Pdf(eq.age, eq.age_err).mean()
                for k, eq in wrightwood_eqs.items()}

ww_mean_ages = OrderedDict(sorted(ww_mean_ages.items(), key=lambda x: x[1]))

pc_mean_ages = {k: cp.stats.Pdf(eq.age, eq.age_err).mean()
                for k, eq in pallet_creek_eqs.items()}

pc_mean_ages = OrderedDict(sorted(pc_mean_ages.items(), key=lambda x: x[1]))


# recurrence interval stuff
def get_rec_ints(eqs, n_quakes=int(5e3), order_check='sort', ravel=False):
    eq_times = np.array([eq.sample_ages(n_quakes) for eq in eqs.values()]).T

    if order_check == 'sort':
        eq_times_sort = np.sort(eq_times, axis=1)

    elif order_check == 'trim':
        eq_times_sort = eq_times.copy()
        for i, row in enumerate(eq_times):
            if ~is_monotonic(row):
                while ~is_monotonic(row):
                    row = np.array([eq.sample_ages(1) for eq in eqs.values()])
            eq_times_sort[i,:] = row.T

    rec_ints = np.diff(eq_times_sort, axis=1)

    if ravel == True:
        rec_ints = rec_ints.ravel()
    else:
        rec_ints = rec_ints.T
    
    return rec_ints

    
def get_rec_pdf(eqs=None, rec_ints=None, n_quakes=int(1e4), **kwargs):
    if rec_ints is None:
        rec_ints = get_rec_ints(eqs, n_quakes, **kwargs)
    rec_int_pdf = RecKDE(rec_ints)

    return rec_int_pdf


def is_monotonic(array):
    return np.all(np.diff(array) >= 0)

print('making recurrence intervals')
ww_rec_ints = get_rec_ints(wrightwood_eqs, order_check='trim', ravel=False)
pc_rec_ints = get_rec_ints(pallet_creek_eqs, order_check='trim', ravel=False)
print('now doing other stuff')

ww_rec_pdfs = [RecKDE(row) for row in ww_rec_ints]
pc_rec_pdfs = [RecKDE(row) for row in pc_rec_ints]

ww_tot_rec_pdf = get_rec_pdf(rec_ints=ww_rec_ints.ravel())
pc_tot_rec_pdf = get_rec_pdf(rec_ints=pc_rec_ints.ravel())


print('Wrightwood mode: ', ww_tot_rec_pdf.mode(), '\n',
      'WrighPuget median: ', ww_tot_rec_pdf.median(), '\n',
      'WrighPuget mean: ', ww_tot_rec_pdf.mean(), '\n',
      )

print('Pallet Creek mode: ', pc_tot_rec_pdf.mode(), '\n',
      'Pallet Creek median: ', pc_tot_rec_pdf.median(), '\n',
      'Pallet Creek mean: ', pc_tot_rec_pdf.mean(), '\n',
      )
# plots to plot

# ages, recurrence intervals
# 3: ages for all, recs for PC, recs for WW
f0, (ax00, ax01, ax03) = plt.subplots(3, 1, figsize=(7,6))
f0.subplots_adjust(top=0.95, bottom=0.1, hspace=0.4)

for k, eq in wrightwood_eqs.items():
    ax00.plot(eq.age, eq.age_err,
              lw=0.5, color='c')

for k, eq in pallet_creek_eqs.items():
    ax00.plot(eq.age, eq.age_err,
              lw=0.5, color='m')

#ax00.set_xlim([16000, 0])
#ax00.set_ylim([0., 0.02])
ax00.set_xlabel('calendar year')
ax00.set_ylabel('eq time probability')
ax00.set_title('a', loc='left', weight='bold')

ax02 = ax01.twinx()
ax02.plot(ww_tot_rec_pdf.x, ww_tot_rec_pdf.y,
          color='grey', lw=1)

ax02.plot(ww_tot_rec_pdf.cdf.x, 
          ww_tot_rec_pdf.cdf.y * ww_tot_rec_pdf.mode()[1],
          color='grey', linestyle='--', lw=1)

# wrightwood recurrence plot
for rec_pdf in ww_rec_pdfs:
    ax01.plot(rec_pdf.x, rec_pdf.y,
              lw=0.5)
ax01.set_ylim([0., 0.025])

ax01.set_title('b', loc='left', weight='bold')
ax01.set_xlabel('recurrence interval')
ax01.set_ylabel('probability')
ax02.set_ylabel('probability')

ax02 = ax01.twinx()
ax02.plot(ww_tot_rec_pdf.x, ww_tot_rec_pdf.y,
          color='grey', lw=1)

ax02.plot(ww_tot_rec_pdf.cdf.x, 
          ww_tot_rec_pdf.cdf.y * ww_tot_rec_pdf.mode()[1],
          color='grey', linestyle='--', lw=1)

# pallet creek recurrence plot
for rec_pdf in pc_rec_pdfs:
    ax03.plot(rec_pdf.x, rec_pdf.y,
              lw=0.5)
ax03.set_ylim([0., 0.025])

ax03.set_title('c', loc='left', weight='bold')
ax03.set_xlabel('recurrence interval')
ax03.set_ylabel('probability')
ax02.set_ylabel('probability')

ax04 = ax03.twinx()
ax04.plot(pc_tot_rec_pdf.x, pc_tot_rec_pdf.y,
          color='grey', lw=1)

ax04.plot(pc_tot_rec_pdf.cdf.x, 
          pc_tot_rec_pdf.cdf.y * pc_tot_rec_pdf.mode()[1],
          color='grey', linestyle='--', lw=1)

f0.savefig('../manuscript/figures/saf_recurrence.pdf')


# hazard
# 3: generic hazard, ww from last date, pc from last date

f2, (ax20, ax21, ax22) = plt.subplots(3, 1, figsize=(7,6))
f2.subplots_adjust(top=0.95, bottom=0.1, hspace=0.4)

ax20.plot(ww_tot_rec_pdf.x, cp.recurrence.hazard(ww_tot_rec_pdf.x,
                                                  ww_tot_rec_pdf),
          label='Wrightwood')

ax20.plot(pc_tot_rec_pdf.x, cp.recurrence.hazard(pc_tot_rec_pdf.x,
                                                  pc_tot_rec_pdf),
          label='Pallet Creek')

ax20.axhline(1 / ww_tot_rec_pdf.mean(), color='C0', linestyle='--', lw=0.5)
ax20.axhline(1 / pc_tot_rec_pdf.mean(), color='C1', linestyle='--', lw=0.5)

ax20.set_ylim([0., 0.04])

ax20.legend(loc='best')
ax20.set_xlabel('years since last event')
ax20.set_ylabel('earthquake hazard λ(t)')
ax20.set_title('a', loc='left', weight='bold')


ax21.plot(ww_tot_rec_pdf.x + list(ww_mean_ages.values())[-1], 
          cp.recurrence.hazard(ww_tot_rec_pdf.x, ww_tot_rec_pdf),
          color='c')

ax21.axvline(list(ww_mean_ages.values())[-1], 
             color='grey', lw=0.5, label='last earthquake')

ax21.axvline(2017, linestyle='--', color='grey', lw=0.5)
ax21.axhline(1 / ww_tot_rec_pdf.mean(), 
             color='c', linestyle='--', lw=0.5, label='mean (Poisson) hazard')

ax21.set_ylim([0., 0.04])

ax21.set_xlabel('calendar year')
ax21.set_ylabel('earthquake hazard λ(t)')
ax21.set_title('b', loc='left', weight='bold')

ax21.legend(loc='best')

ax22.plot(pc_tot_rec_pdf.x + list(pc_mean_ages.values())[-1], 
          cp.recurrence.hazard(pc_tot_rec_pdf.x, pc_tot_rec_pdf),
          color='m')

ax22.axvline(list(pc_mean_ages.values())[-1], 
             color='grey', lw=0.5, label='last earthquake')
ax22.axvline(2017, linestyle='--', color='grey', lw=0.5)
ax22.axhline(1 / pc_tot_rec_pdf.mean(), 
             color='m', linestyle='--', lw=0.5, label='mean (Poisson) hazard')

ax22.set_ylim([0., 0.04])

ax22.set_xlabel('calendar year')
ax22.set_ylabel('earthquake hazard λ(t)')
ax22.set_title('c', loc='left', weight='bold')
ax22.legend(loc='best')

f2.savefig('../manuscript/figures/saf_hazard.pdf')


plt.show()
