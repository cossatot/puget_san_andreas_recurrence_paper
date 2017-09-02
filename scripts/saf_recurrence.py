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


def get_rec_ints(eqs, n_quakes=int(1e4), order_check='sort'):
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
    
    return rec_ints

    
def get_rec_pdf(eqs=None, rec_ints=None, n_quakes=int(1e4), **kwargs):
    if rec_ints is None:
        rec_ints = get_rec_ints(eqs, n_quakes, **kwargs)
    rec_int_pdf = RecKDE(rec_ints.ravel())
    rec_int_pdf.fit()

    return rec_int_pdf


def is_monotonic(array):
    return np.all(np.diff(array) >= 0)

print('making recurrence intervals')
ww_rec_ints = get_rec_ints(wrightwood_eqs, order_check='trim')
pc_rec_ints = get_rec_ints(pallet_creek_eqs, order_check='trim')
print('now doing other stuff')

ww_recs = get_rec_pdf(rec_ints=ww_rec_ints)
pc_recs = get_rec_pdf(rec_ints=pc_rec_ints)


def disp(eqs=None, rec_ints=None, n_quakes=int(1e4)):
    if rec_ints is None:
        rec_ints = get_rec_ints(eqs, n_quakes)
    return np.var(rec_ints) / np.mean(rec_ints)**2


def burstiness(eqs=None, rec_ints=None, n_quakes=int(1e4)):
    if rec_ints is None:
        rec_ints = get_rec_ints(eqs, n_quakes)

    return ((np.std(rec_ints) - np.mean(rec_ints)) 
            / (np.std(rec_ints) + np.mean(rec_ints)))


def memory(eqs=None, rec_ints=None, n_quakes=int(1e4)):
    if rec_ints is None:
        rec_ints = get_rec_ints(eqs, n_quakes)

    n = len(rec_ints)
    m = rec_ints.mean()
    v = rec_ints.var()

    return (1 / (n-1)) * np.sum(((rec_ints[i]-m) * (rec_ints[i+1] - m)
                                 for i in range(n-1))) / v


ww_Bs = [burstiness(eqs=wrightwood_eqs, rec_ints=wri)
         for i, wri in enumerate(ww_rec_ints)]
pc_Bs = [burstiness(eqs=pallet_creek_eqs, rec_ints=pri)
         for i, pri in enumerate(pc_rec_ints)]

ww_Ms = [memory(eqs=wrightwood_eqs, rec_ints=wri)
         for i, wri in enumerate(ww_rec_ints)]
pc_Ms = [memory(eqs=pallet_creek_eqs, rec_ints=pri)
         for i, pri in enumerate(pc_rec_ints)]



plt.style.use('dark_background')
mpl.rcParams['font.size'] = 22



plt.figure()
plt.scatter(ww_Ms, ww_Bs, c='b', alpha=0.1, lw=0, s=10)
plt.scatter(pc_Ms, pc_Bs, c='r', alpha=0.1, lw=0, s=10)
plt.xlim([-1,1])
plt.ylim([-1,1])
plt.xlabel('Memory')
plt.ylabel('Burstiness')





fig = plt.figure(figsize=(15,10))
gs = mpl.gridspec.GridSpec(2,3, height_ratios=[1.5, 2.5])
gs.update(left=0.07, right=0.92, wspace=0.2, hspace=0.4, top=0.95, bottom=0.09)

axs1 = [0, 1, 2]
axs1[0] = plt.subplot(gs[0, :])
axs1[1] = plt.subplot(gs[1, :])
#axs1[2] = plt.subplot(gs[1, 2:3])


saf_rec_y_max = 0.000013
saf_rec_zoom_x_max = 200

axs1[0].set_ylim([0,0.5])
axs1[1].set_xlim([0,500])
#axs1[0].invert_xaxis()
#axs1[1].set_ylim([0.,saf_rec_y_max])
#axs1[2].set_ylim([0.,saf_rec_y_max])
#axs1[2].set_xlim([0,saf_rec_zoom_x_max])

for ev in ww_events:
    axs1[0].plot(wrightwood_eqs[ev].age, wrightwood_eqs[ev].age_err, 'c')

for ev in pc_events:
    axs1[0].plot(pallet_creek_eqs[ev].age, pallet_creek_eqs[ev].age_err, 'm')


axs1[0].set_xlabel('Earthquake time, calendar years BP')
axs1[1].set_xlabel('Earthquake recurrence interval (years)')

axs11y = axs1[1].twinx()
#axs13y = axs1[2].twinx()

# replace rec_pdfs.sum with ww_recs, etc.
axs1[1].plot(ww_recs.x, ww_recs.px,
           #'--', 
           color='c',
           lw=3,
           label='Wrightwood PDF',
           )

#axs1[2].plot(ww_recs.x, ww_recs.px,
#           #'--', 
#           color='c',
#           lw=3,
#           )

axs11y.plot(ww_recs.x, np.cumsum(ww_recs.px) / ww_recs.px.sum(),
           'c--',
           lw=3,
           label='Wrightwood CDF',
           )

#axs13y.plot(ww_recs.x, np.cumsum(ww_recs.px) / ww_recs.px.sum(),
#           'c--',
#           lw=3,
#           )

# pallet creek
axs1[1].plot(pc_recs.x, pc_recs.px,
           #'--', 
           color='m',
           lw=3,
           label='Pallet Creek PDF',
           )

#axs1[2].plot(pc_recs.x, pc_recs.px,
#           #'--', 
#           color='m',
#           lw=3,
#           )

axs11y.plot(pc_recs.x, np.cumsum(pc_recs.px) / pc_recs.px.sum(),
           'm--',
           lw=3,
           label='Pallet Creek CDF',
           )

#axs13y.plot(pc_recs.x, np.cumsum(pc_recs.px) / pc_recs.px.sum(),
#           'm--',
#           lw=3,
#           )

axs1[0].set_title('San Andreas Fault EQ times')
axs1[1].set_title('San Andreas Fault recurrence PDF')#, fontsize=12)

axs1[1].legend(loc='upper right', fontsize=22)


f, (ax0, ax1) = plt.subplots(2, sharex=True, figsize=(15,10))

ax0.plot(ww_recs.x[:-1], cp.recurrence.hazard(ww_recs.x[:-1], ww_recs), 'c',
         lw=3, label='Wrightwood')
ax0.plot(pc_recs.x[:-1], cp.recurrence.hazard(pc_recs.x[:-1], pc_recs), 'm',
         lw=3, label='Pallet Creek')

ax1.plot(ww_recs.x[:-1], cp.recurrence.hazard(ww_recs.x[:-1], ww_recs), 
         'c', lw=3)
ax1.plot(pc_recs.x[:-1], cp.recurrence.hazard(pc_recs.x[:-1], pc_recs), 
         'm', lw=3)

ax0.set_ylim([0., 0.15])
ax1.set_ylim([1e-3, 10])

ax1.set_yscale('log')
ax1.set_yscale('log')

ax0.legend(loc='upper left', fontsize=22)

ax1.set_xlabel('years since last earthquake')
#ax1.set_ylabel('instantaneous earthquake probability')

plt.subplots_adjust(left=0.10, right=0.97, 
                    wspace=0.45, hspace=0.2, 
                    top=0.95, bottom=0.10)


# real years
f, (ax2, ax3) = plt.subplots(2, sharex=True, figsize=(15,10))

ax2.plot(1857+ ww_recs.x[:-1], cp.recurrence.hazard(ww_recs.x[:-1], ww_recs), 'c',
         lw=3, label='Wrightwood')
ax2.plot(1857+ pc_recs.x[:-1], cp.recurrence.hazard(pc_recs.x[:-1], pc_recs), 'm',
         lw=3, label='Pallet Creek')

ax2.axvline(2017, linestyle='--', color='grey')

ax3.plot(1857+ ww_recs.x[:-1], cp.recurrence.hazard(ww_recs.x[:-1], ww_recs), 
         'c', lw=3)
ax3.plot(1857+ pc_recs.x[:-1], cp.recurrence.hazard(pc_recs.x[:-1], pc_recs), 
         'm', lw=3)

ax3.axvline(2017, linestyle='--', color='grey')

ax2.set_ylim([0., 0.15])
ax3.set_ylim([1e-3, 10])

ax3.set_yscale('log')
ax3.set_yscale('log')

ax2.legend(loc='upper left', fontsize=22)

ax3.set_xlabel('year')
#ax3.set_ylabel('instantaneous earthquake probability')

plt.subplots_adjust(left=0.10, right=0.97, 
                    wspace=0.45, hspace=0.2, 
                    top=0.95, bottom=0.10)

plt.show()
