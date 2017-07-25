from collections import OrderedDict

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

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

plt.figure()

for ev in ww_events:
    plt.plot(wrightwood_eqs[ev].age, wrightwood_eqs[ev].age_err, 'b')

for ev in pc_events:
    plt.plot(pallet_creek_eqs[ev].age, pallet_creek_eqs[ev].age_err, 'r')



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

    n = len(eqs) - 1
    m = rec_ints.mean()
    v = rec_ints.var()

    return (1 / (n-1)) * np.sum(((rec_ints[i]-m) * (rec_ints[i+1] - m)
                                 for i in range(n-1))) / v


ww_B = burstiness(eqs=wrightwood_eqs, rec_ints=ww_rec_ints.ravel())
pc_B = burstiness(eqs=pallet_creek_eqs, rec_ints=pc_rec_ints.ravel())

ww_M = memory(eqs=wrightwood_eqs, rec_ints=ww_rec_ints.ravel())
pc_M = memory(eqs=pallet_creek_eqs, rec_ints=pc_rec_ints.ravel())

ww_Bs = [burstiness(eqs=wrightwood_eqs, rec_ints=wri)
         for i, wri in enumerate(ww_rec_ints)]
pc_Bs = [burstiness(eqs=pallet_creek_eqs, rec_ints=pri)
         for i, pri in enumerate(pc_rec_ints)]

ww_Ms = [memory(eqs=wrightwood_eqs, rec_ints=wri)
         for i, wri in enumerate(ww_rec_ints)]
pc_Ms = [memory(eqs=pallet_creek_eqs, rec_ints=pri)
         for i, pri in enumerate(pc_rec_ints)]



plt.figure()
plt.scatter(ww_Ms, ww_Bs, c='b', alpha=0.1, lw=0, s=10)
plt.scatter(pc_Ms, pc_Bs, c='r', alpha=0.1, lw=0, s=10)
plt.xlim([-1,1])
plt.ylim([-1,1])
plt.xlabel('Memory')
plt.ylabel('Burstiness')


print('wrightwood IDI: {}'.format(disp(wrightwood_eqs)))
print('pallet creek IDI: {}'.format(disp(pallet_creek_eqs)))

print('wrightwood B: {}'.format(ww_B))
print('pallet creek B: {}'.format(pc_B))


plt.figure()
plt.plot(ww_recs.x, ww_recs.px, 'b')
plt.plot(pc_recs.x, pc_recs.px, 'r')

plt.show()
