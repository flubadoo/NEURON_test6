from neuron import h
import rxdfast as rxd
from matplotlib import pyplot

h.load_file('stdrun.hoc')

sec = h.Section()
sec.nseg = 1
sec.L = 100

r = rxd.Region([sec])
ca = rxd.Species(r, initial=lambda foo: 2)
na = rxd.Species(r, initial=lambda foo: 3)

# not at all biophysical, but just for a test
na_wave = rxd.Rate(na, na * na * na + 6 * ca)
ca_wave = rxd.Rate(ca, ca + na * ca * 2)

def plot_it():
    pyplot.subplot(2, 1, 1)
    pyplot.plot([seg.x * sec.L for seg in sec], ca.states)
    pyplot.subplot(2, 1, 2)
    pyplot.plot([seg.x * sec.L for seg in sec], na.states)

h.finitialize()
h.fadvance()
"""
for i in xrange(5):
    h.continuerun(i * 25)
    plot_it()

for plot_id in [1, 2]:
    pyplot.subplot(2, 1, plot_id)
    pyplot.xlim([0, sec.L])
    pyplot.ylim([0, 1.1])

print na.states
pyplot.show()
"""