from neuron import h
import rxdfast as rxd
from matplotlib import pyplot

h.load_file('stdrun.hoc')

sec = h.Section()
sec.nseg = 5
sec.L = 100

r = rxd.Region([sec])
na = rxd.Species(r, initial=lambda foo: 1)
ca = rxd.Species(r, initial=lambda foo: 0)

# not at all biophysical, but just for a test
na_wave = rxd.Rate(na, -ca)
ca_wave = rxd.Rate(ca, na)

def plot_it():
    pyplot.subplot(2, 1, 1)
    pyplot.plot([seg.x * sec.L for seg in sec], ca.states)
    pyplot.subplot(2, 1, 2)
    pyplot.plot([seg.x * sec.L for seg in sec], na.states)

h.finitialize()
h.fadvance()
h.continuerun(1.57)

print na.states
print ca.states
