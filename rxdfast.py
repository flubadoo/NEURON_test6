"""
rxdfast module

interconnection between Python and C for NEURON rxd v2
"""

from neuron import nrn_dll_sym, nrn_dll, h
import ctypes
import rxdmath
import itertools
import weakref

set_nonvint_block = nrn_dll_sym('set_nonvint_block')
nrn = nrn_dll()
dll = ctypes.cdll['./rxd.so']

fptr_prototype = ctypes.CFUNCTYPE(None)

set_nonvint_block(dll.rxd_nonvint_block)

set_setup = dll.set_setup
set_setup.argtypes = [fptr_prototype]
set_initialize = dll.set_initialize
set_initialize.argtypes = [fptr_prototype]

setup_solver = dll.setup_solver
setup_solver.argtypes = [ctypes.py_object, ctypes.py_object, ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_double)]

clear_rates = dll.clear_rates
register_rate = dll.register_rate

set_reaction_indices = dll.set_reaction_indices
set_reaction_indices.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_int), ctypes.POINTER(ctypes.c_int)]

states = None

region_sec = None
states = None
my_initializer = None

_all_species = []
_all_rates = []

class Node:
    def __init__(self, x, sec, region, species):
        self.x = x
        self.sec = sec
        self.region = region
        self.species = species

class Region:
    def __init__(self, secs):
        # TODO: remove this requirement (requires modifying do_setup)
        assert(len(secs) == 1)
        # TODO: don't keep hard references; this keeps them from being garbage
        #       collected (difficult because can't do weakrefs either, but doable
        #       because unique internal name)
        # make a copy
        self._secs = list(secs)

class Species(rxdmath._Arithmeticed):
    _ids = itertools.count(0)
    def __init__(self, region, d=0, initial=None):
        self._initial = initial
        self.d = d
        # TODO: change this when a species can live on more than one region
        self._region = [region]
        self._id = self._ids.next()
        self._str = '_species[%d]' % self._id
        # register the species
        # TODO: add callback to redo setup if species goes out of scope
        _all_species.append(weakref.ref(self))
    @property
    def indices(self):
        return range(self._start_i, self._stop_i)
    @property
    def states(self):
        # TODO: make this work even if not previously intialized
        # TODO: this is inefficient in that it gathers all states and then only
        #       returns its own
        return list(states)[self._start_i:self._stop_i]

class Rate:
    def __init__(self, species, rate):
        assert(isinstance(species, Species))
        self._rate = rxdmath._ensure_arithmeticed(rate)
        # TODO: add callback to redo setup if Rate goes out of scope
        _all_rates.append(weakref.ref(self))
        # TODO: this probably shouldn't keep species alive
        self._species = species
    
    def _compile(self):
        self._compiled = self._rate._compile()
        
def _list_to_cint_array(data):
    return (ctypes.c_int * len(data))(*tuple(data))

def _list_to_cdouble_array(data):
    return (ctypes.c_double * len(data))(*tuple(data))


def _num_states(species):
    # TODO: change this when a species can live on more than one region
    return sum(sec.nseg for sec in species._region[0]._secs)

def do_setup():
    # TODO: remove the global neighbor_* variables; replace with a copy in C
    global states, neighbor_index, neighbor_list, neighbor_rates, _all_change_states, _all_reaction_indices
    if not _all_species: return
    num_states = sum(_num_states(species()) for species in _all_species if species() is not None)
    states = h.Vector(num_states)
    # TODO: replace this to support arbitrary geometries, species on different regions
    neighbor_index = [0]
    neighbor_list = []
    neighbor_rates = []
    active_species = []
    node_i = 0
    species_count = 0
    for species in _all_species:
        # dereference and proceed only if species still alive
        species = species()
        if species is None: continue
        active_species.append(species)
        
        species._start_i = node_i
        
        # TODO: change this when supporting != 1 section
        # TODO: change this when allowing != 1 region
        sec = species._region[0]._secs[0]
        dx2 = (sec.L / float(sec.nseg)) ** 2
        # TODO: change this when allowing d to be nonuniform
        # TODO: change this when allowing non-cylindrical geometry
        rate = species.d / dx2;
        for i in xrange(sec.nseg):
            if 0 < i < sec.nseg - 1:
                neighbor_list += [node_i - 1, node_i + 1]
                neighbor_rates += [rate, rate]
            elif i == 0:
                # TODO: change this when allowing nseg = 1
                neighbor_list += [node_i + 1]
                neighbor_rates += [rate]
            else:
                # TODO: change this when allowing nseg = 1
                neighbor_list += [node_i - 1]
                neighbor_rates += [rate]
            # needs to be after we've modified neighbor_list
            neighbor_index.append(len(neighbor_list))
            node_i += 1

        species._stop_i = node_i
        species_count += 1

    neighbor_index = _list_to_cint_array(neighbor_index)
    neighbor_list = _list_to_cint_array(neighbor_list)
    neighbor_rates = _list_to_cdouble_array(neighbor_rates)
    setup_solver(states._ref_x[0], h._ref_dt, len(states), neighbor_index, neighbor_list, neighbor_rates)
    
    # now setup the reactions
    clear_rates()
    reaction_count = 0;
    _all_change_states = []
    _all_reaction_indices = []
    for rate in _all_rates:
        rate = rate()
        if rate is None: continue
        rate._compile()
        register_rate(rate._compiled)
        
        # TODO: this will need changed when more than one section/possibly different geometries
        change_states = _list_to_cint_array(rate._species.indices)
        indices = []
        for i in xrange(len(rate._species.indices)):
            indices += [species._start_i + i for species in active_species]
        
        indices = _list_to_cint_array(indices)
        
        _all_change_states.append(change_states)
        _all_reaction_indices.append(indices)

        set_reaction_indices(reaction_count, len(active_species), len(rate._species.indices), change_states, indices)
        reaction_count += 1

def do_initialize():
    """handle initialization at finitialize time"""
    node_i = 0
    for species in _all_species:
        species = species()
        if species is None: continue
        initializer = species._initial
        for region in species._region:
            for sec in region._secs:
                for seg in sec:
                    states.x[node_i] = initializer(Node(seg.x, sec, region, species))
                    node_i += 1

# register the Python callbacks
do_setup_fptr = fptr_prototype(do_setup)
do_initialize_fptr = fptr_prototype(do_initialize)
set_setup(do_setup_fptr)
set_initialize(do_initialize_fptr)
