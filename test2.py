"""
test-location-based-approach module

prototype for parts of NEURON rxd v2
"""

from neuron import nrn_dll_sym, nrn_dll, h
import ctypes
import rxdmath
import itertools
import weakref

nseg = 5
L = 5 # length of section
num_equations = 2
d = 0 # diffusion rate

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
set_reaction_indices.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.POINTER(ctypes.c_int)]

states = h.Vector(nseg * num_equations)
        
def _list_to_cint_array(data):
    return (ctypes.c_int * len(data))(*tuple(data))

def _list_to_cdouble_array(data):
    return (ctypes.c_double * len(data))(*tuple(data))


def cossinreaction():
    formula = """
    void reaction(double* states, double* rhs) {
        rhs[0] = - _states[1];
        rhs[1] = states[0];
    }
    """
    filename = 'rxd-' + str(uuid.uuid1())
    with open(filename + '.c', 'w') as f:
        f.write(formula)
    os.system('gcc -I/usr/include/python2.7 -lpython2.7 -shared -o %s.so -fPIC %s.c' % (filename, filename))
    dll = ctypes.cdll['./%s.so' % filename]  
    self._filename = filename
    reaction = dll.reaction
    reaction.argtypes = [ctypes.POINTER(ctypes.c_double)] * (num_equations)
    reaction.restype = ctypes.c_double
    return reaction

def do_setup():
    ###################################################################
    #
    # setup the diffusion
    #
    ###################################################################

    dx2 = (L / float(nseg)) ** 2
    # TODO: change this when allowing d to be nonuniform
    # TODO: change this when allowing non-cylindrical geometry
    rate = d / dx2;
    neighbor_list = []
    neighbor_rates = []

    for i in xrange(nseg):
        if 0 < i < nseg - 1:
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

    neighbor_index = _list_to_cint_array(neighbor_index)
    neighbor_list = _list_to_cint_array(neighbor_list)
    neighbor_rates = _list_to_cdouble_array(neighbor_rates)
    setup_solver(states._ref_x[0], h._ref_dt, len(states), neighbor_index, neighbor_list, neighbor_rates)

    ###################################################################
    #
    # setup the reaction
    #
    ###################################################################
    
    # get rid of any old stuff
    clear_rates()

    # compile and transfer the kinetics formula
    register_rate(cossinreaction())

    # now take care of the indices etc
    reaction_index = 0
    num_species_per_location = 2
    num_locations = nseg

    # used to gather all states for a given location together
    location_state_indices = _list_to_cint_array(list(itertools.chain.from_iterable([[i, i + nseg] for i in xrange(nseg)])))

    # finally transfer the reaction metadata
    set_reaction_indices(reaction_index, num_species_per_location, num_locations, location_state_indices)


def do_initialize():
    """handle initialization at finitialize time"""
    for i in xrange(len(states)):
        if i < nseg:
            states.x[i] = 1
        else:
            states.x[i] = 0

# register the Python callbacks
do_setup_fptr = fptr_prototype(do_setup)
do_initialize_fptr = fptr_prototype(do_initialize)
set_setup(do_setup_fptr)
set_initialize(do_initialize_fptr)


h.finitialize()
h.continuerun(1.57)

print 'c:'
print ', '.join(str(v) for i, v in enumerate(states) if i < nseg)
print 's:'
print ', '.join(str(v) for i, v in enumerate(states) if i >= nseg)
