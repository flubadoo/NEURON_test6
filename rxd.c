#include <stdio.h>
#include <assert.h>
#include <Python.h>
#include "matrix/matrix.h"

typedef struct {
	PyObject_HEAD
	void* ho_;
	union {
		double x_;
		char* s_;
		void* ho_;
		double* px_;
	}u;
	void* sym_; // for functions and arrays
	void* iteritem_; // enough info to carry out Iterator protocol
	int nindex_; // number indices seen so far (or narg)
	int* indices_; // one fewer than nindex_
	int type_; // 0 HocTopLevelInterpreter, 1 HocObject
		// 2 function (or TEMPLATE)
		// 3 array
		// 4 reference to number
		// 5 reference to string
		// 6 reference to hoc object
		// 7 forall section iterator
		// 8 pointer to a hoc scalar
		// 9 incomplete pointer to a hoc array (similar to 3)
} PyHocObject;

typedef void (*fptr)(void);
typedef double (*ReactionRate)(double*);

/*
    Globals
*/

/* TODO: remove MAX_REACTIONS; use dynamic storage instead */
#define MAX_REACTIONS 100

fptr _setup, _initialize;
double* states;
double* dt_ptr;
int num_states;
int* neighbor_index;
double* neighbor_rates;
int* neighbor_list;
ReactionRate reactions[MAX_REACTIONS];
int num_reactions;
int _reaction_index;
int _num_species_per_location;
int _num_locations;
int _indices[MAX_REACTIONS];

//int _num_states_involved[MAX_REACTIONS];
//int _num_location[MAX_REACTIONS];
//int* _reaction_indices[MAX_REACTIONS];
//int* _change_states[MAX_REACTIONS];

void clear_rates(void) {
    num_reactions = 0;
}

void register_rate(ReactionRate f) {
    assert(num_reactions < MAX_REACTIONS);
    reactions[num_reactions++] = f;
}

void setup_solver(PyHocObject* my_states, PyHocObject* my_dt_ptr, int my_num_states, int* my_neighbor_index, int* my_neighbor_list, double* my_neighbor_rates) {
    states = my_states -> u.px_;
    num_states = my_num_states;
    dt_ptr = my_dt_ptr -> u.px_;
    neighbor_index = my_neighbor_index;
    neighbor_list = my_neighbor_list;
    neighbor_rates = my_neighbor_rates;
}

void set_reaction_indices(int reaction_id, int num_states_involved, int num_location int* change_states, int* indices) {
    _num_states_involved[reaction_id] = num_states_involved;
    _num_location[reaction_id] = num_location;
    /* TODO: copy? */
    _reaction_indices[reaction_id] = indices;
    _change_states[reaction_id] = change_states;
}

void set_reaction_indices(int reaction_index, int num_species_per_location, int num_locations, int* indices) {
    _reaction_index = reaction_index;
    _num_species_per_location = num_species_per_location;
    _num_locations = _num_locations;
    _indices = indices;
}

void _fadvance(void) {
    double* old_states;
    double* states_for_reaction;
    double* states_for_reaction_dx;
    int i, j, k, l, m, max_j, j_count;
    double current;
    double shift;
    double dt = *dt_ptr;
    double dx = FLT_EPSILON;
    ReactionRate current_reaction;
    MAT *jacobian;
    MAT *jacobian_copy = NULL;
    VEC *x;
    VEC *b;
    PERM *pivot;

    old_states = (double*) malloc(sizeof(double) * num_states);
    assert(old_states); 
    for (i = 0; i < num_states; i++) {
        old_states[i] = states[i];
    }

    for (i = 0; i < num_states; i++) {
        current = states[i];
        shift = 0;
        for (j = neighbor_index[i]; j < neighbor_index[i + 1]; j++) {

            shift += neighbor_rates[j] * (old_states[neighbor_list[j]] - current);
        }
        states[i] += dt * shift;
    }

    for (i = 0; i < _num_locations; ++i) {
        jacobian = m_get(num_reactions, _num_species_per_location);
        b = v_get(num_reactions);
        x = v_get(num_states);
        states_for_reaction = (double*) malloc(sizeof(double) * _num_species_per_location);
        states_for_reaction_dx = (double*) malloc(sizeof(double) * _num_species_per_location);

        for (j = 0; j < num_reactions; ++j){
            for (k = 0; k < _num_species_per_location; ++k){
                states_for_reaction[k] = old_states[k+_num_species_per_location];
                states_for_reaction_dx[k] = states_for_reaction[k];
            }
            for (l = 0; l < _num_species_per_location; ++l){
                states_for_reaction_dx[l] = states_for_reaction[l] + dx;
                double pd = (current_reaction(states_for_reaction_dx) - current_reaction(states_for_reaction))/dx;
                m_set_val(jacobian, i, l, pd);
                states_for_reaction_dx[l] -= dx;
            }
        }
    }
}

/*
void _fadvance(void) {
    double* old_states;
    double* states_for_reaction;
    double* states_for_reaction_dx;
    int i, j, k, l, m, max_j, j_count;
    double current;
    double shift;
    double dt = *dt_ptr;
    double dx = FLT_EPSILON;
    int num_states_involved;
    ReactionRate current_reaction;
    MAT *jacobian;
    MAT *jacobian_copy = NULL;
    VEC *x;
    VEC *b;
    PERM *pivot;

    jacobian = m_get(num_reactions, num_states);
    b = v_get(num_reactions);
    x = v_get(num_states);
    

    old_states = (double*) malloc(sizeof(double) * num_states);
    assert(old_states); 
    for (i = 0; i < num_states; i++) {
        old_states[i] = states[i];
    }

    for (i = 0; i < num_states; i++) {
        current = states[i];
        shift = 0;

        for (j = neighbor_index[i]; j < neighbor_index[i + 1]; j++) {

            shift += neighbor_rates[j] * (old_states[neighbor_list[j]] - current);
        }
        states[i] += dt * shift;
    }


    for (i = 0; i < num_reactions; i++) {

        num_states_involved = _num_states_involved[i];
        states_for_reaction = (double*) malloc(sizeof(double) * num_states_involved);
        states_for_reaction_dx = (double*) malloc(sizeof(double) * num_states_involved);
        assert(states_for_reaction);
        assert(states_for_reaction_dx);

        max_j = _num_location[i] * num_states_involved;
        current_reaction = reactions[i];

        for (j_count = j = 0; j < max_j; j += num_states_involved, j_count++) {

            for (k = 0; k < num_states_involved; k++) {
                states_for_reaction[k] = old_states[_reaction_indices[i][j + k]];
                states_for_reaction_dx[k] = states_for_reaction[k];
            }

            for (l = 0; l < num_states_involved; l++) {

                states_for_reaction_dx[l] = states_for_reaction[l] + dx;
                double pd = (current_reaction(states_for_reaction_dx) - current_reaction(states_for_reaction))/dx;
                m_set_val(jacobian, i, l, pd);
                states_for_reaction_dx[l] -= dx;
            }

            v_set_val(b, i, dt * current_reaction(states_for_reaction));
        }
        free(states_for_reaction);
        free(states_for_reaction_dx);
    }


    
    int ct, total_locations = 0;
    for (ct = 0; ct < num_reactions; ++ct){
        total_locations += _num_location[ct];
    }

    for (j_count = 0; j_count < total_locations; j_count++) {
        for (i = 0; i < num_reactions; i++){
            num_states_involved = _num_states_involved[i];
            states_for_reaction = (double*) malloc(sizeof(double) * num_states_involved);
            states_for_reaction_dx = (double*) malloc(sizeof(double) * num_states_involved);
            current_reaction = reactions[i];

            for (k = 0; k < num_states_involved; l++){
                states_for_reaction[k] = old_states[_reaction_indices[i][k]];
                states_for_reaction_dx[k] = states_for_reaction[k];
            }

            for (l = 0; l < num_states_involved){
                states_for_reaction_dx[l] = states_for_reaction[l] + dx;
                double pd = (current_reaction(states_for_reaction_dx) - current_reaction(states_for_reaction))/dx;
                m_set_val(jacobian, i, l, pd);
                states_for_reaction_dx[l] -= dx;
            }
        }

        jacobian_copy = m_copy(jacobian, jacobian_copy);

        for (i = 0; i < num_reactions; ++i){
            for (j = 0; j < num_states; ++j){
                if (i == j) m_set_val(jacobian_copy, i, j, 1-m_get_val(jacobian, i, j)*dt);
                else m_set_val(jacobian_copy, i, j, -m_get_val(jacobian, i, j)*dt);
            }
        }

        pivot = px_get(jacobian_copy->m);
        LUfactor(jacobian_copy, pivot);
        LUsolve(jacobian_copy, pivot, b, x);

        printf("Changed state: \n");
        for (m = 0; m < num_states; ++m){
            states[m] += v_get_val(x, m);
            printf("%f\t", states[m]);
        }
        free(old_states);
        free(jacobian_copy);
        free(jacobian);
        free(b);
        free(pivot);
    }

    we loop through each location. 
    each iteration, we go through each reaction, each one being an rxd.rate and loop thru each state
    to generate a partial derivative which we cumulatively add to ones in the same location
    to generate a change for each location. As such, we would have a jacobian for each location.
    
    
    jacobian_copy = m_copy(jacobian, jacobian_copy);

    for (i = 0; i < num_reactions; ++i){
        for (j = 0; j < num_states; ++j){
            if (i == j) m_set_val(jacobian_copy, i, j, 1-m_get_val(jacobian, i, j)*dt);
            else m_set_val(jacobian_copy, i, j, -m_get_val(jacobian, i, j)*dt);
        }
    }

    printf("\nJacobian \n");
    m_output(jacobian);
    printf("\n\n");

    printf("I - dtJ\n");
    m_output(jacobian_copy);
    printf("\n\n");

    pivot = px_get(jacobian_copy->m);
    LUfactor(jacobian_copy, pivot);
    LUsolve(jacobian_copy, pivot, b, x);

    printf("Y Change: \n");
    v_output(x);
    printf("\n");

    printf("Changed state: \n");
    for (m = 0; m < num_states; ++m){
        states[m] += v_get_val(x, m);
        printf("%f\t", states[m]);
    }

    printf("\n");

    free(old_states);
    free(jacobian_copy);
    free(jacobian);
    free(b);
    free(pivot);
}
*/

int rxd_nonvint_block(int method, int size, double* p1, double* p2, int thread_id) {
    switch (method) {
        case 0:
            _setup();
            break;
        case 1:
            _initialize();
            break;
        case 2:
            /* compute outward current to be subtracted from rhs */
            break;
        case 3:
            /* conductance to be added to d */
            break;
        case 4:
            /* fixed step solve */
            _fadvance();
            break;
        case 5:
            /* ode_count */
            break;
        case 6:
            /* ode_reinit(y) */
            break;
        case 7:
            /* ode_fun(t, y, ydot); from t and y determine ydot */
            break;
        case 8:
            /* ode_solve(dt, t, b, y); solve mx=b replace b with x */
        case 9:
            /* ode_jacobian(dt, t, ypred, fpred); optionally prepare jacobian for fast ode_solve */
            break;
        case 10:
            /* ode_abs_tol(y_abs_tolerance); fill with cvode.atol() * scalefactor */
            break;
        default:
            printf("Unknown rxd_nonvint_block call: %d\n", method);
            break;
    }
	/* printf("method=%d, size=%d, thread_id=%d\n", method, size, thread_id);	 */
}

void set_setup(fptr setup_fn) {
	_setup = setup_fn;
}

void set_initialize(fptr initialize_fn) {
	_initialize = initialize_fn;
}

