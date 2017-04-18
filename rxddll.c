
    void reaction(double* states, double* rhs) {
        rhs[0] = -states[1];
        rhs[1] = states[0];
    }


    