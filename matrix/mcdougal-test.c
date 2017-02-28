#include "matrix.h"
#include <stdlib.h>

/*
// works but unnecessary
void set_val(MAT* m, int i, int j, double v) {
    m->me[i][j] = v;
}
*/

int main(void) {
    MAT* a;
    MAT* acopy = NULL;
    VEC* x;
    VEC* b;
    PERM* pivot;
    printf("Hello from mcdougal-test.c\n");

    /* create a 2x2 matrix */
    a = m_get(2, 2);

    /* store some data */
    m_set_val(a, 0, 0, 2);
    m_set_val(a, 0, 1, 2);
    m_set_val(a, 1, 0, 2);
    m_set_val(a, 1, 1, 2);

    /* create a 2-vector RHS */
    b = v_get(2);
    v_set_val(b, 0, 2);
    v_set_val(b, 1, 2);

    /* and storage for x */
    x = v_get(2);

    /* solve the matrix problem; work on a copy of the matrix because the solve is destructive */
    acopy = m_copy(a, acopy);
    pivot = px_get(acopy->m);
    LUfactor(acopy, pivot);
    LUsolve(acopy, pivot, b, x);
    free(acopy);
    acopy = NULL;
    free(pivot);
    pivot = NULL;

    /* print the matrix */
    printf("A:\n");
    m_output(a);
    printf("\n\n");

    /* print x */
    printf("x:\n");
    v_output(x);
    printf("\n\n");

    /* print b (aka the rhs) */
    printf("b:\n");
    v_output(b);
    printf("\n\n");

    printf("Check:\n");
    printf("%g   %g\n", m_get_val(a, 0, 0) * v_get_val(x, 0) + m_get_val(a, 0, 1) * v_get_val(x, 1), m_get_val(a, 1, 0) * v_get_val(x, 0) + m_get_val(a, 1, 1) * v_get_val(x, 1));

    free(a);
    free(x);
    free(b);

    return EXIT_SUCCESS;
}