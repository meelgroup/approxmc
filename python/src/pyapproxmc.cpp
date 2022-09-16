// Python bindings for ApproxMC, heavily based on the Python bindings written for CryptoMiniSat.

#include <Python.h>
#include "../cryptominisat/src/cryptominisat.h"
#include "../../src/approxmc.h"

#include <limits>

#define MODULE_NAME "pyapproxmc"
#define MODULE_DOC "ApproxMC approximate model counter."

typedef struct {
    PyObject_HEAD
    ApproxMC::AppMC* appmc;
    std::vector<CMSat::Lit> tmp_cl_lits;

    int verbosity;
    uint32_t seed;
    double epsilon;
    double delta;
    std::vector<uint32_t> sampling_set;
} Counter;

static const char counter_create_docstring[] = \
"Counter(verbosity=0, seed=1, epsilon=0.8, delta=0.2, sampling_set=None)\n\
Create Counter object.\n\
\n\
:param verbosity: Verbosity level: 0: nothing printed; 15: very verbose.\n\
:param seed: Random seed\n\
:param epsilon: epsilon parameter as per PAC guarantees\n\
:param delta: delta parameter as per PAC guarantees\n\
:param sampling_set: (Optional) If provided, the number of solutions counted\n\
    is over the variables in sampling_set.";

/********** Internal Functions **********/

/* Helper functions */

static int parse_sampling_set(Counter *self, PyObject *sample_set_obj)
{
    PyObject *iterator = PyObject_GetIter(sample_set_obj);
    if (iterator == NULL) {
        PyErr_SetString(PyExc_TypeError, "iterable object expected");
        return 1;
    }

    PyObject *lit;
    while ((lit = PyIter_Next(iterator)) != NULL) {
        long val = PyLong_AsLong(lit);
        if (val == 0) {
            PyErr_SetString(PyExc_ValueError, "non-zero integer expected");
            return 1;
        }
        if (val > std::numeric_limits<int>::max()/2
            || val < std::numeric_limits<int>::min()/2
        ) {
            PyErr_Format(PyExc_ValueError, "integer %ld is too small or too large", val);
            return 1;
        }

        long var = val - 1;
        self->sampling_set.push_back(var);
        Py_DECREF(lit);
    }
    Py_DECREF(iterator);

    return 0;
}

static void setup_counter(Counter *self, PyObject *args, PyObject *kwds)
{
    static char const* kwlist[] = {"verbosity", "seed", "epsilon", "delta", "sampling_set", NULL};

    // All parameters have the same default as the command line defaults
    // except for verbosity which is 0 by default.
    self->verbosity = 0;
    self->seed = 1;
    self->epsilon = 0.8;
    self->delta = 0.2;

    PyObject* sample_set_obj = NULL;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|iIddO", const_cast<char**>(kwlist),
        &self->verbosity, &self->seed, &self->epsilon, &self->delta, &sample_set_obj))
    {
        return;
    }

    if (sample_set_obj != NULL && parse_sampling_set(self, sample_set_obj)) {
        return;
    }

    if (self->verbosity < 0) {
        PyErr_SetString(PyExc_ValueError, "verbosity must be at least 0");
        return;
    }
    if (self->epsilon <= 0) {
        PyErr_SetString(PyExc_ValueError, "epsilon must be greater than 0");
        return;
    }
    if (self->delta < 0 || self->delta >= 1) {
        PyErr_SetString(PyExc_ValueError, "delta must be greater than or equal to 0 and less than 1");
        return;
    }

    self->appmc = new ApproxMC::AppMC;
    self->appmc->set_verbosity(self->verbosity);
    self->appmc->set_seed(self->seed);
    self->appmc->set_epsilon(self->epsilon);
    self->appmc->set_delta(self->delta);
    self->appmc->set_projection_set(self->sampling_set);

    return;
}

static int convert_lit_to_sign_and_var(PyObject* lit, long& var, bool& sign)
{
    if (!PyLong_Check(lit))  {
        PyErr_SetString(PyExc_TypeError, "integer expected !");
        return 0;
    }

    long val = PyLong_AsLong(lit);
    if (val == 0) {
        PyErr_SetString(PyExc_ValueError, "non-zero integer expected");
        return 0;
    }
    if (val > std::numeric_limits<int>::max()/2
        || val < std::numeric_limits<int>::min()/2
    ) {
        PyErr_Format(PyExc_ValueError, "integer %ld is too small or too large", val);
        return 0;
    }

    sign = (val < 0);
    var = std::abs(val) - 1;

    return 1;
}

static int parse_clause(Counter *self, PyObject *clause, std::vector<CMSat::Lit>& lits)
{
    PyObject *iterator = PyObject_GetIter(clause);
    if (iterator == NULL) {
        PyErr_SetString(PyExc_TypeError, "iterable object expected");
        return 0;
    }

    PyObject *lit;
    long int max_var = 0;
    while ((lit = PyIter_Next(iterator)) != NULL) {
        long var;
        bool sign;
        int ret = convert_lit_to_sign_and_var(lit, var, sign);
        Py_DECREF(lit);
        if (!ret) {
            Py_DECREF(iterator);
            return 0;
        }
        max_var = std::max(var, max_var);

        lits.push_back(CMSat::Lit(var, sign));
    }

    if (!lits.empty() && max_var >= (long int)self->appmc->nVars()) {
        self->appmc->new_vars(max_var-(long int)self->appmc->nVars()+1);
    }

    Py_DECREF(iterator);
    if (PyErr_Occurred()) {
        return 0;
    }

    return 1;
}

static int _add_clause(Counter *self, PyObject *clause)
{
    self->tmp_cl_lits.clear();
    if (!parse_clause(self, clause, self->tmp_cl_lits)) {
        return 0;
    }
    self->appmc->add_clause(self->tmp_cl_lits);

    return 1;
}


/* add_clause function */

PyDoc_STRVAR(add_clause_doc,
"add_clause(clause)\n\
Add a clause to the solver.\n\
\n\
:param clause: An iterator contains literals (ints)");

static PyObject* add_clause(Counter *self, PyObject *args, PyObject *kwds)
{
    static char const* kwlist[] = {"clause", NULL};
    PyObject *clause;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", const_cast<char**>(kwlist), &clause)) {
        return NULL;
    }

    if (_add_clause(self, clause) == 0 ) {
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;

}

/* count function */

PyDoc_STRVAR(count_doc,
"count()\n\
Approximately count the number of solutions for the clauses that have been \n\
added with add_clause().\n\
\n\
:return: A tuple. The first part of the tuple is the cell solution count and\n\
    the second part is the hash count."
);

static PyObject* count(Counter *self, PyObject *args, PyObject *kwds)
{
    PyObject *result = PyTuple_New((Py_ssize_t) 2);
    if (result == NULL) {
        PyErr_SetString(PyExc_SystemError, "failed to create a tuple");
        return NULL;
    }

    auto res = self->appmc->count();

    PyTuple_SET_ITEM(result, 0, PyLong_FromLong((long)res.cellSolCount));
    PyTuple_SET_ITEM(result, 1, PyLong_FromLong((long)res.hashCount));

    return result;
}

/********** Python Bindings **********/
static PyMethodDef Counter_methods[] = {
    {"count",     (PyCFunction) count,       METH_VARARGS | METH_KEYWORDS, count_doc},
    {"add_clause",(PyCFunction) add_clause,  METH_VARARGS | METH_KEYWORDS, add_clause_doc},
    {NULL, NULL}  // Sentinel
};

static void Counter_dealloc(Counter* self)
{
    delete self->appmc;
    Py_TYPE(self)->tp_free ((PyObject*) self);
}

static int Counter_init(Counter *self, PyObject *args, PyObject *kwds)
{
    if (self->appmc != NULL) {
        delete self->appmc;
    }

    setup_counter(self, args, kwds);

    if (!self->appmc) {
        return -1;
    }

    return 0;
}

static PyTypeObject pyapproxmc_CounterType = 
{
    PyVarObject_HEAD_INIT(NULL, 0)  /*ob_size*/
    "pyapproxmc.Counter",           /*tp_name*/
    sizeof(Counter),                /*tp_basicsize*/
    0,                              /*tp_itemsize*/
    (destructor)Counter_dealloc,    /*tp_dealloc*/
    0,                              /*tp_print*/
    0,                              /*tp_getattr*/
    0,                              /*tp_setattr*/
    0,                              /*tp_compare*/
    0,                              /*tp_repr*/
    0,                              /*tp_as_number*/
    0,                              /*tp_as_sequence*/
    0,                              /*tp_as_mapping*/
    0,                              /*tp_hash */
    0,                              /*tp_call*/
    0,                              /*tp_str*/
    0,                              /*tp_getattro*/
    0,                              /*tp_setattro*/
    0,                              /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
    counter_create_docstring,       /* tp_doc */
    0,                              /* tp_traverse */
    0,                              /* tp_clear */
    0,                              /* tp_richcompare */
    0,                              /* tp_weaklistoffset */
    0,                              /* tp_iter */
    0,                              /* tp_iternext */
    Counter_methods,                /* tp_methods */
    0,                              /* tp_members */
    0,                              /* tp_getset */
    0,                              /* tp_base */
    0,                              /* tp_dict */
    0,                              /* tp_descr_get */
    0,                              /* tp_descr_set */
    0,                              /* tp_dictoffset */
    (initproc)Counter_init,         /* tp_init */
};

PyMODINIT_FUNC PyInit_pyapproxmc(void) 
{
    PyObject* m;

    pyapproxmc_CounterType.tp_new = PyType_GenericNew;
    if (PyType_Ready(&pyapproxmc_CounterType) < 0) {
        // Return NULL on Python3 and on Python2 with MODULE_INIT_FUNC macro
        // In pure Python2: return nothing.
        return NULL;
    }

    static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,  /* m_base */
        MODULE_NAME,            /* m_name */
        MODULE_DOC,             /* m_doc */
        -1,                     /* m_size */
        NULL,                   /* m_methods */
        NULL,                   /* m_reload */
        NULL,                   /* m_traverse */
        NULL,                   /* m_clear */
        NULL,                   /* m_free */
    };

    m = PyModule_Create(&moduledef);

    if (!m) {
        return NULL;
    }

    // Add the Counter type
    Py_INCREF(&pyapproxmc_CounterType);
    if (PyModule_AddObject(m, "Counter", (PyObject *)&pyapproxmc_CounterType)) {
        Py_DECREF(m);
        return NULL;
    }

    return m;
}
