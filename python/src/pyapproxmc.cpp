/*
 ApproxMC

 Copyright (c) 2019-2020, Mate Soos and Kuldeep S. Meel. All rights reserved
 Copyright (c) 2009-2018, Mate Soos. All rights reserved.
 Copyright (c) 2015, Supratik Chakraborty, Daniel J. Fremont,
 Kuldeep S. Meel, Sanjit A. Seshia, Moshe Y. Vardi
 Copyright (c) 2014, Supratik Chakraborty, Kuldeep S. Meel, Moshe Y. Vardi

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in
 all copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 THE SOFTWARE.
 */

#include <Python.h>
#include "../cryptominisat/src/cryptominisat.h"
#include "../../src/approxmc.h"
#include "../arjun/src/arjun.h"

#include <limits>
#include <vector>
#include <set>

#define MODULE_NAME "pyapproxmc"
#define MODULE_DOC "ApproxMC approximate model counter."

typedef struct {
    PyObject_HEAD
    ApproxMC::AppMC* appmc = NULL;
    ArjunNS::Arjun* arjun = NULL;
    std::vector<CMSat::Lit> tmp_cl_lits;
    bool count_called = false;

    int verbosity;
    uint32_t seed;
    double epsilon;
    double delta;
} Counter;

static const char counter_create_docstring[] = \
"Counter(verbosity=0, seed=1, epsilon=0.8, delta=0.2)\n\
Create Counter object.\n\
\n\
:param verbosity: Verbosity level: 0: nothing printed; 15: very verbose.\n\
:param seed: Random seed\n\
:param epsilon: epsilon parameter as per PAC guarantees\n\
:param delta: delta parameter as per PAC guarantees";

/********** Internal Functions **********/

/* Helper functions */

static void setup_counter(Counter *self, PyObject *args, PyObject *kwds)
{
    static char const* kwlist[] = {"verbosity", "seed", "epsilon", "delta", NULL};

    // All parameters have the same default as the command line defaults
    // except for verbosity which is 0 by default.
    self->verbosity = 0;
    self->seed = 1;
    self->epsilon = 0.8;
    self->delta = 0.2;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|iIdd", const_cast<char**>(kwlist),
        &self->verbosity, &self->seed, &self->epsilon, &self->delta))
    {
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
        PyErr_SetString(PyExc_ValueError, "delta must be greater or equal to 0, and less than 1");
        return;
    }

    self->appmc = new ApproxMC::AppMC;
    self->appmc->set_verbosity(self->verbosity);
    self->appmc->set_seed(self->seed);
    self->appmc->set_epsilon(self->epsilon);
    self->appmc->set_delta(self->delta);

    self->arjun = new ArjunNS::Arjun;
    self->arjun->set_seed(self->seed);
    self->arjun->set_verbosity(self->verbosity);

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

static int parse_clause(Counter *self, PyObject *clause, std::vector<CMSat::Lit>& lits, bool allow_more_vars = true)
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

    if (!lits.empty() && max_var >= (long int)self->arjun->nVars()) {
        if (allow_more_vars)
            self->arjun->new_vars(max_var-(long int)self->arjun->nVars()+1);
        else {
            PyErr_SetString(PyExc_ValueError,
                    "ERROR: Sampling vars contain variables that are not in the original clauses!");
            return 0;
        }
    }

    Py_DECREF(iterator);
    if (PyErr_Occurred()) return 0;
    return 1;
}

static int _add_clause(Counter *self, PyObject *clause)
{
    self->tmp_cl_lits.clear();
    if (!parse_clause(self, clause, self->tmp_cl_lits)) {
        return 0;
    }
    self->arjun->add_clause(self->tmp_cl_lits);

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

template <typename T>
static int _add_clauses_from_array(Counter *self, const size_t array_length, const T *array)
{
    if (array_length == 0) {
        return 1;
    }
    if (array[array_length - 1] != 0) {
        PyErr_SetString(PyExc_ValueError, "last clause not terminated by zero");
        return 0;
    }
    size_t k = 0;
    long val = 0;
    std::vector<CMSat::Lit>& lits = self->tmp_cl_lits;
    for (val = (long) array[k]; k < array_length; val = (long) array[++k]) {
        lits.clear();
        long int max_var = 0;
        for (; k < array_length && val != 0; val = (long) array[++k]) {
            long var;
            bool sign;
            if (val > std::numeric_limits<int>::max()/2
                || val < std::numeric_limits<int>::min()/2
            ) {
                PyErr_Format(PyExc_ValueError, "integer %ld is too small or too large", val);
                return 0;
            }

            sign = (val < 0);
            var = std::abs(val) - 1;
            max_var = std::max(var, max_var);

            lits.push_back(CMSat::Lit(var, sign));
        }
        if (!lits.empty()) {
            if (max_var >= (long int)self->arjun->nVars()) {
                self->arjun->new_vars(max_var-(long int)self->arjun->nVars()+1);
            }
            self->arjun->add_clause(lits);
        }
    }
    return 1;
}

static int _add_clauses_from_buffer(Counter *self, Py_buffer *view)
{
    if (view->ndim != 1) {
        PyErr_Format(PyExc_ValueError, "invalid clause array: expected 1-D array, got %d-D", view->ndim);
        return 0;
    }
    if (strcmp(view->format, "i") != 0 && strcmp(view->format, "l") != 0 && strcmp(view->format, "q") != 0) {
        PyErr_Format(PyExc_ValueError, "invalid clause array: invalid format '%s'", view->format);
        return 0;
    }

    void * array_address = view->buf;
    size_t itemsize = view->itemsize;
    size_t array_length = view->len / itemsize;

    if (itemsize == sizeof(int)) {
        return _add_clauses_from_array(self, array_length, (const int *) array_address);
    }
    if (itemsize == sizeof(long)) {
        return _add_clauses_from_array(self, array_length, (const long *) array_address);
    }
    if (itemsize == sizeof(long long)) {
        return _add_clauses_from_array(self, array_length, (const long long *) array_address);
    }
    PyErr_Format(PyExc_ValueError, "invalid clause array: invalid itemsize '%ld'", itemsize);
    return 0;
}

PyDoc_STRVAR(add_clauses_doc,
"add_clauses(clauses)\n\
Add iterable of clauses to the solver.\n\
\n\
:param clauses: List of clauses. Each clause contains literals (ints)\n\
    Alternatively, this can be a flat array.array or other contiguous\n\
    buffer (format 'i', 'l', or 'q') of zero separated and terminated\n\
    clauses of literals (ints).\n\
:type clauses: <list> or <array.array>\n\
:return: None\n\
:rtype: <None>"
);

static PyObject* add_clauses(Counter *self, PyObject *args, PyObject *kwds)
{
    static char const* kwlist[] = {"clauses", NULL};
    PyObject *clauses;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", const_cast<char**>(kwlist), &clauses)) {
        return NULL;
    }

    if (PyObject_CheckBuffer(clauses)) {
        Py_buffer view;
        memset(&view, 0, sizeof(view));
        if (PyObject_GetBuffer(clauses, &view, PyBUF_CONTIG_RO | PyBUF_FORMAT) != 0) {
            return NULL;
        }

        int ret = _add_clauses_from_buffer(self, &view);
        PyBuffer_Release(&view);

        if (ret == 0 || PyErr_Occurred()) {
            return 0;
        }
        Py_INCREF(Py_None);
        return Py_None;
    }

    PyObject *iterator = PyObject_GetIter(clauses);
    if (iterator == NULL) {
        PyErr_SetString(PyExc_TypeError, "iterable object expected");
        return NULL;
    }

    PyObject *clause;
    while ((clause = PyIter_Next(iterator)) != NULL) {
        _add_clause(self, clause);
        /* release reference when done */
        Py_DECREF(clause);
    }

    /* release reference when done */
    Py_DECREF(iterator);
    if (PyErr_Occurred()) {
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}


static void get_cnf_from_arjun(Counter* self)
{
    const uint32_t orig_num_vars = self->arjun->get_orig_num_vars();
    self->appmc->new_vars(orig_num_vars);
    self->arjun->start_getting_small_clauses(
        std::numeric_limits<uint32_t>::max(),
        std::numeric_limits<uint32_t>::max(),
        false);
    std::vector<CMSat::Lit> clause;

    bool ret = true;
    while (ret) {
        ret = self->arjun->get_next_small_clause(clause);
        if (!ret) {
            break;
        }

        bool ok = true;
        for(auto l: clause) {
            if (l.var() >= orig_num_vars) {
                ok = false;
                break;
            }
        }

        if (ok) self->appmc->add_clause(clause);
    }
    self->arjun->end_getting_small_clauses();
}

static void transfer_unit_clauses_from_arjun(Counter* self)
{
    std::vector<CMSat::Lit> cl(1);
    auto units = self->arjun->get_zero_assigned_lits();
    for(const auto& unit: units) {
        if (unit.var() < self->appmc->nVars()) {
            cl[0] = unit;
            self->appmc->add_clause(cl);
        }
    }
}

static uint32_t set_up_sampling_set(Counter* self, const std::vector<uint32_t>& sampling_vars)
{
    uint32_t orig_sampling_set_size;
    orig_sampling_set_size = self->arjun->set_starting_sampling_set(sampling_vars);
    return orig_sampling_set_size;
}

/* count function */

PyDoc_STRVAR(count_doc,
"count(projection)\n\
Approximately count the number of solutions to the formula. It can only be called ONCE.\n\
\n\
:param projection: the projection over which to count the solutions over\n\
\n\
:return: A tuple. The first part of the tuple is the cell solution count and\n\
    the second part is the hash count."
);

static PyObject* count(Counter *self, PyObject *args, PyObject *kwds)
{
    if (self->count_called) {
        PyErr_SetString(PyExc_ValueError, "ERROR: Counter.count() may only be called once!");
        return NULL;
    } else {
        self->count_called = true;
    }
    static char const* kwlist[] = {"projection", NULL};
    PyObject *py_sampling_vars = NULL;
    if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", const_cast<char**>(kwlist), &py_sampling_vars)) {
        return NULL;
    }

    // Get sampling vars from user
    std::vector<uint32_t> sampling_vars;
    if (py_sampling_vars == NULL) {
        for(uint32_t i = 0; i < self->arjun->nVars(); i++) sampling_vars.push_back(i);
    } else {
        std::vector<CMSat::Lit> sampling_lits;
        if (parse_clause(self, py_sampling_vars, sampling_lits, false) == 0 ) {
            return NULL;
        }
        for(const auto& l: sampling_lits) {
            if (l.var() > self->arjun->nVars()) {
                PyErr_SetString(PyExc_ValueError,
                        "ERROR: Sampling vars contain variables that are not in the original clauses!");
                return NULL;
            }
            sampling_vars.push_back(l.var());
        }
    }

   //print_orig_sampling_vars(sampling_vars, self->arjun);
   uint32_t orig_sampling_set_size = set_up_sampling_set(self, sampling_vars);
   sampling_vars = self->arjun->get_indep_set();
   std::vector<uint32_t> empty_occ_sampl_vars = self->arjun->get_empty_occ_sampl_vars();
   //print_final_indep_set(sampling_vars , orig_sampling_set_size, empty_occ_sampl_vars);

    std::set<uint32_t> sampl_vars_set;
    sampl_vars_set.insert(sampling_vars.begin(), sampling_vars.end());
    for(auto const& v: empty_occ_sampl_vars) {
        assert(sampl_vars_set.find(v) != sampl_vars_set.end()); // this is guaranteed by arjun
        sampl_vars_set.erase(v);
    }
    const size_t offset_count_by_2_pow = empty_occ_sampl_vars.size();
    sampling_vars.clear();
    sampling_vars.insert(sampling_vars.end(), sampl_vars_set.begin(), sampl_vars_set.end());

    // Now do ApproxMC
    get_cnf_from_arjun(self);
    transfer_unit_clauses_from_arjun(self);
    ApproxMC::SolCount sol_count;
    if (!sampling_vars.empty()) {
        self->appmc->set_projection_set(sampling_vars);
        sol_count = self->appmc->count();
    } else {
        bool ret = self->appmc->find_one_solution();
        sol_count.hashCount = 0;
        if (ret) sol_count.cellSolCount = 1;
        else sol_count.cellSolCount = 0;
    }

    // Fill return value
    PyObject *result = PyTuple_New((Py_ssize_t) 2);
    if (result == NULL) {
        PyErr_SetString(PyExc_SystemError, "failed to create a tuple");
        return NULL;
    }
    PyTuple_SET_ITEM(result, 0, PyLong_FromLong((long)sol_count.cellSolCount));
    PyTuple_SET_ITEM(result, 1, PyLong_FromLong((long)sol_count.hashCount+offset_count_by_2_pow));
    return result;
}

/********** Python Bindings **********/
static PyMethodDef Counter_methods[] = {
    {"count",     (PyCFunction) count,       METH_VARARGS | METH_KEYWORDS, count_doc},
    {"add_clause",(PyCFunction) add_clause,  METH_VARARGS | METH_KEYWORDS, add_clause_doc},
    {"add_clauses", (PyCFunction) add_clauses,  METH_VARARGS | METH_KEYWORDS, add_clauses_doc},
    {NULL, NULL}  // Sentinel
};

static void Counter_dealloc(Counter* self)
{
    delete self->appmc;
    delete self->arjun;
    Py_TYPE(self)->tp_free ((PyObject*) self);
}

static int Counter_init(Counter *self, PyObject *args, PyObject *kwds)
{
    if (self->appmc != NULL) delete self->appmc;
    if (self->arjun != NULL) delete self->arjun;

    setup_counter(self, args, kwds);

    if (!self->appmc) return -1;

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

    // Add the version string
    // they're using.
#if defined(_MSC_VER)
#else
    if (PyModule_AddStringConstant(m, "__version__", APPMC_FULL_VERSION) == -1) {
        Py_DECREF(m);
        return NULL;
    }
    if (PyModule_AddStringConstant(m, "VERSION", APPMC_FULL_VERSION) == -1) {
        Py_DECREF(m);
        return NULL;
    }
#endif

    // Add the Counter type
    Py_INCREF(&pyapproxmc_CounterType);
    if (PyModule_AddObject(m, "Counter", (PyObject *)&pyapproxmc_CounterType)) {
        Py_DECREF(m);
        return NULL;
    }

    return m;
}
