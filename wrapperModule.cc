#include <Python.h>
#include "test.h"

static char module_docstring[] =
    "This module runs C-based Robospect";
static char test_docstring[] =
    "Run a test script.";

static PyObject *test_wrapperModule(PyObject *self, PyObject *args);

static PyMethodDef module_methods[] = {
    {"test", test_wrapperModule, METH_VARARGS, test_docstring},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC init_wrapperModule(void)
{
    PyObject *m = Py_InitModule3("_wrapperModule", module_methods, module_docstring);
    if (m == NULL)
        return;
}
