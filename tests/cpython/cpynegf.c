#include <Python.h>
#include "libnegf.h"

static PyObject* cpynegf_gethandlersize(PyObject* self, PyObject* args)
{
  int hsize;
  negf_gethandlersize_(&hsize);
  return Py_BuildValue("i",hsize);
}

/*
* Bind Python function names to our C functions
*/
static PyMethodDef CpynegfMethods[] = {
    {"gethandlersize", cpynegf_gethandlersize, METH_VARARGS},
        {NULL, NULL}
};


PyMODINIT_FUNC
initcpynegf(void)
{
    (void) Py_InitModule("cpynegf", CpynegfMethods);
}
