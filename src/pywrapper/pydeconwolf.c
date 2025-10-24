/*    Copyright (C) 2025 Erik L. G. Wernersson
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

/* A thin wrapper around some of the functionality in deconwolf
 * TODO:
 * store docstrings in normal text files and use `xxd -i` to generate header files
 */


#include <Python.h>

#define NPY_NO_DEPRECATED_API NPY_2_0_API_VERSION
#include <numpy/ndarrayobject.h>
#include "numpy/npy_3kcompat.h"

// TODO: for sure we have something like this in the dw codebase
#ifdef __GNUC__
#  define UNUSED(x) UNUSED_ ## x __attribute__((__unused__))
#else
#  define UNUSED(x) UNUSED_ ## x
#endif

#include <deconwolf/fim_tiff.h>
#include <deconwolf/dw_version.h>
#include <deconwolf/dw_bwpsf.h>
#include <deconwolf/dw.h>
#include <deconwolf/fim.h>

typedef int64_t i64;
typedef float f32;

#include "docstrings.h"

// TODO: use standard error types, this clutters the documentation and does
// not add anything useful
static PyObject *pydw_error;

static int
check_matrix_argument(PyObject * _obj, int mindim, int maxdim)
{
    assert(mindim <= maxdim);

    if(!PyArray_Check(_obj))
    {
        printf("Not an numpy array");
        return EXIT_FAILURE;
    }

    PyArrayObject * obj = (PyArrayObject*) _obj;

    if(0)
    {
        int flags = PyArray_FLAGS(obj);
        printf("   NPY_ARRAY_C_CONTIGUOUS: %d\n", (flags & NPY_ARRAY_C_CONTIGUOUS) > 0);
        printf("   NPY_ARRAY_F_CONTIGUOUS: %d\n", (flags & NPY_ARRAY_F_CONTIGUOUS) > 0);
        printf("        NPY_ARRAY_OWNDATA: %d\n", (flags & NPY_ARRAY_OWNDATA) > 0);
        printf("        NPY_ARRAY_ALIGNED: %d\n", (flags & NPY_ARRAY_ALIGNED) > 0);
        printf("      NPY_ARRAY_WRITEABLE: %d\n", (flags & NPY_ARRAY_WRITEABLE) > 0);
        printf("NPY_ARRAY_WRITEBACKIFCOPY: %d\n", (flags & NPY_ARRAY_WRITEBACKIFCOPY) > 0);
    }

    int ndim = PyArray_NDIM(obj);
    int type = PyArray_TYPE(obj);

    // https://numpy.org/doc/stable/reference/c-api/dtype.html
    if(type != NPY_FLOAT32)
    {
        printf("Incorrect data type, only NPY_FLOAT32 supported at the moment\n");
        return EXIT_FAILURE;
    }

    if(ndim < mindim || ndim > maxdim)
    {
        if(mindim == maxdim)
        {
            printf("Wrong number of dimensions, should be %d\n", mindim);
        } else {
            printf("Wrong number of dimensions, should be %d to %d!\n", mindim, maxdim);
        }
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}


static PyObject *
x_imread(PyObject * UNUSED(self), PyObject *args, PyObject *keywds)
{

    char * filename = NULL;
    int verbose = 0;

    static char *kwlist[] = {"filename", "verbose", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, keywds, "s|i", kwlist,
                                     &filename, &verbose))
    {
        return NULL;
    }
    fim_tiff_init();

    i64 M = 0;
    i64 N = 0;
    i64 P = 0;
    float * V = fim_tiff_read(filename, NULL, &M, &N, &P, verbose);

    if(V == NULL)
    {
        PyErr_SetString(pydw_error,
                        "Can't open that file\n");
        return NULL;
    }

    const npy_intp outsize[3] = {P, N, M};
    PyObject * result_obj =
        PyArray_SimpleNewFromData(3, outsize, NPY_FLOAT32, V);
    return result_obj;
}

static PyObject *
x_bw(PyObject * UNUSED(self), PyObject *args, PyObject *keywds)
{
    bw_conf * conf = bw_conf_new();

    static char *kwlist[] = {"NA", "ni",
                             "emission", "dx", "dz",
                             "nplanes", "size", "verbose", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, keywds, "fffff|iii", kwlist,
                                     &conf->NA, &conf->ni,
                                     &conf->lambda, &conf->resLateral, &conf->resAxial,
                                     // Optional
                                     &conf->P, &conf->M, &conf->verbose))
    {
        PyErr_SetString(pydw_error,
                        "Something weird with the arguments\n");
        return NULL;
    }

    conf->N = conf->M;

    if(bw_conf_validate(conf))
    {
        PyErr_SetString(pydw_error,
                        "The arguments does not make sense, please see the documentation.\n");
        return NULL;
    }
    if(conf->verbose > 1)
    {
        bw_conf_printf(stdout, conf);
    }

    BW(conf);

    if(conf->V == NULL)
    {
        bw_conf_free(&conf);
        PyErr_SetString(pydw_error,
                        "Something went wrong\n");
        return NULL;
    }

    const npy_intp outsize[3] = {conf->P, conf->N, conf->M};
    PyObject * result_obj =
        PyArray_SimpleNewFromData(3, outsize, NPY_FLOAT32, conf->V);
    bw_conf_free(&conf);
    return result_obj;
}


static PyObject *
x_imwrite(PyObject * UNUSED(self), PyObject *args, PyObject *keywds)
{
    const char * filename = NULL;
    int verbose = 0;
    PyObject * _array;
    static char *kwlist[] = {"array", "filename", "verbose", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, keywds, "Os|i", kwlist,
                                     &_array, &filename, &verbose))
    {
        PyErr_SetString(pydw_error,
                        "Error parsing the arguments\n");
        return NULL;
    }
    if(check_matrix_argument(_array, 2, 3))
    {
        PyErr_SetString(pydw_error,
                        "Not a valid image array\n");
        return NULL;
    }

    fim_tiff_init();

    PyArrayObject * array = (PyArrayObject*) _array;

    int M = PyArray_DIM(array, 1); // Py_ssize_t
    int N = PyArray_DIM(array, 0);
    int P = 1;

    if(PyArray_NDIM(array) == 3)
    {
        M = PyArray_DIM(array, 2);
        N = PyArray_DIM(array, 1);
        P = PyArray_DIM(array, 0);
    }

    fim_tiff_write_float(filename, (f32*) PyArray_DATA(array), NULL, M, N, P);
    Py_RETURN_NONE;
}

static PyObject *
x_deconvolve(PyObject * UNUSED(self), PyObject *args, PyObject *keywds)
{

    dw_opts * config = dw_opts_new();
    config->write_result_to_buffer = 1;
    config->iter_type = DW_ITER_FIXED; // disable autostop

    PyObject * _image;
    PyObject * _psf;
    int use_gpu = 0;
    static char *kwlist[] = {"image", "psf", "iterations", "bq", "verbose", "threads", "gpu", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, keywds, "OOi|iiip", kwlist,
                                     &_image, &_psf, &config->nIter,
                                     &config->borderQuality,
                                     &config->verbosity,
                                     &config->nThreads_FFT,
                                     &use_gpu))
    {
        dw_opts_free(&config);
        PyErr_SetString(pydw_error,
                        "Error parsing the arguments\n");
        return NULL;
    }

    if(use_gpu)
    {
        dw_opts_enable_gpu(config);
    }

    config->nThreads_OMP = config->nThreads_FFT;
    config->log = stdout;
    config->leave_log_open = 1;

    if(PyUnicode_Check(_image))
    {
        const char* s = PyUnicode_AsUTF8(_image);
        printf("image: %s\n", s);
        config->imFile = strdup(s);
    } else {
        if(check_matrix_argument(_image, 3, 3) == 0)
        {
            PyArrayObject * image = (PyArrayObject*) _image;
            i64 M = PyArray_DIM(image, 2);
            i64 N = PyArray_DIM(image, 1);
            i64 P = PyArray_DIM(image, 0);
            config->input_image_size.M = M;
            config->input_image_size.N = N;
            config->input_image_size.P = P;
            config->input_image = fim_malloc(M*N*P*sizeof(float));
            memcpy(config->input_image, PyArray_DATA(image), M*N*P*sizeof(float));
            printf("Got a %ld x %ld x %ld image\n", M, N, P);
        } else
        {
            dw_opts_free(&config);
            PyErr_SetString(pydw_error,
                            "Not a valid image array\n");
            return NULL;
        }

    }

    if(PyUnicode_Check(_psf))
    {
        const char* s = PyUnicode_AsUTF8(_psf);
        printf("psf: %s\n", s);
        config->psfFile = strdup(s);
    } else {
        if(check_matrix_argument(_psf, 3, 3) == 0)
        {
            PyArrayObject * psf = (PyArrayObject*) _psf;
            i64 M = PyArray_DIM(psf, 2);
            i64 N = PyArray_DIM(psf, 1);
            i64 P = PyArray_DIM(psf, 0);
            config->input_psf_size.M = M;
            config->input_psf_size.N = N;
            config->input_psf_size.P = P;
            config->input_psf = fim_malloc(M*N*P*sizeof(float));
            memcpy(config->input_psf, PyArray_DATA(psf), M*N*P*sizeof(float));
            printf("Got a %ld x %ld x %ld image\n", M, N, P);
        } else
        {
            dw_opts_free(&config);
            PyErr_SetString(pydw_error,
                            "Not a valid image array\n");
            return NULL;
        }
    }

    if(dw_opts_validate_and_init(config))
    {
        PyErr_SetString(pydw_error,
                        "Something failed during dw_opts_validate_and_init\n");
        return NULL;
    }

    dw_run(config);

    if(config->result == NULL)
    {
        PyErr_SetString(pydw_error,
                        "Something failed here, I did not get any output back from deconwolf\n");
        return NULL;
    }

    const npy_intp outsize[3] = {config->result_size.P, config->result_size.N, config->result_size.M};
    PyObject * result_obj =
        PyArray_SimpleNewFromData(3, outsize, NPY_FLOAT32, config->result);

    dw_opts_free(&config);
    return result_obj;
}


// The 2nd element is allowed to be a PyCFunctionWithKeywords
static PyMethodDef pydeconwolf_methods[] = {
    {"imread",      (PyCFunction) x_imread,    METH_VARARGS | METH_KEYWORDS, imread__doc__},
    {"imwrite",     (PyCFunction) x_imwrite,   METH_VARARGS | METH_KEYWORDS, imwrite__doc__},
    {"gen_psf_bw",  (PyCFunction) x_bw,        METH_VARARGS | METH_KEYWORDS, bw__doc__},
    {"deconvolve",  (PyCFunction) x_deconvolve, METH_VARARGS | METH_KEYWORDS, deconvolve__doc__},
    {NULL, NULL, 0, NULL},
};

/* Module definition structure */
static struct PyModuleDef
pydeconwolf_module_info = {
    PyModuleDef_HEAD_INIT,
    "pydeconwolf",
    pydeconwolf__doc__,
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    pydeconwolf_methods,
    NULL,
    NULL,
    NULL,
    NULL
};


PyMODINIT_FUNC
PyInit_pydeconwolf(void)
{
    PyObject* module = PyModule_Create(&pydeconwolf_module_info);
    if (module == NULL)
    {
        fprintf(stderr, "Unable to initialize the module!");
        fprintf(stderr, "%s/%s():%d\n", __FILE__, __FUNCTION__, __LINE__);
        return NULL;
    }

    PyModule_AddIntConstant(module, "VERSION_MAJOR", atoi(DW_VERSION_MAJOR));
    PyModule_AddIntConstant(module, "VERSION_MINOR", atoi(DW_VERSION_MINOR));
    PyModule_AddIntConstant(module, "VERSION_PATCH", atoi(DW_VERSION_PATCH));

    import_array();

    pydw_error = PyExc_RuntimeError;

    return module;
}
