// #include <python.h>
#include "python.h"
#include <windows.h>
#include <stdlib.h>
// #include "longobject.c"

// ctrl + shift + B : run tasks.json


////////////
// source from : https://stackoverflow.com/questions/21888140/de-bruijn-algorithm-binary-digit-count-64bits-c-sharp
// used generator from http://chessprogramming.wikispaces.com/De+Bruijn+Sequence+Generator
Py_ssize_t DeBruijnMSB64table[] = 
{
    0 , 47, 1 , 56, 48, 27, 2 , 60,
    57, 49, 41, 37, 28, 16, 3 , 61,
    54, 58, 35, 52, 50, 42, 21, 44,
    38, 32, 29, 23, 17, 11, 4 , 62,
    46, 55, 26, 59, 40, 36, 15, 53,
    34, 51, 20, 43, 31, 22, 10, 45,
    25, 39, 14, 33, 19, 30, 9 , 24,
    13, 18, 8 , 12, 7 , 6 , 5 , 63,
};
// the cyclc number has to be in the last 16th of all possible values
// any beyond the 62914560th(0x03C0_0000) should work for this purpose
unsigned long long DeBruijnMSB64multi = 0x03F79D71B4CB0A89uL; // the last one
Py_ssize_t get_bit_size(Py_ssize_t value)
{
    Py_ssize_t val = (unsigned long long) value;
    val |= val >> 1;
    val |= val >> 2;
    val |= val >> 4;
    val |= val >> 8;
    val |= val >> 16;
    val |= val >> 32;
    return DeBruijnMSB64table[val * DeBruijnMSB64multi >> 58];
}

int PyObject_Print(PyObject* obj) {
    PyObject* repr = PyObject_Repr(obj);
    if (repr == NULL) {
        printf("Error calling repr() on the object!\n");
        return;
    }
    PyObject* str = PyUnicode_AsUTF8String(repr);
    if (str == NULL) {
        printf("Error converting repr() result to string!\n");
        Py_DECREF(repr);
        return;
    }
    const char* c_str = PyBytes_AsString(str);
    printf("%s\n", c_str);
    Py_DECREF(repr);
    Py_DECREF(str);
}

Py_ssize_t get_size_power_of_2(Py_ssize_t sizeu, Py_ssize_t sizev){
    printf("size %d, %d\n",sizeu,sizev);
    Py_ssize_t thres = sizeu;
    if (sizeu < sizev) thres = sizev;
    assert(thres<=(1<<30)); // arbitrary limitation
    Py_ssize_t cnt = 1;
    Py_ssize_t size = 1;
    while(size < thres){
        size = size<<1;
        cnt++;
    }
    return cnt+1;   
}

PyObject* mod_by_b(PyObject *U,PyObject *W, Py_ssize_t b, Py_ssize_t log_b){
    PyObject* Y_mod_b = PyList_New(b);
    PyObject* U_with_pad = PyLong_FromLong(0);
    PyObject* W_with_pad = PyLong_FromLong(0);
    PyObject* bit_mask = PyNumber_Subtract(PyLong_FromLong(b),PyLong_FromLong(1));
    PyObject* shift_amt = PyLong_FromLong(3*log_b);
    PyObject* Y_with_pad;
    for(Py_ssize_t i=b-1;i>=0;i--){
        PyObject* U_with_pad_shifted = PyNumber_Lshift(U_with_pad, shift_amt);
        PyObject* Ui = PyList_GET_ITEM(U,i);        
        PyObject* Ui_bit_masked = PyNumber_And(Ui,bit_mask);
        U_with_pad = PyNumber_Or(U_with_pad_shifted, Ui_bit_masked);

        PyObject* W_with_pad_shifted = PyNumber_Lshift(W_with_pad, shift_amt);
        PyObject* Wi = PyList_GET_ITEM(W,i);        
        PyObject* Wi_bit_masked = PyNumber_And(Wi,bit_mask);
        W_with_pad = PyNumber_Or(W_with_pad_shifted, Wi_bit_masked);
    }

    Y_with_pad = PyNumber_Multiply(U_with_pad, W_with_pad); 
    // python multiply is a built-in karasuba_mul

    for(Py_ssize_t i=0;i<b;i++){
        PyList_SET_ITEM(Y_mod_b,i, PyNumber_And(Y_with_pad, bit_mask));
        Y_with_pad = PyNumber_Rshift(Y_with_pad, shift_amt);
    }
    return Y_mod_b;
}

PyObject* mod_bit_recursive(PyObject *a, Py_ssize_t l){
    printf("----mod_bit_recursive\n");
    PyObject* ret = PyLong_FromLong(0);
    PyObject* one = PyLong_FromLong(1);
    PyObject* zero = PyLong_FromLong(0);
    PyObject* shift_size = PyLong_FromLongLong(2*l);
    PyObject* thres = PyNumber_Lshift(one, shift_size);
    thres = PyNumber_Add(thres,one);    
    // PyNumber_InPlaceAdd(thres,one); not working as  expected
    if (PyObject_RichCompareBool(a,zero,Py_EQ)){
        return ret;
    }
    else if (PyObject_RichCompareBool(a,zero,Py_LT)){
        a = PyNumber_Add(a,thres);
        // PyNumber_InPlaceAdd(a,thres);
    }
    PyObject* bit_mask = PyNumber_Subtract(PyNumber_Lshift(one, shift_size), one);
    PyObject* a_back;
    
    boolean is_positive = TRUE;
    while (PyObject_RichCompareBool(a,zero,Py_GT)){
        a_back = PyNumber_And(a,bit_mask);
        a = PyNumber_Rshift(a, shift_size);
        // PyNumber_InPlaceRshift(a, shift_size);
        if (is_positive){
            ret = PyNumber_Add(ret,a_back);
            // PyNumber_InPlaceAdd(ret,a_back);
            if (PyObject_RichCompareBool(ret,thres,Py_GE)){
                ret = PyNumber_Subtract(ret,thres);
                // PyNumber_InPlaceSubtract(ret,thres);
            }
        }
        else{
            ret = PyNumber_Subtract(ret,a_back);
            // PyNumber_InPlaceSubtract(ret,a_back);
            if (PyObject_RichCompareBool(ret,zero,Py_LT)){
                ret = PyNumber_Add(ret,thres);
                // PyNumber_InPlaceAdd(ret,thres);
            }
        }
        is_positive = !is_positive;
    }
    PyObject_Print(ret,stdout,1);
    printf("  ret\n");
    return ret;
    
}

void bit_reversal(PyObject *x){
    Py_ssize_t rev;
    Py_ssize_t bit;
    PyObject* tmp;
    printf("%u\n",PyList_Size(x));
    for(Py_ssize_t i=0;i<PyList_Size(x);i++){
        Py_ssize_t i_backup = i;
        bit = i;
        rev = 0;
        // printf("%u\n",get_bit_size(PyList_Size(x))-1);
        for(Py_ssize_t _=0;_<get_bit_size(PyList_Size(x))-1;_++){
            rev = (rev<<1) + (i&1);
            i = i >> 1;
        }
        if (bit<rev){
            tmp = PyList_GET_ITEM(x, bit);
            PyList_SET_ITEM(x, bit, PyList_GET_ITEM(x, rev));
            PyList_SET_ITEM(x, rev, tmp);
        }
        i = i_backup;
    }
    Py_DECREF(tmp);
    printf("bit reversal end--\n");
    PyObject_Print(x);
}

void fft_in_place(PyObject *U, PyObject *omegas,Py_ssize_t l, Py_ssize_t b){
    printf("------------Starting fft\n");
    Py_ssize_t iterations = get_bit_size(PyList_Size(U))-1;
    bit_reversal(U);
    Py_ssize_t M = 2;
    Py_ssize_t m = 1;
    Py_ssize_t k;
    PyObject *even,*odd,*k_pos_num;
    for(Py_ssize_t _=0;_<iterations;_++){
        for(Py_ssize_t i=0;i<PyList_Size(U);i+=M){
            Py_ssize_t g = 0;
            printf("g\n");
            // PyObject_Print(PyNumber_Multiply(PyList_GET_ITEM(U, 1), PyList_GET_ITEM(omegas,g)));
            for(Py_ssize_t j=0;j<(M>>1);j++){
                k = i + j + (M >> 1);
                printf("k : %d\n",k);
                even = PyList_GET_ITEM(U, i+j);
                k_pos_num = PyList_GET_ITEM(U, k);
                printf("befoer odd\n");
                PyObject_Print(even);
                printf("k_pos_num : %d\n",Py_REFCNT(k_pos_num));
                printf("k_pos_num : %d\n",Py_REFCNT(k_pos_num));
                odd = mod_bit_recursive(PyNumber_Multiply(k_pos_num, PyList_GET_ITEM(omegas,g)), l);
                printf("k_pos_num : %d\n",k_pos_num->ob_refcnt);
                // need to change PyNumber_Multiply function to SSA later
                // PyObject_Print(even);
                printf("i+j %d\n",i+j);
                printf("%d\n",k);
                PyList_SET_ITEM(U, i+j, mod_bit_recursive(PyNumber_Add(even,odd),l));
                PyList_SET_ITEM(U,   k, mod_bit_recursive(PyNumber_Subtract(even,odd),l));
                g = g + (PyList_Size(U) >> m);
                g = g&(b-1);
            }
        }    
        m++;
        M <<= 1;
    }
    printf("fft resulet\n");
    PyObject_Print(U);
    PyObject_Print(U);
    PyObject_Print(U);
    PyObject_Print(U);
    
}

PyObject* pointwise_multiplication(PyObject *U, PyObject *W, Py_ssize_t l){
    PyObject *Y = PyList_New(PyList_Size(U));
    for(Py_ssize_t i=0;i<PyList_Size(Y);i++){
        PyObject *Yi = mod_bit_recursive(PyNumber_Multiply(PyList_GET_ITEM(U,i),PyList_GET_ITEM(W,i)),1);
        // TODO: need SSA instead PyNumber_Multiply
        PyList_SET_ITEM(Y,i,Yi);
    }
    return Y;
}

void divide_by_b(PyObject *Y, Py_ssize_t l, Py_ssize_t b){
    Py_ssize_t shift_amt = (l<<2)-(get_bit_size(b)-1);
    PyObject *divide_amt = _PyLong_Lshift(PyLong_FromLong(1) ,shift_amt);
    for(Py_ssize_t i=0;i<PyList_Size(Y);i++){
        PyObject *Yi = PyList_GET_ITEM(Y,i);
        PyList_SET_ITEM(Y,i, mod_bit_recursive(PyNumber_Multiply(Yi, divide_amt),l));
    }
}

void chinese_remainder(PyObject *Y, PyObject *Y_mod_b, Py_ssize_t b, Py_ssize_t l){
    printf("Chinese remainder\n");
    for(Py_ssize_t i=0;i<PyList_Size(Y);i++){
        PyObject *Yi = PyList_GET_ITEM(Y,i);
        PyObject *Y_mod_bi = PyList_GET_ITEM(Y_mod_b,i);
        PyObject *delta = PyNumber_Subtract(Y_mod_bi, Yi);
        PyObject_Print(delta);
        delta = PyNumber_And(delta, PyLong_FromLong(b-1));
        PyObject_Print(delta);
        Yi = PyNumber_Add(Yi, delta);
        Yi = PyNumber_Add(Yi, _PyLong_Lshift(delta,2*l));
        PyList_SET_ITEM(Y, i, Yi);
    }
    printf("Chinese remainder\n");
}

void print_PyLongObject(PyLongObject *obj, boolean endLine) {
    Py_ssize_t i, size = Py_ABS(Py_SIZE(obj));
    printf("(");
    for (i = 0; i < size; i++) {
        printf("%d ", obj->ob_digit[i]);
    }
    printf(")");
    if (endLine){
        printf("\n");
    }
}

void print_PyList_PyLongObject(PyObject* obj){
    printf("[");
    PyLongObject *ele;
    for (Py_ssize_t i = 0; i < Py_SIZE(obj); i++) {
        ele = PyList_GET_ITEM(obj, i);
        print_PyLongObject(ele,FALSE);
        printf(" ");
    }
    Py_DECREF(ele);
    printf("]\n");
}



PyObject* PyLong_Copy(PyObject* obj) {
    long value = PyLong_AsLong(obj);
    if (value == -1 && PyErr_Occurred()) {
        return NULL;
    }
    return PyLong_FromLong(value);
}

void pylongobject_shift(PyLongObject *num,  PyObject *shift_amount){
    PyObject *shifted = Py_TYPE(num)->tp_as_number->nb_lshift((PyObject *) num, shift_amount);
}

void cal_omega_list(PyObject *omegas, PyObject *inverse_omegas, Py_ssize_t k_rem, Py_ssize_t b, Py_ssize_t l){
    // Elements of omega and inverse_omegas are smaller than 2**4l
    // omega = 2**(4l/b) = 1<<(4*(1<<(log_l - log_b)))
    PyObject *log_omega = PyLong_FromLong(0);
    PyObject *log_inverse_omega;
    PyObject *omega = PyLong_FromLong(1);
    PyObject *log_omega_add = PyLong_FromLong(PyLong_SHIFT*4/(k_rem+1));

    PyObject *log_omega_inverse_sub = PyLong_FromLong(4*l);
    PyObject* inverse_omega = PyLong_FromLong(1);

    PyList_SET_ITEM(omegas, 0, omega);
    PyList_SET_ITEM(inverse_omegas, 0, inverse_omega);

    for (Py_ssize_t i = 1; i < b; i++) {
        log_omega = PyNumber_Add(log_omega,log_omega_add); // inplace
        omega = PyNumber_Lshift(omega, log_omega); // inplace
        // omegas[i] = omegas[i-1]<<(log_omega)
        
        PyList_SET_ITEM(omegas, i, mod_bit_recursive(omega,l));

        log_inverse_omega = PyNumber_Subtract(log_omega_inverse_sub,log_omega);
        inverse_omega = PyNumber_Lshift(PyLong_FromLong(1), log_inverse_omega);
        // need mod_bit_op here
        PyList_SET_ITEM(inverse_omegas, i, mod_bit_recursive(inverse_omega,l));
    }    
}

void numToArray(PyObject *U, PyLongObject *u,Py_ssize_t k_half,Py_ssize_t b){
    Py_ssize_t size_u = Py_ABS(Py_SIZE(u));
    printf("----------------------------\n");
    printf("size %d\n",size_u);
    for (Py_ssize_t i = 0; i < size_u; i++) {
        printf("i=%d, v=%d\n",i, u->ob_digit[i]);
    }
    digit *pu = u->ob_digit;
    digit *p_last = u->ob_digit + size_u;
    for (Py_ssize_t i = 0; i < b; i++) {
        PyLongObject *num_piece = _PyLong_New(1<<k_half);
        // PyLong_Copy(num_piece, PyLong_FromLong(0)); // initialize num_piece to zero, TODO : memset?
        memset(num_piece->ob_digit, 0, (1<<k_half) * sizeof(digit));
        print_PyLongObject(num_piece,TRUE);
        
        // num_piece->ob_base.ob_type = &PyLong_Type;
        for (Py_ssize_t j = 0; j < (1<<k_half); j++) {
            if(pu >= p_last){
                break;
            }
            printf("pu %d\n",*pu);
            num_piece->ob_digit[j] = *pu++;
        }
        PyList_SET_ITEM(U, i, (PyObject*)num_piece);
    }
}

PyObject* array_to_num(PyObject *Y, Py_ssize_t l){
    printf("--array to num--\n");
    PyObject *y = PyLong_FromLong(0);
    for(Py_ssize_t i=PyList_GET_SIZE(Y)-1;i>=0;i--){
        PyObject *Yi = PyList_GET_ITEM(Y,i);
        y = PyNumber_Add(Yi,_PyLong_Lshift(y,l));
        PyObject_Print(y);
    }
}

int mai2n(){
    Py_Initialize();
    ///////////////////////////////////////////
    PyObject *tmp1 = PyLong_FromLong((1<<29)*3);
    PyObject *tmp2 = PyLong_FromLong((1<<20)*41+1264);
    if (PyObject_RichCompareBool(tmp1,tmp2,Py_GT)){
        return 0;
    }
    printf("wfe\n");
    PyObject *u;
    u = PyNumber_Multiply(tmp1,tmp2);
    PyObject_Print(u,stdout,Py_PRINT_RAW); // Segmentation fault here

    const char* uu = PyUnicode_AsUTF8(PyObject_Repr(u));
    printf("%s\n", uu); // it works
    Py_Finalize();
    return 0;
}

static PyObject* SSA(PyObject* self, PyObject* args){
// int main(){
    // putenv("PYTHONHOME=");
    // putenv("PYTHONIOENCODING=utf-8");

    
    if (!Py_IsInitialized()) {
    fprintf(stderr, "Failed to initialize Python interpreter\n");
    return 1;
    }
    
    int num1, num2, sts;
    PyLongObject *u;
    PyLongObject *v;
    
    // PyLongObject *U, *V;

    // 
    u = PyTuple_GET_ITEM(args, 0);
    v = PyTuple_GET_ITEM(args, 1);
    
    ///////////////////////////////////////////

    // v->ob_base.ob_size = 2; // 
    // u->ob_base.ob_size = 5;
    // v = (PyLongObject *)tmp1;
    // u = (PyLongObject *)tmp2;
    // if (!PyArg_ParseTuple(args,"ii",&num1,&num2)){ // ii means integer, integer
    //     return NULL;
    // }
    // u = PyLong_FromLong(2^33+1);
    // v = PyLong_FromLong(2^33+3);

}

static PyObject* _SSA(PyLongObject *u, PyLongObject *w){
    printf("---staring _SSA");
    Py_ssize_t usize = Py_ABS(Py_SIZE(u));
    Py_ssize_t vsize = Py_ABS(Py_SIZE(w));
    Py_ssize_t k = get_size_power_of_2(usize,vsize); // IGNORE size of ob_digit;
    Py_ssize_t k_half = k>>1;
    Py_ssize_t k_rem = k&1;
    Py_ssize_t log_l = k_half; 
    Py_ssize_t log_b = k_half+k_rem;
    Py_ssize_t b = 1<<(log_b); // b=2**[k/2]_ceiling
    Py_ssize_t l = (1<<(log_l))*PyLong_SHIFT; // b=2**[k/2]_ceiling  // restore size of ob_digit here

    printf("log_l : %d\n", log_l);
    printf("log_b : %d\n", log_b);

    PyObject* U = PyList_New(b);
    PyObject* W = PyList_New(b);

    numToArray(U, u,k_half,b);
    numToArray(W, w,k_half,b);
    PyObject_Print(U,stdout,1);
    PyObject_Print(W,stdout,1);
    printf("\n");

    // return 0;
    PyObject* Y_mod_b = mod_by_b(U,W,b,log_b);

    PyObject* omegas = PyList_New(b);
    // omega = 2**(4l/b)
    PyObject* inverse_omegas = PyList_New(b);

    cal_omega_list(omegas, inverse_omegas, k_rem, b, l);
    printf("omegas");
    PyObject_Print(omegas);
    printf("%d\n", W->ob_refcnt);
    fft_in_place(W, omegas, l, b);
    printf("%d\n", W->ob_refcnt);
    fft_in_place(U, omegas, l, b);    
    
    printf("------------Starting pointwise mulitplicatoi\n");
    PyObject* Y = pointwise_multiplication(U,W,l);
    fft_in_place(Y,inverse_omegas,l,b);
    divide_by_b(Y,l,b);
    
    chinese_remainder(Y,Y_mod_b,b,l);
    return array_to_num(Y,l);
}

int main(){
    printf("---main started---------\n");
    Py_Initialize();
    printf("---Py_Initialized---------\n");
    // PyObject *test1 = PyLong_FromLong(128);
    // PyObject_Print(test1,stdout,1);
    // PyObject *test2 = PyLong_FromLong(1);
    // PyNumber_InPlaceRshift(test1,test2);
    // PyObject_Print(test1,stdout,1);
    // return 0;
    PyLongObject *u;
    PyLongObject *v;
    PyObject *tmp1 = PyLong_FromLong((1<<29)+(1<<40));
    PyObject *tmp2 = PyLong_FromLong((1<<20));    
    v = _PyLong_New(2);
    u = _PyLong_New(5);
    memset(v->ob_digit, 1, 2 * sizeof(digit));
    memset(u->ob_digit, 2, 5 * sizeof(digit));
    PyObject_Print(PyNumber_Multiply((PyObject *)u,(PyObject *)v),stdout,0);
    PyObject *result = _SSA(u,v);
    PyObject_Print(result,stdout,0);
    PyObject_Print(PyNumber_Multiply((PyObject *)u,(PyObject *)v),stdout,0);
    return 0;
}

static PyObject* c_add(PyObject* self, PyObject* args){
    int num1, num2, sts;
    PyLongObject *a = NULL; 
    PyLongObject *b = NULL;
    
    a = PyTuple_GET_ITEM(args, 0);
    b = PyTuple_GET_ITEM(args, 1);
    // if (!PyArg_ParseTuple(args,"ii",&num1,&num2)){ // ii means integer, integer
    //     return NULL;
    // }
    Py_ssize_t asize = Py_ABS(Py_SIZE(a));
    Py_ssize_t bsize = Py_ABS(Py_SIZE(b));
    // printf(a->ob_base.ob_size);

    // if a=12+(14)*2**30 -> a->ob_digit[0]=12; a->ob_digit[1]=14
    printf("%I64d \n",a->ob_digit[0]);
    printf("%I64d \n",a->ob_digit[1]);
    printf("%I64d \n",asize);
    printf("%I64d \n",bsize);

    return PyLong_FromLong(1+3);
}

static PyObject* version(PyObject* self){
    return Py_BuildValue("s","Version 0.01");
}

static PyMethodDef Examples[]={
    // METH_VARARGS : This is the typical calling convention,
    // where the methods have the type PyCFunction. The function expects two PyObject* values.
    // The first one is the self object for methods; for module functions, it is the module object.
    // The second parameter (often called args) is a tuple object representing all arguments.
    // This parameter is typically processed using PyArg_ParseTuple() or PyArg_UnpackTuple().
    {"c_add", c_add, METH_VARARGS, "add two numbers"},
    {"version", (PyCFunction) version, METH_NOARGS, "returns the version of the module"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef Example = {
    PyModuleDef_HEAD_INIT,
    "Example",
    "add Module",
    -1, // global setting
    Examples
};

//INITIALIZER FUNCTION
PyMODINIT_FUNC PyInit_Example(void){
    return PyModule_Create(&Example);
}