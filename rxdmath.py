# TODO: don't misue str() like this... make my own thing (e.g. _semi_compile)

class _Arithmeticed:
    def __str__(self):
        return self._str
    def __abs__(self, obj):
        result = _Arithmeticed()
        result._str = 'fabs(%s)' % obj
        return result
    def __neg__(self):
        result = _Arithmeticed()
        result._str = '(-%s)' % self
        return result
    def __add__(self, obj):
        result = _BinaryOp()
        result.left = self
        result.right = obj
        result.op = '+'
        return result
    def __sub__(self, obj):
        result = _BinaryOp()
        result.left = self
        result.right = obj
        result.op = '-'
        return result
    def __mul__(self, obj):
        result = _BinaryOp()
        result.left = self
        result.right = obj
        result.op = '*'
        return result
    def __div__(self, obj):
        result = _BinaryOp()
        result.left = self
        result.right = obj
        result.op = '/'
        return result
    def __radd__(self, obj):
        return _ensure_arithmeticed(obj) + self
    def __rmul__(self, obj):
        return _ensure_arithmeticed(obj) * self
    def __rdiv__(self, obj):
        return _ensure_arithmeticed(obj) / self
    def __rsub__(self, obj):
        return _ensure_arithmeticed(obj) - self
    def _to_c(self, function_name, args):
        return 'double %s(%s) {return %s;}' % (function_name, args, self)
    def _compile_given_args(self, args):
        import uuid
        import os
        import ctypes
        filename = 'rxd-' + str(uuid.uuid1())
        with open(filename + '.c', 'w') as f:
            f.write(self._to_c('myfunc', args))
        os.system('gcc -I/usr/include/python2.7 -lpython2.7 -shared -o %s.so -fPIC %s.c' % (filename, filename))
        dll = ctypes.cdll['./%s.so' % filename]  
        self._filename = filename
        myfunc = dll.myfunc
        myfunc.argtypes = [ctypes.c_double] * (1 + args.count(','))
        myfunc.restype = ctypes.c_double
        return dll.myfunc
    def _compile(self):
        return self._compile_given_args('double* _species')

class _BinaryOp(_Arithmeticed):
    def __str__(self):
        return '(%s %s %s)' % (self.left, self.op, self.right)


def _ensure_arithmeticed(obj):
    if isinstance(obj, _Arithmeticed):
        return obj
    else:
        result = _Arithmeticed()
        result._str = str(obj)
        return result
