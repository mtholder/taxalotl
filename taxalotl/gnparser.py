#!/usr/env python
"""Adapter for calling gnparser

Adapted based on instructions on
    https://medium.com/learning-the-go-programming-language/calling-go-functions-from-other-languages-4c7d8bcc69bf
"""
import ctypes
import json
import os

_notfound_message = """LIBGNPARSER_FILEPATH must be in your env (or 
the gnparse_shared_object_filepath argument must be supplied).

See https://gitlab.com/gogna/gnparser#command-line for instructions for building the 
library.
"""

_gnparse_lib = None

def get_gnparse_lib(gnparse_shared_object_filepath=None):
    global _gnparse_lib
    if _gnparse_lib is None:
        if not gnparse_shared_object_filepath:
            gnparse_shared_object_filepath = os.environ.get('LIBGNPARSER_FILEPATH')
            if not gnparse_shared_object_filepath:
                raise RuntimeError(_notfound_message)
        if not os.path.isfile(gnparse_shared_object_filepath):
            m = 'LIBGNPARSER_FILEPATH "{}" is not a file '.format(gnparse_shared_object_filepath)
            raise RuntimeError(m)
        _gnparse_lib = ctypes.cdll.LoadLibrary(gnparse_shared_object_filepath)
        _gnparse_lib.ParseToString.argtypes = [ctypes.c_char_p, ctypes.c_char_p]
        _gnparse_lib.ParseToString.restype = ctypes.c_char_p
    return _gnparse_lib
_compact_arg = ctypes.c_char_p("compact".encode('utf-8'))


def parse_name_to_dict(name, gnparse_shared_object_filepath=None):
    lib = get_gnparse_lib(gnparse_shared_object_filepath=gnparse_shared_object_filepath)
    cn = ctypes.c_char_p(name.encode('utf-8'))
    x = lib.ParseToString(cn, _compact_arg).decode('utf-8')
    return json.loads(x)


def main():
    import sys
    for arg in sys.argv[1:]:
        x = parse_name_to_dict(arg)
        print(json.dumps(x, ensure_ascii=True, indent=2))

if __name__ == '__main__':
    main()
