#!/usr/bin/python

import sys, re, os

if len(sys.argv) > 1:

    fileName = sys.argv[1]

    with open(fileName, 'r') as file:
        data = file.read()

    # Equation blocks: $$ ... $$ -> \f[ ... \f]

    data = re.sub(r'(?<!\\)\$(?<!\\)\$(.+?)(?<!\\)\$(?<!\\)\$', r'\\f[\1\\f]', data, flags=re.S)

    # Inline equations: $ ... $ -> \f$ ... \f$

    data = re.sub(r'(?<!\\)(?<!\\f)\$(.+?)(?<!\\)(?<!\\f)\$', r'\\f$\1\\f$', data)

    # Relative links to relative links w.r.t. path of fileName

    path = re.sub(r'^(.*[^/])$', r'\1/', os.path.dirname(fileName))

    data = re.sub(r'(?<!\!)\[(.+?)\]\(((?!http|\/).+?)\)', r'[\1]('+path+r'\2)', data)

    print(data)
