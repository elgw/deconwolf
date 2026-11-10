#!/bin/env/python
import sys

outfile = 'pydocstrings.h'
infile = 'pydocstrings.txt'

class fdoc:

    def __init__(self):
        self.name : str = ''
        self.header : str = ''
        self.body : list[str] = []
        self.state : int = 0

    def print(self, fid):
        while (len(self.body) > 0) and (self.body[-1] == ''): self.body.pop()

        print(f'PyDoc_STRVAR({self.name}__doc__,', file=fid)
        print(rf'"{self.header}\n"', file=fid)
        print(r'"--\n\n"', file=fid)
        for line in self.body:
            tabbed = False
            if (len(line) > 0) and (line[0] == '\t'):
                tabbed = True
            if (len(line) > 1) and (line[0:2] == '  '):
                tabbed = True
            if tabbed:
                print(rf'"   {line.lstrip()}\n"', file=fid)
            else:
                print(rf'"{line}\n"', file=fid)
        print(");", file=fid)

    def add(self, line):
        line = line.rstrip()
        if self.state == 0:
            self.name = line
            self.state += 1
            return
        if self.state == 1:
            self.header = line
            self.state += 1
            return
        if self.state == 2:
            if len(line) == 0 and len(self.body) == 0:
                return
            self.body.append(line)
            return

if __name__ == '__main__':
    with open('pydocstrings.txt', 'r') as fid:
        funs = []
        fun = None
        for line in fid:
            if line[0] == '[' and line[-2] == ']':
                if not fun is None:
                    funs.append(fun)
                fun = fdoc()
                fun.add(line[1:-2])
            else:
                fun.add(line)
        funs.append(fun)


    print(f"Found {len(funs)} documented function")
    print(f'Writing to {outfile}')
    with open(outfile, 'w') as fid:
        print(f"// this file was generated from {infile} using genpydocstrings.py", file=fid)
        for fun in funs:
            print(fun.name)
            fun.print(fid)
