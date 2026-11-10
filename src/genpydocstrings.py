#!/bin/env python
import sys
import argparse


class fdoc:

    def __init__(self):
        self.is_function = True
        self.name : str = ''
        self.header : str = ''
        self.body : list[str] = []
        # 0: write to name
        # 1: write to header
        # 2: write to body
        self.state : int = 0

    def valid(self):
        if self.state == 2:
            return True
        return False

    def print(self, fid):
        # remove trailing empty lines
        while (len(self.body) > 0) and (self.body[-1] == ''): self.body.pop()

        print(f'PyDoc_STRVAR({self.name}__doc__,', file=fid)
        if self.is_function:
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
        if (len(line) > 1) and (line[0:2] == '//'):
            return
        if self.state == 0:
            if (len(line) > 1) and (line[0] == '[') and (line[-1] == ']'):
                self.name = line[1:-1]
                self.state = +2
                self.is_function = False
            else:
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

def write_h_file(conf, funs, fid):
        print(f"// this file was generated from {conf.infile} using genpydocstrings.py", file=fid)
        print("", file=fid)
        for fun in funs:
            fun.print(fid)

def parse_txt_file(conf, fid):
    funs = []
    fun = fdoc()
    for line in fid:
        if (len(line) > 1) and (line[0] == '[') and (line[-2] == ']'):
            if fun.valid():
                funs.append(fun)
            fun = fdoc()
            fun.add(line[1:-2])
        else:
            fun.add(line)

    funs.append(fun)
    return funs

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="convert pydocstrings.txt into pydocstrings.h",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument('--infile', required=False,
                        help='specify input file', default='pydocstrings.txt');

    parser.add_argument('--outfile', required=False,
                        help='specify input file', default='pydocstrings.h');
    parser.add_argument('--stdout', required=False,
                        help='Write to standard output instead of any file',
                        action='store_true')

    conf = parser.parse_args()

    with open(conf.infile, 'r') as fid:
        funs = parse_txt_file(conf, fid)

    if conf.stdout:
        write_h_file(conf, funs, sys.stdout);
        sys.exit(0)

    print(f"Found {len(funs)} documented function")
    for fun in funs:
        print(f"{fun.name} {len(fun.header)} {len(fun.body)}")

    print(f'Writing to {conf.outfile}')
    with open(conf.outfile, 'w') as fid:
        write_h_file(conf, funs, fid)
