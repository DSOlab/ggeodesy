#! /usr/bin/python
#-*- coding: utf-8 -*-

from __future__ import print_function
import sys, os, re
import argparse
from shutil import copyfile
from shutil import move 

def printer(*msg): print(*msg, file = sys.stderr)

def setCppStandard(makefile, std, overwrite=False):
    mk_out = open(makefile+'.scratch', 'w') if overwrite else sys.stdout
    std_pattern = '.*-std=([^ ]*) *.*'
    with open(makefile, 'r') as mk_in:
        for line in mk_in.readlines():
            regex_match = re.match(std_pattern, line)
            if not regex_match: print(line.strip(), file=mk_out)
            else: print(line.strip().replace(regex_match.group(1), std), file=mk_out)
    if overwrite:
        mk_out.close()
        move(makefile+'.scratch', makefile)
    return
"""
def setMakefileSymbol(makefile, symbol_dct={}, symbol_lst=[]):
    ## example symbol_dct={'std':'c++17'}
    symbol_match = {}
    for key in symbol_dct: symbol_match[key] = False
    for key in symbol_lst: symbol_match[key] = False
    if not os.path.isfile(makefile):
        msg=("[ERROR] Failed to locate directory path {:}".format(
          args.project_dir))
        raise RuntimeError(msg)
    with open(makefile, 'r') as mk_in:
        line = mk_in.readline()
        ## read untill the first on-comment occurance of [CM]XXFLAGS
        while not re.match(' *[^#].*[CM]XXFLAGS.*', line):
            print(line.strip())
            line = mk_in.readline()
        ## read untill first blank line (or whitespace only line)
        flag_lines = []
        while not re.match(' *$', line):
            line_matched = False
            for symbol in symbol_dct:
                regex_match = re.match('.*{:} *= *([^ ]*) *.*'.format(symbol), line)
                if regex_match:
                    # print(line.strip().replace(regex_match.group(1), symbol_dct[symbol]))
                    flag_lines.append(line.replace(regex_match.group(1), symbol_dct[symbol]))
                    symbol_match[symbol] = line_matched = True
            for symbol in symbol_lst:
                regex_match = re.match('.*({:}).*'.format(symbol), line)
                if regex_match:
                    symbol_match[symbol] = line_matched = True
            if not line_matched: flag_lines.append(line)
            line = mk_in.readline()
        for key in symbol_dct:
            if not symbol_match[key]:
                flag_lines.append('{:}={:} \\\n'.format(key, symbol_dct[key]))
        for key in symbol_lst:
            if not symbol_match[key]:
                flag_lines.append('{:} \\\n'.format(key))
        for line in flag_lines: print(line.strip())
        ## read rest of file
        while line:
            print(line.strip())
            line = mk_in.readline()
        """

class myFormatter(
    argparse.ArgumentDefaultsHelpFormatter,
    argparse.RawTextHelpFormatter):
  pass
parser = argparse.ArgumentParser(
    formatter_class=myFormatter,
    description=('Setup ggeodesy Project for autotools'),
    epilog=('National Technical University of Athens\n'
      'Dionysos Satellite Observatory\n'
      'Send bug reports to:\n'
      'Xanthos Papanikolaou, xanthos@mail.ntua.gr'
      'Dimitris Anastasiou,danast@mail.ntua.gr\n'
      'November, 2020'))

parser.add_argument('-c', '--compile-mode',
    default='debug',
    metavar='COMPILE_MODE',
    dest='compile_mode',
    required=False,
    choices=['debug', 'production'],
    help=('Choose compile mode (i.e. debug or production)'))

parser.add_argument('-s', '--cpp-standard',
    default='17',
    metavar='CPP_STANDARD',
    dest='cpp_std',
    required=False,
    choices=['17', 'c++17', '14', 'c++14', '11', 'c++11'],
    help=('Choose C++ standard for building (i.e. argument for \"-std=\")'))

parser.add_argument('-d', '--project-dir',
    default=os.path.abspath("."),
    metavar='PROJECT_DIR',
    dest='project_dir',
    required=False,
    help=('Path of the project.'))

parser.add_argument('--include-boost',
    dest='include_boost',
    action='store_true',
    help=('Include the folder /boost and relevent sources in the build.'))

##  parse cmd
args  = parser.parse_args()

## folders where we need a makefile.am
mk_folders = ['src', 'test']
if args.include_boost: mk_folders.append('boost')

## check that the Project Directory is ok
if not os.path.isdir(args.project_dir):
    printer("[ERROR] Failed to locate directory path {:}".format(
      args.project_dir))
    sys.exit(1)

## check that we can locate all needed folders and Makefiles
for dir in [ os.path.join(args.project_dir, d) for d in mk_folders ]:
    if not os.path.isdir(dir):
        printer("[ERROR] Failed to locate directory path {:}".format(dir))
    for mk_file in ['Makefile.am.debug', 'Makefile.am.production']:
        if not os.path.isfile(os.path.join(args.project_dir, dir, mk_file)):
            printer("[ERROR] Failed to locate file {:}".format(os.path.join(
              args.project_dir, dir, mk_file)))

## cool, let's copy the needed makefiles
mk_target = 'Makefile.am.production' if args.compile_mode == 'production' else 'Makefile.am.debug'
for dir in [ os.path.join(args.project_dir, d) for d in mk_folders ]:
    copyfile(os.path.join(dir, mk_target), os.path.join(dir, 'Makefile.am'))

## nice, lets now create the top-directory Makefile.am
with open(os.path.join(args.project_dir, 'Makefile.am'), 'w') as mk_out:
    print('SUBDIRS = {:}'.format(' '.join(mk_folders)), file=mk_out)

## set the cpp standard
try:
  int(args.cpp_std)
  cpp_std = 'c++' + args.cpp_std
except:
  cpp_std = args.cpp_std
for dir in [ os.path.join(args.project_dir, d) for d in mk_folders ]:
    setCppStandard(os.path.join(dir, 'Makefile.am'), cpp_std, True)

## all done
print('All done! Makefile.am created in all needed folders. Now run the '
  'following to build the project:\nautoreconf -if\n./configure\nmake\n'
  'and you should be good to go!')
