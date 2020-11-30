#! /usr/bin/python
#-*- coding: utf-8 -*-

from __future__ import print_function
import sys, os
import argparse
from shutil import copyfile

def printer(*msg): print(*msg, file = sys.stderr)

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

