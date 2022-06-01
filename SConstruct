from __future__ import print_function
import os, sys, glob

## Prefic for install(ed) files
prefix="/usr/local"
if not os.path.isdir(prefix):
    print('[ERROR] Cannot find \'prefix\' directory, aka {:}; aborting'.format(prefix), file=sys.stderr)
    sys.exit(1)

## Library version
lib_version="0.1.0"
## Library name
lib_name="geodesy"
## Include dir (following prefix) if any
inc_dir="geodesy"
## the rootdir of the project
root_dir=os.path.abspath(os.getcwd())

## get number of CPUs and use for parallel builds
num_cpu = int(os.environ.get('NUM_CPU', 2))
SetOption('num_jobs', num_cpu)
print("running with -j %s" % GetOption('num_jobs'))

## Source files (for lib)
lib_src_files = glob.glob(r"src/*.cpp")

## Headers (for lib)
hdr_src_files = glob.glob(r"src/*.hpp")

## Environments ...
denv = Environment(CXXFLAGS='-std=c++17 -g -pg -Wall -Wextra -Werror -pedantic -W -Wshadow -Winline -Wdisabled-optimization -DDEBUG')
## g++ complains about failing to inline functions if we use the '-Winline' here ... 
penv = Environment(CXXFLAGS='-std=c++17 -Wall -Wextra -Werror -pedantic -W -Wshadow -O2 -march=native')

## Command line arguments ...
debug = ARGUMENTS.get('debug', 0)
boostg = ARGUMENTS.get('boost', 0)
eigen = ARGUMENTS.get('eigen', 0)

## Construct the build enviroment
env = denv.Clone() if int(debug) else penv.Clone()

## (shared) library ...
vlib = env.SharedLibrary(source=lib_src_files, target=lib_name, CPPPATH=['.'], SHLIBVERSION=lib_version)

## Build ....
#env.Alias(target='install', source=env.Install(dir=os.path.join(prefix, 'lib'), source=vlib))
env.Alias(target='install', source=env.Install(dir=os.path.join(prefix, 'include', inc_dir), source=hdr_src_files))
env.Alias(target='install', source=env.InstallVersionedLib(dir=os.path.join(prefix, 'lib'), source=vlib))

if eigen:
    math_lib = ''
else:
    math_lib = 'matvec'

## Tests ...
tests_sources = glob.glob(r"test/*.cpp")
env.Append(RPATH=root_dir)
for tsource in tests_sources:
  ttarget = tsource.replace('_', '-').replace('.cpp', '.out')
  env.Program(target=ttarget, source=tsource, CPPPATH='src/', LIBS=vlib+['datetime', math_lib], LIBPATH='.')

## Boost test executables
if boostg:
  print('>> note that we\'ll be building boost executables ...')
  boost_sources = ['boost/test_geodtest.cpp', 'boost/test_geodtime.cpp']
  for bsource in boost_sources:
    btarget = bsource.replace('_', '-').replace('.cpp', '.out')
    env.Program(target=btarget, source=bsource, CPPPATH='src/',
                LIBS=vlib+['datetime', math_lib], LIBPATH='.')
