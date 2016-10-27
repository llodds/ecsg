import os

RSFROOT=os.getenv('RSFROOT')
RSFROOT_LIB=os.path.join(RSFROOT, 'lib')
RSFROOT_INC=os.path.join(RSFROOT, 'include')


env = Environment(CXX='icpc',
                  CXXFLAGS=['-O3', '-qopenmp'],
                  CPPPATH=[RSFROOT_INC],
                  LIBS=['rsf++','rsf', 'm'],
                  LIBPATH=[RSFROOT_LIB],
                  LINKFLAGS=['-qopenmp'])
#   omp no longer works on OSX El Caption and later versions when using scons
#env = Environment(CXX='icpc',
#                  CXXFLAGS=['-O3'],
#                  CPPPATH=[RSFROOT_INC],
#                  LIBS=['rsf++','rsf', 'm'],
#                  LIBPATH=[RSFROOT_LIB],
#                  LINKFLAGS=[])


BuiltTarget = env.Program(target = 'main',
                          source = Glob('*.cpp'))
#   modify @rpath
env.AddPostAction(BuiltTarget,'install_name_tool -change @rpath/libiomp5.dylib /opt/intel/compilers_and_libraries_2017.0.102/mac/compiler/lib/libiomp5.dylib main')
