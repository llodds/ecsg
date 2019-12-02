import os

RSFROOT=os.getenv('RSFROOT')
RSFROOT_LIB=os.path.join(RSFROOT, 'lib')
RSFROOT_INC=os.path.join(RSFROOT, 'include')


env = Environment(CXX='/usr/local/opt/llvm/bin/clang++',
                  CXXFLAGS=['-O3', '-fopenmp'],
                  CPPPATH=[RSFROOT_INC],
                  LIBS=['rsf++','rsf', 'm'],
                  LIBPATH=[RSFROOT_LIB, '/usr/local/opt/llvm/lib'],
                  LINKFLAGS=['-fopenmp'])


BuiltTarget = env.Program(target = 'main',
                          source = Glob('*.cpp'))
#   modify @rpath
env.AddPostAction(BuiltTarget,'install_name_tool -change @rpath/libiomp5.dylib /opt/intel/compilers_and_libraries_2017.0.102/mac/compiler/lib/libiomp5.dylib main')
