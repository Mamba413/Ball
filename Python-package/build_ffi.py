from cffi import FFI
import os

# Initialize CFFI builder
ffibuilder = FFI()

# Manually declare the C function prototypes to export
ffibuilder.cdef("""
    void bcor_test(double *, double *, double *, int *, int *, int *, int *, int *, int *, int *, int *, int *);
    void bcov_test(double *, double *, double *, double *, int *, int *, int *, int *);
    void kbcov_test(double *, double *, double *, int *, int *, int *, int *, int *);
    void bd_test(double *, double *, double *, int *, int *, int *, int *, int *, int *);
    void SRCT_new(double *, int *, int *, double *, int *, double *);
    void bd_gwas_screening(double *, double *, double *, int *, int *, double *, int *, int *, int *, int *, int *, int *, int *, int *);
    void bd_gwas_refining_single(double *, double *, double *, int *, int *, double *, int *, int *, int *, int *, int *, int *, int *, int *);
    void bdd_matrix_bias(double *, double *, int *, int *, int *);
    void bdd_matrix_bias_two_group(double *, double *, int *, int *, int *);
    void joint_kernel_matrix_bias(double *, double *, int *, int *, int *, int *);
    void cross_kernel_matrix_bias_crude(double *, double *, int *, int *, int *, int *);
""")

# Locate top-level src directory
here = os.path.dirname(os.path.abspath(__file__))
root_src = os.path.abspath(os.path.join(here, '..', 'src'))
if not os.path.isdir(root_src):
    raise FileNotFoundError(f"Cannot locate src directory at {root_src}")

# Collect all .c source files
sources = [os.path.join(root_src, f) for f in os.listdir(root_src) if f.endswith('.c')]

# Collect all .h header files for inclusion
headers = [f for f in os.listdir(root_src) if f.endswith('.h')]
include_block = '\n'.join(f'#include "{h}"' for h in headers)


# Configure extension module
ffibuilder.set_source(
    'Ball.cball',         # Python module name
    include_block,         # Generated include directives
    sources=sources,
    include_dirs=[root_src],
    extra_compile_args=['-O3'],
)

if __name__ == '__main__':
    ffibuilder.compile(verbose=True)