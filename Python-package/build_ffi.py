from cffi import FFI
import os
import shutil

# Initialize CFFI builder
ffibuilder = FFI()

# Determine paths
here = os.path.dirname(__file__)
# Ensure a local copy of src/ inside Python-package
dest_src = os.path.join(here, 'src')
root_src = os.path.abspath(os.path.join(here, '..', 'src'))
# Copy top-level src to python-package/src, overwriting if exists
if os.path.exists(dest_src):
    shutil.rmtree(dest_src)
shutil.copytree(root_src, dest_src)

# Use the copied src for building
src_root = dest_src

# 1) Load C function declarations from header
hdr_path = os.path.join(src_root, 'ball.h')
with open(hdr_path, 'r') as hdr_file:
    ffibuilder.cdef(hdr_file.read())

# 2) Automatically include all .c source files in src_root
source_files = []
for fname in os.listdir(src_root):
    if fname.endswith('.c'):
        source_files.append(os.path.join(src_root, fname))

ffibuilder.set_source(
    'ball._cball',            # Python module name
    '#include "ball.h"',     # Header inclusion
    sources=source_files,
    include_dirs=[src_root],
    extra_compile_args=['-O3'],
)

if __name__ == '__main__':
    ffibuilder.compile(verbose=True)