AC_PREREQ([2.61])
AC_INIT([Ball], [0.1.0])

AC_PROG_CC

# Copy C sources and headers from top-level src into R-package/src
cd src
echo '*****************************************************************************************'
if test -f Ball_omp.h; then
  echo "Complier code exists in src/ directory. Install R package according to the complier code."
  cd ..
else 
  echo "* Copying C sources and headers from top-level src into R-package/src/..."
  cd ..
  PKGDIR="."
  SRCDIR="../src"
  DESTDIR="src"
  # Remove old copies (except ball_init.c)
  find "$DESTDIR" -maxdepth 1 -type f \( -name "*.c" ! -name "ball_init.c" -o -name "*.h" -o -name "*.o" -o -name ".so" \) -exec rm -f {} +
  # Copy new sources and headers
  cp "$SRCDIR"/*.c "$SRCDIR"/*.h "$DESTDIR"
  echo "* Copy complete."
fi
echo '*****************************************************************************************'

# Generate Makevars after copy
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT