AC_PREREQ([2.69])
AC_INIT([Ball], [0.1.0])

# Ensure a representative C source exists in sibling src/
AC_CONFIG_SRCDIR([../src/bcor.c])

AC_PROG_CC

# Copy C sources and headers from top-level src into R-package/src
AC_CONFIG_COMMANDS([copy-sources],
  [echo "* Copying C sources and headers from top-level src into R-package/src/..."
   PKGDIR="."
   SRCDIR="../src"
   DESTDIR="src"
   # Remove old copies (except ball_init.c)
   find "$DESTDIR" -maxdepth 1 -type f \( -name "*.c" ! -name "ball_init.c" -o -name "*.h" -o -name "*.o" -o -name ".so" \) -exec rm -f {} +
   # Copy new sources and headers
   cp "$SRCDIR"/*.c "$SRCDIR"/*.h "$DESTDIR"
   echo "* Copy complete."
  ],
  []
)

# Generate Makevars after copy
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT