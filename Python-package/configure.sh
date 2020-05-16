#!/bin/bash
cp distutils.cfg $1"/Lib/distutils"
cp $1"/vcruntime140.dll" $1"/libs"
cp cygwinccompiler.py $1"/Lib/distutils"
cd $1
gendef.exe "python"$2".dll"
dlltool --dllname "python"$2".dll" --def "python"$2".def" --output-lib "libpython"$2".a"
mv "libpython"$2".a" $1"/libs"
cp "python"$2".dll" $1"/libs"
mv "python"$2".def" $1"/libs"
gendef.exe vcruntime140.dll
dlltool -D vcruntime140.dll -d vcruntime140.def -l libvcruntime140.a
mv libvcruntime140.a $1"/libs"