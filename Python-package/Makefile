cp README.md README.md
pandoc --from=markdown --to=rst --output=README.rst README.md

# renew cball.py when c code is updated
cd ./src
swig -python cball.i
cd ..
cp ./src/cball.py ./Ball/cball.py
python setup.py build

# generate .tar.gz package
python setup.py sdist
cp ./dist/Ball-0.1.0.tar.gz ./path_xxx