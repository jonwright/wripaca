
cd ..
mingw32-make
cd python
python setup.py build install --home=test
cd test
set PYTHONPATH=lib/python
python basic_test.py

