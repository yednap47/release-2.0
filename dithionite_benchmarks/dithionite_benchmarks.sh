#!/bin/bash

rm */*.out
rm */*.tec
rm */*.dat

cd dithionite_s2o4_disp
../../src/pflotran/chrotran -pflotranin s2o4_disp.in >/dev/null 2>/dev/null
python s2o4_disp.py

cd ../dithionite_s2o4_o2
../../src/pflotran/chrotran -pflotranin s2o4_o2.in >/dev/null 2>/dev/null
python s2o4_o2.py

cd ../dithionite_s2o4_fe3
../../src/pflotran/chrotran -pflotranin s2o4_fe3.in >/dev/null 2>/dev/null
python s2o4_fe3.py

cd ../dithionite_fe2_o2
../../src/pflotran/chrotran -pflotranin fe2_o2.in >/dev/null 2>/dev/null
python fe2_o2.py

cd ../dithionite_fe2_cr6
../../src/pflotran/chrotran -pflotranin fe2_cr6.in >/dev/null 2>/dev/null
python fe2_cr6.py

cd ../dithionite_1d
mpirun -np 8 ../../src/pflotran/chrotran -pflotranin dithionite_1d.in >/dev/null 2>/dev/null
python dithionite_1d.py

cd ..

if [ -d "results_benchmark" ]; then
  rm -r results_benchmark
fi

mkdir results_benchmark
mv dithionite_s2o4_disp/s2o4_disp.png results_benchmark
mv dithionite_s2o4_o2/s2o4_o2.png results_benchmark
mv dithionite_s2o4_fe3/s2o4_fe3.png results_benchmark
mv dithionite_fe2_o2/fe2_o2.png results_benchmark
mv dithionite_fe2_cr6/fe2_cr6.png results_benchmark
mv dithionite_1d/dithionite_1d.png results_benchmark


echo results are in results_benchmark directory
