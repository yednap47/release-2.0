#!/bin/bash

rm */*.out
rm */*.tec
rm */*.dat
rm microbe_clogging/microbe_clogging.h5
rm -r results_benchmark

# Chrotran v1.0 benchmarks
cd abiotic
../../src/pflotran/chrotran -pflotranin abiotic.in >/dev/null 2>/dev/null
python abiotic.py

cd ../abiotic_mimt
../../src/pflotran/chrotran -pflotranin abiotic_mimt.in >/dev/null 2>/dev/null
python abiotic_mimt.py

cd ../microbe_biocide
../../src/pflotran/chrotran -pflotranin microbe_biocide.in >/dev/null 2>/dev/null
python microbe_biocide.py

cd ../microbe_enzymatic
../../src/pflotran/chrotran -pflotranin microbe_enzymatic.in >/dev/null 2>/dev/null
python microbe_enzymatic.py

cd ../microbe_growth_decay
../../src/pflotran/chrotran -pflotranin microbe_growth_decay.in >/dev/null 2>/dev/null
python microbe_growth_decay.py

cd ../microbe_growth_inhibitor
../../src/pflotran/chrotran -pflotranin microbe_growth_noinhibitor.in >/dev/null 2>/dev/null
../../src/pflotran/chrotran -pflotranin microbe_growth_inhibitor.in >/dev/null 2>/dev/null
python microbe_growth_inhibitor.py

cd ../microbe_respiration
../../src/pflotran/chrotran -pflotranin microbe_respiration.in >/dev/null 2>/dev/null
python microbe_respiration.py

cd ../microbe_clogging
../../src/pflotran/chrotran -pflotranin microbe_clogging.in >/dev/null 2>/dev/null
python microbe_clogging.py

# Chrotran v2.0 benchmarks
cd ../dithionite_s2o4_disp
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
mv ./*/*.png results_benchmark
echo results are in results_benchmark directory
