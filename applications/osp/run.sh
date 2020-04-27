cd ../../epidemics/cantons/py
python3 run_osp_cases.py --days 30 --samples 100 --TMCMC 0
cd ../../../applications/osp
python3 main.py --path '../../epidemics/cantons/py' --nMeasure 1 --nSensors 1 --Ntheta 100
#cp result.npy ~/covid19_project/covid19/epidemics/cantons/py
