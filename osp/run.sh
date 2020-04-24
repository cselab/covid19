cd ../epidemics/cantons/py
python3 run_osp_cases.py 30 100 --TMCMC 0
cd ../../../osp
python3 main.py --path '../epidemics/cantons/py' --nMeasure 1 --nSensors 1
#cp result.npy ~/covid19_project/covid19/epidemics/cantons/py
