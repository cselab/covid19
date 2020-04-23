cd ../epidemics/cantons/py
#python3 run_osp_cases.py 100 1000
cd ../../../osp
python3 main.py --path '../epidemics/cantons/py' --nMeasure 1 --nSensors 3
#cp result.npy ~/covid19_project/covid19/epidemics/cantons/py
