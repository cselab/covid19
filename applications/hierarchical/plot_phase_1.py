import sys
sys.path.append('../../')
from epidemics.data.files.canton_population import CANTON_LIST
import os
from subprocess import call

def create_folder(name):
    if not os.path.exists(name):
        os.makedirs(name)

def plot(model):
    shell_filename = 'plot_phase_1.sh'

    if os.path.exists(shell_filename):
        os.remove(shell_filename)

    with open(shell_filename, 'a') as file:
        file.write('#!/bin/sh \n')

    for canton in CANTON_LIST:
        samples_dir = './data/'+str(canton)+'/'+str(model)+'/_korali_samples/'
        output_dir  = './data/phase_1/'+str(model)+'/'
        output_file = output_dir +'posteriors/' +'/'+str(canton)+'.png'
        command = 'python -m korali.plotter --dir ' + str(samples_dir) + ' --output '+str(output_file)

        with open(shell_filename, 'a') as file:
            file.write(command+'\n')

        create_folder(output_dir)
        create_folder(output_dir+'posteriors/')
        create_folder(output_dir+'predictions/')

        call(['cp', './data/'+str(canton)+'/'+str(model)+'/figures/prediction.png',
                    output_dir+'/predictions/'+'/'+str(canton)+'.png'])

    call(['chmod', '744', shell_filename])
    rx = call('./'+shell_filename)
    os.remove(shell_filename)

if __name__ == "__main__":  

    if len(sys.argv) == 1:
        print('Default model: sir_altone_nbin')
        model = 'sir_altone_nbin'
    else:
        model = sys.argv[1] 

    plot(model)