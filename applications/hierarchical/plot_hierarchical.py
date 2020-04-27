import sys
sys.path.append('../../')
from epidemics.data.files.canton_population import CANTON_LIST
import os
from subprocess import call
import argparse

def create_folder(name):
    if not os.path.exists(name):
        os.makedirs(name)

def create_shell_script(shell_filename):

    if os.path.exists(shell_filename):
        os.remove(shell_filename)

    with open(shell_filename, 'a') as file:
        file.write('#!/bin/sh \n')
        
def plot_phase_1(model):

    shell_filename = 'plot_phase_1.sh'
    create_shell_script(shell_filename)

    for canton in CANTON_LIST:
        phase_1_dir = './data/'+str(model)+'/phase_1_results/'
        samples_dir = phase_1_dir+str(canton)+'/'+str(model)+'/_korali_samples/'
        output_file = phase_1_dir+str(canton)+'/'+str(model)+'/figures/posterior.png'
        
        command = 'python -m korali.plotter --dir ' + str(samples_dir) + ' --output '+str(output_file)
        with open(shell_filename, 'a') as file:
            file.write(command+'\n')

    call(['chmod', '744', shell_filename])
    rx = call('./'+shell_filename)
    os.remove(shell_filename)

def plot_phase_2(model):

    shell_filename = 'plot_phase_2.sh'
    create_shell_script(shell_filename)

    phase_2_dir = './data/'+str(model)+'/phase_2_results/'
    
    samples_dir = phase_2_dir+'/_korali_samples/'
    output_file = phase_2_dir+'/figures/posterior.png'
    create_folder(phase_2_dir+'/figures')

    command = 'python -m korali.plotter --dir ' + str(samples_dir) + ' --output '+str(output_file)
    with open(shell_filename, 'a') as file:
        file.write(command+'\n')

    call(['chmod', '744', shell_filename])
    rx = call('./'+shell_filename)
    os.remove(shell_filename)

def plot_phase_3(model):

    shell_filename = 'plot_phase_3.sh'
    create_shell_script(shell_filename)

    for canton in CANTON_LIST:
        phase_3_dir = './data/'+str(model)+'/phase_3_results/'
        samples_dir = phase_3_dir+str(canton)+'/'+str(model)+'/_korali_samples/'
        output_file = phase_3_dir+str(canton)+'/'+str(model)+'/figures/posterior.png'
        create_folder(phase_3_dir+str(canton)+'/'+str(model)+'/figures')

        command = 'python -m korali.plotter --dir ' + str(samples_dir) + ' --output '+str(output_file)
        with open(shell_filename, 'a') as file:
            file.write(command+'\n')

    call(['chmod', '744', shell_filename])
    rx = call('./'+shell_filename)
    os.remove(shell_filename)

if __name__ == "__main__":  

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--phases', '-ph', default='all', help='Which phases to plor.')
    parser.add_argument('--model', '-m', default='sir_altone_nbin', help='Model type')
    args = parser.parse_args()

    if args.phases == 'all':
        plot_phase_1(args.model)
        plot_phase_2(args.model) 
        plot_phase_3(args.model)  
    elif args.phases == '1':
        plot_phase_1(args.model)
    elif args.phases == '2':
        plot_phase_2(args.model)
    elif args.phases == '3':
        plot_phase_3(args.model)





