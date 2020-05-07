import sys
sys.path.append('../../')
from epidemics.data.files.canton_population import CANTON_LIST, CANTON_LIST_SHORT
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
        
def plot_phase_1(model,regions):

    shell_filename = 'plot_phase_1.sh'
    create_shell_script(shell_filename)

    for canton in regions:
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

def plot_phase_3(model,regions):

    shell_filename = 'plot_phase_3.sh'
    create_shell_script(shell_filename)

    for canton in regions:
        phase_3_dir = './data/'+str(model)+'/phase_3_results/'
        samples_dir = phase_3_dir+str(canton)+'/'+str(model)+'/_korali_samples_phase_3/'
        output_file = phase_3_dir+str(canton)+'/'+str(model)+'/figures/posterior.png'
        create_folder(phase_3_dir+str(canton)+'/'+str(model)+'/figures')

        command = 'python -m korali.plotter --dir ' + str(samples_dir) + ' --output '+str(output_file)
        with open(shell_filename, 'a') as file:
            file.write(command+'\n')

    call(['chmod', '744', shell_filename])
    rx = call('./'+shell_filename)
    os.remove(shell_filename)

def plot_move_files(model,regions):

    for canton in regions:
        out_dir = './data/'+str(model)+'/figures/'+canton
        create_folder(out_dir)

        phase_1_dir = './data/'+str(model)+'/phase_1_results/'+canton+'/'+model+'/figures/'

        call(['cp', phase_1_dir+'prediction.png', out_dir+'/prediction_phase_1.png'])
        call(['cp', phase_1_dir+'posterior.png', out_dir+'/posterior_phase_1.png'])

        phase_3_dir = './data/'+str(model)+'/phase_3_results/'+canton+'/'+model

        call(['cp', phase_3_dir+'/figures/posterior.png', out_dir+'/posterior_phase_3.png'])
        call(['cp', phase_3_dir+'/'+canton+'/'+model+'/figures/prediction.png', out_dir+'/prediction_phase_3.png'])

if __name__ == "__main__":  

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--phases', '-ph', default='all', help='Which phases to plor.')
    parser.add_argument('--model', '-m', default='sir_altone_nbin', help='Model type')
    parser.add_argument('--regions', '-r', default='cantons_short', help='Model type')

    args = parser.parse_args()

    if args.regions == 'cantons':
        regions = CANTON_LIST
    elif args.regions == 'cantons_short':
        regions = CANTON_LIST_SHORT

    if args.phases == 'all':
        plot_phase_1(args.model,regions)
        plot_phase_2(args.model) 
        plot_phase_3(args.model,regions)  
    elif args.phases == '1':
        plot_phase_1(args.model,regions)
    elif args.phases == '2':
        plot_phase_2(args.model)
    elif args.phases == '3':
        plot_phase_3(args.model,regions)
    elif args.phases == 'move':
        plot_move_files(args.model,regions)




