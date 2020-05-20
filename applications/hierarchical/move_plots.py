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

def plot_move_files(data_path,model,regions):

    for canton in regions:
        out_dir = data_path+str(model)+'/figures/'+canton
        create_folder(out_dir)

        phase_1_dir = data_path+str(model)+'/phase_1_results/'+canton+'/'+model+'/figures/'

        call(['cp', phase_1_dir+'prediction.png', out_dir+'/prediction_phase_1.png'])
        call(['cp', phase_1_dir+'posterior.png', out_dir+'/posterior_phase_1.png'])

        phase_3_dir = data_path+str(model)+'/phase_3_results/'+canton+'/'+model

        call(['cp', phase_3_dir+'/figures/posterior.png', out_dir+'/posterior_phase_3.png'])
        call(['cp', phase_3_dir+'/figures/prediction.png', out_dir+'/prediction_phase_3.png'])

if __name__ == "__main__":  

    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--model', '-m', default='sir_altone_nbin', help='Model type')
    parser.add_argument('--regions', '-r', default='cantons_short', help='Model type')
    parser.add_argument('--dir', '-dir', default='./data/', help='Model type')

    args = parser.parse_args()

    if args.regions == 'cantons':
        regions = CANTON_LIST
    elif args.regions == 'cantons_short':
        regions = CANTON_LIST_SHORT

    if args.phases == 'all':
        plot_phase_1(args.dir,args.model,regions)
        plot_phase_2(args.dir,args.model) 
        plot_phase_3(args.dir,args.model,regions)  
    elif args.phases == '1':
        plot_phase_1(args.dir,args.model,regions)
    elif args.phases == '2':
        plot_phase_2(args.dir,args.model)
    elif args.phases == '3':
        plot_phase_3(args.dir,args.model,regions)
    elif args.phases == 'move':
        plot_move_files(args.dir,args.model,regions)




