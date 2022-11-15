#!/usr/bin/env python
import os
import subprocess
import sys


'''
Helper functions to help run programs with the python subprocess module.
program - program filename. Expected to be in PATH.
arguments - a dictionary where the key value pairs define the arguments and options.

stdout and stderr are file handles
'''

def run(program, paired_arguments, single_arguments, stdout=None, stderr=None):
    ls = [program]
    if paired_arguments != None:
        for key, value in paired_arguments.items():
            ls.append(key)
            ls.append(value)
    
    if single_arguments != None:
        for argument in single_arguments:
            ls.append(argument)
    subprocess.run(ls, stdout=stdout, stderr=stderr)
'''
Template:
def run_xx():
    paired_arguments = {
        '-o' : output,
    }

    single_arguments = [
        in1,
        in2
    ]

    run(xx, paired_arguments, single_arguments)
'''

def run_minimap2_reads2ref(reads_path, ref_path, num_threads, output_path, minimap2_path='minimap2'):
    if os.path.exists(output_path):
        print(f'{output_path} already exists! Using existing.', file=sys.stderr)
        return
    with open(output_path, 'w') as sam_file:
            subprocess.run([minimap2_path, '-a', '--eqx', '-t', \
                str(num_threads), ref_path, reads_path], stdout=sam_file)

def run_hifiasm_hetero_reads_only(output_prefix, num_threads, list_of_reads_paths, reuse, bin_path='hifiasm'):
    hap1_fasta = output_prefix + '.hap1.fasta'
    hap2_fasta = output_prefix + '.hap2.fasta'
    if reuse and os.path.isfile(hap1_fasta) and os.path.isfile(hap2_fasta):
        return hap1_fasta, hap2_fasta
    named_arguments = {
        '-o' : output_prefix,
        '-t' : str(num_threads)
    }
    print('Assembling with hifiasm...', file=sys.stderr)
    run(bin_path, named_arguments, list_of_reads_paths)
    hap1_gfa = output_prefix + '.bp.hap1.p_ctg.gfa'
    hap2_gfa = output_prefix + '.bp.hap2.p_ctg.gfa'
    # awk '/^S/{print ">"$2;print $3}' test.p_ctg.gfa > test.p_ctg.fa
    with open(hap1_fasta, 'w') as hap1_fasta_handle, open(hap2_fasta, 'w') as hap2_fasta_handle:
        run('awk', None, ['/^S/{print ">"$2;print $3}', hap1_gfa], stdout=hap1_fasta_handle)
        run('awk', None, ['/^S/{print ">"$2;print $3}', hap2_gfa], stdout=hap2_fasta_handle)
    return hap1_fasta, hap2_fasta

def run_hifiasm_homozygous_reads_only(output_prefix, num_threads, paths_to_reads, reuse, bin_path='hifiasm'):
    
    contigs_fasta = output_prefix + '.contigs.fasta'
    if reuse and os.path.isfile(contigs_fasta):
        return contigs_fasta
    named_arguments = {
        '-o' : output_prefix,
        '-t' : str(num_threads)
    }

    print('Assembling with hifiasm...', file=sys.stderr)
    args = ['-l0']
    args += paths_to_reads
    #args.append(paths_to_reads)
    run(bin_path, named_arguments, args)
    contigs_gfa = output_prefix + '.bp.p_ctg.gfa'

    # awk '/^S/{print ">"$2;print $3}' test.p_ctg.gfa > test.p_ctg.fa
    with open(contigs_fasta, 'w') as contigs_fasta_handle:
        run('awk', None, ['/^S/{print ">"$2;print $3}', contigs_gfa], stdout=contigs_fasta_handle)
    return contigs_fasta


def run_pomoxis_assess_assm(assm_path, ref_path, num_threads, output_prefix='assm'):
    named_arguments = {
        '-i' : assm_path,
        '-r' : ref_path,
        '-d' : 'map-ont',
        '-c' : str(0),
        '-a' : '',
        '-t' : str(num_threads),
        '-p' : output_prefix
    }
    run('assess_assembly', named_arguments, None)

def run_quast(output_dir, assm_paths, ref_path, num_threads):
    named_arguments = {
        '-o' : output_dir,
        '-r' : ref_path,
        '-t' : str(num_threads)
    }

    run('quast.py', named_arguments, assm_paths)

