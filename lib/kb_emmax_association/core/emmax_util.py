import glob
import re
import time
import json
import os
import traceback
import uuid
import errno
import subprocess
import zipfile
import shutil
import csv

from pprint import pprint

import pandas as pd 
import numpy as np 

from DataFileUtil.DataFileUtilClient import DataFileUtil 
from KBaseReport.KBaseReportClient import KBaseReport 

PLINK_PRE = 'emmax_test_run'
TOP_SNP_FN = 'top_emmax_snps.ps'
TEMPLATE_DIRECTORY = '/kb/module/lib/kb_emmax_association/core/manhattan_plot/'

def log(message, prefix_newline=False):
    """Logging function, provides a hook to suppress or redirect log messages."""
    print(('\n' if prefix_newline else '') + '{0:.2f}'.format(time.time()) + ': ' + str(message))

class EmmaxUtil:

    def __init__(self, config):
        self.ws_url = config["workspace-url"]
        self.callback_url = config['SDK_CALLBACK_URL']
        self.token = config['token']
        self.shock_url = config['shock-url']
        self.scratch = os.path.join(config['scratch'], 'emmax_assoc_'+str(uuid.uuid4()))
        os.mkdir(self.scratch)
        self.config = config
        self.dfu = DataFileUtil(self.callback_url)
        self.kbr = KBaseReport(self.callback_url, token=self.config['token'])
    
    def _validate_phenotype_file(self, phenotype_filepath, tfam_filepath, case_control=False):
        #TODO: check if case/control.  If so, recode to 2/1 if necessary.  
        #TODO: verify the number of samples in pheno file matches tfam count
        #TODO: verify the order of the two files matches.  If not, reorder pheno
        #TODO: verify the FID and wiFID match the tfam
        pass


    def _validate_emmax_params(self, params):
        #TODO: All manner of param validation 
        pass

    def _create_tsv_file(self, top_snp_filepath, tsv_filename):
        log("Generating tsv file from {}".format(top_snp_filepath))
        cols = ['SNP', 'CHR', 'BP', 'P']
        snpData = pd.read_csv(top_snp_filepath, delimiter='\t', header=None, names=['ID', 'SE', 'p'])

        tsvData = pd.DataFrame(columns = cols)

        tsvData['SNP'] = snpData.ID 
        tsvData['CHR'] = snpData.ID.str[1:2]
        tsvData['P'] = snpData.p

        tsvData.sort_values(by='CHR', inplace=True)

        bps = []
        current_chr = tsvData.iloc[0]['CHR']
        count = 1
        for idx, row in tsvData.iterrows():
            if current_chr != row['CHR']:
                current_chr = row['CHR']
                count = 1
            bps.append(count)
            count = count + 1
        tsvData['BP'] = bps
        tsv_filepath = os.path.join(self.scratch, tsv_filename)
        tsvData.to_csv(tsv_filepath, sep='\t', index=False)        
        return tsv_filepath
        
    def _copyDirectory(self, src, dest):
            try:
                 shutil.copytree(src, dest)
                 # Directories are the same
            except shutil.Error as e:
                 print('Directory not copied. Error: %s' % e)
            # Any error saying that the directory doesn't exist
            except OSError as e:
                 print('Directory not copied. Error: %s' % e)

    def _run_subprocess(self, command, print_output=False, use_shell=False):
        log("Executing command:\n{}\n".format(command))
        p = subprocess.Popen(command,
                             cwd=self.scratch,
                             stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT,
                             shell=use_shell)
        output, errors = p.communicate()
        if print_output:
            print("Command output:\n{}".format(output))
        if errors:
            print("Command error output\n{}".format(errors))
        if p.returncode != 0:
            error_msg = "Error running command:\n{}\nReturn code: {}".format(command, p.returncode)
            raise ValueError(error_msg)

    def _download_variation_file(self, variation_obj_ref):
        log("Retrieving Variation Object {}...".format(variation_obj_ref))
        
        try:
            variation_shock_id = self.dfu.get_objects({
                'object_refs': [variation_obj_ref]
            })['data'][0]['data']['variation_file_reference']

            self.dfu.shock_to_file({
                'shock_id': variation_shock_id,
                'file_path': self.scratch,
                'unpack' : 'unpack'
            })
        except Exception as e:
            log("Error while retrieving Variation Object {}".format(variation_obj_ref))
            log(e)
            raise ValueError(e)

        variation_filename = [f for f in os.listdir(self.scratch) if f.endswith('.vcf')][0]
        variation_filepath = os.path.join(self.scratch, variation_filename)
        log("Variation file successfully downloaded to {}".format(variation_filepath))
        return variation_filepath 

    def _move_phenotype_data(self, pheno_filename):
        """ This is here until we get a Phenotype/Trait object working """
        pheno_filepath = os.path.join(self.scratch, pheno_filename)
        shutil.copy('/kb/module/data/' + pheno_filename, pheno_filepath)
        return pheno_filepath

    def _convert_vcf_to_plink(self, variation_filepath, fam_id = '--double-id', plink_file_prefix='plink_out'):
        #FIXME: This function should probably be in VCF utils
        log("Generating PLINK .tfam and .tped from VCF...")
        plink_cmd = ['plink']
        plink_cmd.append('--vcf')
        plink_cmd.append(variation_filepath)
        plink_cmd.append('--recode12')
        plink_cmd.append('transpose')
        plink_cmd.append('--output-missing-genotype')
        plink_cmd.append('0')
        plink_cmd.append(fam_id)
        plink_cmd.append('--out')
        plink_cmd.append(plink_file_prefix)

        self._run_subprocess(plink_cmd, print_output=True)

    def _generate_kinship_matrix(self, tped_prefix, matrix_type='BN'):
        #TODO: check for existence of files before attempting matrix generation

        log("Generating {} Kinship Matrix...".format(matrix_type))
        kinship_cmd = ['emmax-kin', '-v', '-d', '10']
        if (matrix_type == 'BN'):
            # kinship_cmd.append('-d')
            # kinship_cmd.append('10')
            pass
        elif (matrix_type == 'IBS'):
            kinship_cmd.append('-s')
            # kinship_cmd.append('-d')
            # kinship_cmd.append('10')
        else:
            log("Invalid matrix type specified.  Aborting")
            raise ValueError("Invalid matrix type specified")
        kinship_cmd.append(tped_prefix)
        self._run_subprocess(kinship_cmd, print_output=True)

        kinship_matrix_filename = [f for f in os.listdir(self.scratch) if f.endswith('.kinf')][0]
        kinship_matrix_filepath = os.path.join(self.scratch, kinship_matrix_filename)
        log("Variation file successfully downloaded to {}".format(kinship_matrix_filepath))
        return kinship_matrix_filepath

    def _emmax_association(self, plink_prefix, pheno_filepath, kinship_filepath, emmax_params):
        log("Running EMMAX association analysis...")
        emmax_cmd = ['emmax']
        for param in emmax_params:
            emmax_cmd.append(param)
        emmax_cmd.append('-t')
        emmax_cmd.append(plink_prefix)
        emmax_cmd.append('-p')
        emmax_cmd.append(pheno_filepath)
        emmax_cmd.append('-k')
        emmax_cmd.append(kinship_filepath)
        emmax_cmd.append('-o')
        emmax_cmd.append(plink_prefix)

        self._run_subprocess(emmax_cmd)
        emmax_filenames = [f for f in os.listdir(self.scratch) if f.endswith('.reml') or f.endswith('.ps')]
        return emmax_filenames
        
    def _select_top_snps(self, count, ps_filepath, output_filepath):
        #select_cmd = ['awk', '{print $NF, $0}', ps_filepath, "|", "sort", "-n", "|", "cut", "-f2-", "-d' '"]
        select_cmd = ["awk '{{print $NF,$0}}' {} | sort -n | cut -f2- -d' ' | sed -n -e '1,{}p' > {}".format(ps_filepath, str(count), output_filepath)]
        self._run_subprocess(select_cmd, print_output=False, use_shell=True)

    def _generate_output_files(self):
        log('Ziping EMMAX .reml and .ps files...')
        output_files = list()
        allowed_extensions = ['.ps', '.reml', '.tsv']
        result_file = os.path.join(self.scratch, 'emmax_results.zip')
        with zipfile.ZipFile(result_file, 'w',
                                zipfile.ZIP_DEFLATED,
                                allowZip64=True) as zip_file:
            for root, dirs, files in os.walk(self.scratch):
                for file in files:
                    if (file.endswith(tuple(allowed_extensions))):
                        if file in zip_file.namelist():
                            continue
                        zip_file.write(os.path.join(root, file), file)

        output_files.append({'path': result_file,
                             'name': os.path.basename(result_file),
                             'label': os.path.basename(result_file),
                             'description': 'File(s) generated by EMMAX'})
        log("Importer output generated: {}".format(output_files))
        return output_files

    def _generate_html_report(self, template_dir, tsv_filepath):
        log("Generating HTML report...")
        html_report = []

        output_dir = os.path.join(self.scratch, 'html')

        self._copyDirectory(template_dir, output_dir)

        result_file_path = os.path.join(output_dir, 'index.html')

        shutil.copyfile(tsv_filepath, os.path.join(output_dir, 'emmax_top.tsv'))
        report_shock_id = self.dfu.file_to_shock({'file_path': output_dir,
                'pack': 'zip'})['shock_id']

        html_report.append({'shock_id': report_shock_id,
                'name': os.path.basename(result_file_path),
                'label': os.path.basename(result_file_path),
                'description': 'Manhattan plot'})
        return html_report

    def run_emmax_association(self, params):
        variation_filepath = self._download_variation_file(params['variation_obj_ref'])
        self._convert_vcf_to_plink(variation_filepath, '--double-id', params['output_file_prefix'])
        
        pheno_filepath = self._move_phenotype_data('flcReordered.pheno')
        kinship_filepath = self._generate_kinship_matrix(params['output_file_prefix'])
        
        emmax_params = ['-v', '-d', '10']
        emmax_assoc_files = self._emmax_association(params['output_file_prefix'], pheno_filepath, kinship_filepath, emmax_params)
        
        full_result_filepath = os.path.join(self.scratch, params['output_file_prefix'] + '.ps')
        top_snp_filepath = os.path.join(self.scratch, TOP_SNP_FN)
        self._select_top_snps(params['snp_return_count'], full_result_filepath, top_snp_filepath)
        tsv_filepath = self._create_tsv_file(top_snp_filepath, 'emmax_top.tsv')
        output_html_files = self._generate_html_report(TEMPLATE_DIRECTORY, tsv_filepath)
        output_emmax_files = self._generate_output_files()
        report_params = {
            'message': '',
            'workspace_name' : params.get('workspace_name'),
            'file_links' : output_emmax_files,
            'html_links' : output_html_files,
            'direct_html_link_index' : 0,
            'html_window_height' : 333,
            'report_object_name' : 'emmax_assoc_html_report_' + str(uuid.uuid4())
        }
        kbr_output = self.kbr.create_extended_report(report_params)
        report_output = {
            'report_name': kbr_output['name'],
            'report_ref': kbr_output['ref'],
        }
        log("EMMAX report generated successfully!")
        return report_output