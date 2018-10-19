# -*- coding: utf-8 -*-
#BEGIN_HEADER
import os 
import subprocess
import uuid
import errno
import json

from kb_emmax_association.core import emmax_util
#END_HEADER


class kb_emmax_association:
    '''
    Module Name:
    kb_emmax_association

    Module Description:
    A KBase module: kb_emmax_association
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.scratch = config['scratch']
        self.callback_url = os.environ['SDK_CALLBACK_URL']
        #END_CONSTRUCTOR
        pass

    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
    def emmax_association(self, ctx, emmax_association_params):
        """
        emmax_association: run EMMAX association on SNP variation data
        :param emmax_association_params: instance of type "emmax_association_params" (required
            params: variation_obj_ref: KBase.Variations object reference
            output_file_prefix: User specified prefix for EMMAX output files
            snp_return_count: integer value that determines the top N SNPs found via EMMAX)
        """        
        #BEGIN emmax_association
        returnVal = {}

        print '--->\nRunning kb_emmax_association.emmax_association\nparams:'
        print json.dumps(emmax_association_params, indent=1)

        for key, value in emmax_association_params.iteritems():
            if isinstance(value, basestring):
                emmax_association_params[key] = value.strip()

        emmax_util_params = self.config
        emmax_util_params['token'] = ctx['token']
        emmax_util_params['SDK_CALLBACK_URL'] = self.callback_url

        emmax_runner = emmax_util.EmmaxUtil(emmax_util_params)

        try:
            returnVal = emmax_runner.run_emmax_association(emmax_association_params)
        except Exception as e:
            print("Error running EMMAX association")
            raise ValueError(e)

        #END emmax_association

        # At some point might do deeper type checking...
        if not isinstance(returnVal, dict):
            raise ValueError('Method import_variation return value ' +
                             'returnVal is not type dict as required.')
        # return the results
        return [returnVal]