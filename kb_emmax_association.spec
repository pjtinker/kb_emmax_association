/*
A KBase module: kb_emmax_association
*/

module kb_emmax_association {

    /* An X/Y/Z style reference
    */
    typedef string obj_ref;

    /*
        required params:
            variation_obj_ref: reference to Variation Object
            output_file_prefix: prefix for EMMAX output files
    */

    typedef structure {
        string workspace_name;
        obj_ref variation_obj_ref;
        string output_file_prefix;
        int snp_return_count;
    } emmax_association_params;

    typedef structure {
        string report_name;
        string report_ref;
    } emmax_association_results;

    funcdef emmax_association(emmax_association_params)
        returns (emmax_association_results) authentication required;
};
