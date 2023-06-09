include { check_mandatory_parameter; check_optional_parameters; check_parameter_value } from './params_utilities.nf'

def default_params(){
    /***************** Setup inputs and channels ************************/
    def params = [:] as nextflow.script.ScriptBinding$ParamsMap
    // Defaults for configurable variables
    params.help = false
    params.version = false
    params.reads = false
    params.database = false
    params.taxon = false
    params.sequencing_date = false
    params.output_dir = false
    return params
}

def check_params(Map params) { 
    final_params = params
    
    // set up reads files
    final_params.reads = check_mandatory_parameter(params, 'reads')
     
    // set up output directory
    final_params.output_dir = check_mandatory_parameter(params, 'output_dir') - ~/\/$/
     
    // set up database
    final_params.database = check_mandatory_parameter(params, 'database') - ~/\/$/
    
    // set up sequencing date
    final_params.sequencing_date = check_mandatory_parameter(params, 'sequencing_date')
	
    // set up taxon
    final_params.taxon = check_mandatory_parameter(params, 'taxon') - ~/\/$/
        
    // check taxon is valid
    final_params.taxon = check_parameter_value('taxon', final_params.taxon, ['S', 'G', 'F', 'O', 'C', 'P', 'K', 'D', 'U'])
    
    return final_params
}

