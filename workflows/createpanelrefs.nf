/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_createpanelrefs_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { GENS_PON                    } from '../subworkflows/local/gens_pon'
include { GERMLINECNVCALLER_COHORT    } from '../subworkflows/local/germlinecnvcaller_cohort'
include { BAM_CREATE_SOM_PON_GATK     } from '../subworkflows/nf-core/bam_create_som_pon_gatk'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { CNVKIT_BATCH                } from '../modules/nf-core/cnvkit/batch'
include { MULTIQC                     } from '../modules/nf-core/multiqc'

// Initialize file channels based on params, defined in the params.genomes[params.genome] scope
ch_dict                         = params.dict                        ? Channel.fromPath(params.dict).map { dict -> [[id:dict.baseName], dict]}.collect()
                                : Channel.empty()
ch_fai                          = params.fai                         ? Channel.fromPath(params.fai).map { fai -> [[id:fai.baseName], fai]}.collect()
                                : Channel.empty()
ch_fasta                        = params.fasta                       ? Channel.fromPath(params.fasta).map { fasta -> [[id:fasta.baseName], fasta]}.collect()
                                : Channel.empty()
// Initialize cnvkit specific parameters
ch_cnvkit_targets               = params.cnvkit_targets              ? Channel.fromPath(params.cnvkit_targets).map { targets -> [[id:targets.baseName], targets]}.collect()
                                : Channel.value([[id:'null'], []])
// Initialize germlinecnvcaller specific parameters
ch_gcnv_exclude_bed             = params.gcnv_exclude_bed            ? Channel.fromPath(params.gcnv_exclude_bed).map { exclude -> [[id:exclude.baseName], exclude]}.collect()
                                : Channel.value([[id:'null'], []])
ch_gcnv_exclude_interval_list   = params.gcnv_exclude_interval_list  ? Channel.fromPath(params.gcnv_exclude_interval_list).map { exclude -> [[id:exclude.baseName], exclude]}.collect()
                                : Channel.value([[id:'null'], []])
ch_gcnv_mappable_regions        = params.gcnv_mappable_regions       ? Channel.fromPath(params.gcnv_mappable_regions).collect()
                                : Channel.value([[id:'null'], []])
ch_gcnv_ploidy_priors           = params.gcnv_ploidy_priors          ? Channel.fromPath(params.gcnv_ploidy_priors).collect()
                                : Channel.empty()
ch_gcnv_target_bed              = params.gcnv_target_bed             ? Channel.fromPath(params.gcnv_target_bed).map { targets -> [[id:targets.baseName], targets]}.collect()
                                : Channel.value([[id:'null'], []])
ch_gcnv_target_interval_list    = params.gcnv_target_interval_list   ? Channel.fromPath(params.gcnv_target_interval_list).map { targets -> [[id:targets.baseName], targets]}.collect()
                                : Channel.value([[id:'null'], []])
ch_gcnv_segmental_duplications  = params.gcnv_segmental_duplications ? Channel.fromPath(params.gcnv_segmental_duplications).collect()
                                : Channel.value([[id:'null'], []])
// Initialize mutect2 specific parameters
ch_mutect2_target_bed           = params.mutect2_target_bed          ? Channel.fromPath(params.mutect2_target_bed).collect()
                                : Channel.value([])

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow CREATEPANELREFS {

    take:
    ch_samplesheet // channel: samplesheet read in from --input
    main:
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    if (params.tools && params.tools.split(',').contains('cnvkit')) {

        ch_samplesheet
            .map{ meta, bam, bai, cram, crai -> [meta + [id:'panel'], bam]}
            .groupTuple()
            .map {meta, bam -> [ meta, [], bam ]}
            .set { ch_cnvkit_input }

        CNVKIT_BATCH ( ch_cnvkit_input, ch_fasta, [[:],[]], ch_cnvkit_targets, [[:],[]], true )
        ch_versions = ch_versions.mix(CNVKIT_BATCH.out.versions)
    }

    if (params.tools && params.tools.split(',').contains('germlinecnvcaller')) {

        ch_samplesheet
            .map{meta, bam, bai, cram, crai ->
                if (bam)  return [ meta + [data_type:'bam'], bam, bai ]
                if (cram) return [ meta + [data_type:'cram'], cram, crai ]
            }
            .set { ch_germlinecnvcaller_input }

        GERMLINECNVCALLER_COHORT (  ch_dict,
                                    ch_fai,
                                    ch_fasta,
                                    ch_germlinecnvcaller_input,
                                    ch_gcnv_ploidy_priors,
                                    ch_gcnv_mappable_regions,
                                    ch_gcnv_segmental_duplications,
                                    ch_gcnv_target_bed,
                                    ch_gcnv_target_interval_list,
                                    ch_gcnv_exclude_bed,
                                    ch_gcnv_exclude_interval_list,
                                    params.gcnv_model_name )

        ch_versions = ch_versions.mix(GERMLINECNVCALLER_COHORT.out.versions)
    }

    if (params.tools && params.tools.split(',').contains('mutect2')) {

        ch_mutect2_input = ch_samplesheet.map{meta, bam, bai, cram, crai ->
            if (bam)    return [ meta + [data_type:'bam'], bam, bai, [] ]
            if (cram)   return [ meta + [data_type:'cram'], cram, crai, [] ]
        }

        BAM_CREATE_SOM_PON_GATK(ch_mutect2_input,
            ch_fasta,
            ch_fai,
            ch_dict,
            params.mutect2_pon_name,
            ch_mutect2_target_bed)

        ch_versions = ch_versions.mix(BAM_CREATE_SOM_PON_GATK.out.versions)

    }

    if (params.tools && params.tools.split(',').contains('gens')) {

        ch_samplesheet
            .map{meta, bam, bai, cram, crai ->
                if (bam)  return [ meta + [data_type:'bam'], bam, bai ]
                if (cram) return [ meta + [data_type:'cram'], cram, crai ]
            }
            .set { ch_gens_input }

        GENS_PON(ch_dict,
                ch_fai,
                ch_fasta,
                ch_gens_input,
                params.gens_pon_name )


        ch_versions = ch_versions.mix(GENS_PON.out.versions)
    }

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_ceatepanelrefs_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()


    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
