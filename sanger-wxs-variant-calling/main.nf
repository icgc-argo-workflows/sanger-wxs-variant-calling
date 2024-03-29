#!/usr/bin/env nextflow
nextflow.enable.dsl = 2
name = 'sanger-wxs-variant-calling'
short_name = 'sanger-wxs'
version = '3.1.6-3.7.0'


/*
========================================================================================
                    ICGC ARGO Sanger WXS Variant Calling Workflow
========================================================================================
#### Homepage / Documentation
https://github.com/icgc-argo/sanger-wxs-variant-calling
#### Authors
Junjun Zhang @junjun-zhang <junjun.zhang@oicr.on.ca>
Linda Xiang @lindaxiang <linda.xiang@oicr.on.ca>
----------------------------------------------------------------------------------------

Required Parameters (no default):
--study_id                              SONG study ID
--tumour_aln_analysis_id                Tumour WXS sequencing_alignment SONG analysis ID
--normal_aln_analysis_id                Normal WXS sequencing_alignment SONG analysis ID
--ref_genome_fa                         Reference genome '.fa' file, secondary file ('.fa.fai') is expected to be under the same folder
--sanger_ref_genome_tar                 Tarball containing reference genome files from the same genome build
--sanger_vagrent_annot                  Tarball containing VAGrENT annotation reference
--sanger_ref_cnv_sv_tar                 Tarball containing CNV/SV reference
--sanger_ref_snv_indel_tar              Tarball containing SNV/Indel reference
--sanger_qcset_tar                      Tarball containing QC Genotype reference
--song_url                              SONG server URL
--score_url                             SCORE server URL
--api_token                             SONG/SCORE API Token

General Parameters (with defaults):
--cpus                                  cpus given to all process containers (default 1)
--mem                                   memory (GB) given to all process containers (default 1)

Download Parameters (object):
--download
{
    song_container_version              song docker container version, defaults set below
    score_container_version             score docker container version, defaults set below
    song_url                            song url for download process (defaults to main song_url param)
    score_url                           score url for download process (defaults to main score_url param)
    api_token                           song/score API token for download process (defaults to main api_token param)
    song_cpus
    song_mem
    score_cpus
    score_mem
    score_transport_mem
}

generateBas Parameters (object):
--generateBas
{
    container_version                   docker container version, defaults to unset
    cpus                                cpus for seqDataToLaneBam container, defaults to cpus parameter
    mem                                 memory (GB) for seqDataToLaneBam container, defaults to mem parameter
}

sangerWxsVariantCaller Parameters (object):
--sangerWxsVariantCaller
{
    container_version                   docker container version, defaults to unset
    cpus                                cpus for bwaMemAligner container, defaults to cpus parameter
    mem                                 memory (GB) for bwaMemAligner container, defaults to mem parameter
    exclude                             reference contigs to exclude, default: 'chrUn%,HLA%,%_alt,%_random,chrM,chrEBV'
    ploidy                              ploidy estimate of the genome, default: 2.0
    purity                              purity estimate of the genome, default: 1.0
    seq_format                          input aligned sequence format: default: cram
}

repackSangerResults Parameters (object):
--repackSangerResults
{
    container_version                   docker container version, defaults to unset
    cpus                                cpus for bamMergeSortMarkdup container, defaults to cpus parameter
    mem                                 memory (GB) for bamMergeSortMarkdup container, defaults to mem parameter
    library_strategy                    library strategy of input reads, default: WXS
}

payloadGenVariantCall (object):
--payloadGenVariantCall
{
    container_version                   docker container version, defaults to unset
    cpus                                cpus for align container, defaults to cpus parameter
    mem                                 memory (GB) for align container, defaults to mem parameter
}

Upload Parameters (object):
--upload
{
    song_container_version              song docker container version, defaults set below
    score_container_version             score docker container version, defaults set below
    song_url                            song url for upload process (defaults to main song_url param)
    score_url                           score url for upload process (defaults to main score_url param)
    api_token                           song/score API token for upload process (defaults to main api_token param)
    song_cpus                           cpus for song container, defaults to cpus parameter
    song_mem                            memory (GB) for song container, defaults to mem parameter
    score_cpus                          cpus for score container, defaults to cpus parameter
    score_mem                           memory (GB) for score container, defaults to mem parameter
    score_transport_mem                 memory (GB) for score_transport, defaults to mem parameter
    extract_cpus                        cpus for extract container, defaults to cpus parameter
    extract_mem                         memory (GB) extract score container, defaults to mem parameter
}

*/

params.study_id = ""
params.tumour_aln_analysis_id = ""
params.normal_aln_analysis_id = ""

// the following params if provided local files will be used
params.tumour_aln_metadata = "NO_FILE1"
params.tumour_aln_cram = "NO_FILE2"
params.tumour_extra_info = "NO_FILE3"
params.normal_aln_metadata = "NO_FILE4"
params.normal_aln_cram = "NO_FILE5"
params.normal_extra_info = "NO_FILE6"

// dir for outputs, must be set when running in local mode
params.publish_dir = ""
params.cleanup = true
params.cpus = 1
params.mem = 1
params.max_retries = 5  // set to 0 will disable retry
params.first_retry_wait_time = 1  // in seconds

// song/score setting
params.song_url = ""
params.song_container = "ghcr.io/overture-stack/song-client"
params.song_container_version = "5.0.2"
params.score_url = ""
params.score_container = "ghcr.io/overture-stack/score"
params.score_container_version = "5.9.0"
params.score_mem = 20
params.score_cpus = 8
params.score_force = false
params.api_token = ""

params.download = [:]
params.generateBas = [:]
params.sangerWxsVariantCaller = [:]
params.repackSangerResults = [:]
params.prepSangerSupplement = [:]
params.prepSangerQc = [:]
params.cavemanVcfFix = [:]
params.payloadGenVariantCall = [:]
params.uploadVariant = [:]
params.extractSangerCall = [:]

download_params = [
    'max_retries': params.max_retries,
    'first_retry_wait_time': params.first_retry_wait_time,
    'song_url': params.song_url,
    'song_container': params.song_container,
    'song_container_version': params.song_container_version,
    'song_cpus': params.cpus,
    'song_mem': params.mem,
    'score_url': params.score_url,
    'score_container': params.score_container,
    'score_container_version': params.score_container_version,
    'score_cpus' : params.score_cpus,
    'score_mem' : params.score_mem,
    'score_transport_mem' : params.score_mem, 
    'api_token': params.api_token,
    *:(params.download ?: [:])
]

generateBas_params = [
    'ref_genome_fa': '',
    *:(params.generateBas ?: [:])
]

sangerWxsVariantCaller_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'species': 'human',
    'assembly': 'GRCh38',
    'skipannot': true,
    'exclude': 'chrUn%,HLA%,%_alt,%_random,chrM,chrEBV',
    'ref_genome_tar': '',
    'vagrent_annot': '',
    'ref_snv_indel_tar': '',
    *:(params.sangerWxsVariantCaller ?: [:])
]

repackSangerResults_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'library_strategy': 'WXS',
    *:(params.repackSangerResults ?: [:])
]

prepSangerSupplement_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    *:(params.prepSangerSupplement ?: [:])
]
prepSangerQc_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    *:(params.prepSangerQc ?: [:])
]
cavemanVcfFix_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    *:(params.cavemanVcfFix ?: [:])
]

extractSangerCall_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    *:(params.extractSangerCall ?: [:])
]

payloadGenVariantCall_params = [
    'cpus': params.cpus,
    'mem': params.mem,
    'publish_dir': params.publish_dir,
    *:(params.payloadGenVariantCall ?: [:])
]

upload_params = [
    'max_retries': params.max_retries,
    'first_retry_wait_time': params.first_retry_wait_time,
    'song_url': params.song_url,
    'song_container': params.song_container,
    'song_container_version': params.song_container_version,
    'song_cpus': params.cpus,
    'song_mem': params.mem,
    'score_url': params.score_url,
    'score_container': params.score_container,
    'score_container_version': params.score_container_version,
    'score_force' : params.score_force,
    'score_cpus' : params.score_cpus,
    'score_mem' : params.score_mem,
    'score_transport_mem' : params.score_mem,
    'api_token': params.api_token,
    *:(params.upload ?: [:])
]


// Include all modules and pass params
include { SongScoreDownload as dnldT; SongScoreDownload as dnldN } from './wfpr_modules/github.com/icgc-argo-workflows/nextflow-data-processing-utility-tools/song-score-download@2.9.0/main.nf' params(download_params)
include { generateBas as basT; generateBas as basN; } from './modules/raw.githubusercontent.com/icgc-argo/variant-calling-tools/generate-bas.0.2.1.0/tools/generate-bas/generate-bas' params(generateBas_params)
include { sangerWxsVariantCaller as sangerWxs } from './modules/raw.githubusercontent.com/icgc-argo/variant-calling-tools/sanger-wxs-variant-caller.3.1.6-3/tools/sanger-wxs-variant-caller/sanger-wxs-variant-caller.nf' params(sangerWxsVariantCaller_params)
include { repackSangerResults as repack } from './modules/raw.githubusercontent.com/icgc-argo/variant-calling-tools/repack-sanger-results.0.2.0.0/tools/repack-sanger-results/repack-sanger-results' params(repackSangerResults_params)
include { cavemanVcfFix as cavemanFix } from './modules/raw.githubusercontent.com/icgc-argo/variant-calling-tools/caveman-vcf-fix.0.1.0.0/tools/caveman-vcf-fix/caveman-vcf-fix' params(cavemanVcfFix_params)
include { prepSangerSupplement as prepSupp } from './modules/raw.githubusercontent.com/icgc-argo/variant-calling-tools/prep-sanger-supplement.0.1.2.0/tools/prep-sanger-supplement/prep-sanger-supplement' params(prepSangerSupplement_params)
include { prepSangerQc as prepQc } from './modules/raw.githubusercontent.com/icgc-argo/variant-calling-tools/prep-sanger-qc.0.1.3.0/tools/prep-sanger-qc/prep-sanger-qc' params(prepSangerQc_params)
include { extractFilesFromTarball as extractVarSnv; extractFilesFromTarball as extractVarIndel; extractFilesFromTarball as extractQC } from './modules/raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/extract-files-from-tarball.0.2.0.0/tools/extract-files-from-tarball/extract-files-from-tarball' params(extractSangerCall_params)
include { payloadGenVariantCalling as pGenVarSnv; payloadGenVariantCalling as pGenVarIndel; payloadGenVariantCalling as pGenVarSupp; payloadGenVariantCalling as pGenQc } from "./wfpr_modules/github.com/icgc-argo-workflows/data-processing-utility-tools/payload-gen-variant-calling@0.7.0/main" params(payloadGenVariantCall_params)
include { SongScoreUpload as upSnv; SongScoreUpload as upIndel; SongScoreUpload as upQc; SongScoreUpload as upSupp} from './wfpr_modules/github.com/icgc-argo-workflows/nextflow-data-processing-utility-tools/song-score-upload@2.9.3/main.nf' params(upload_params)
include { cleanupWorkdir as cleanup } from './wfpr_modules/github.com/icgc-argo/data-processing-utility-tools/cleanup-workdir@1.0.0/main'
include { getSecondaryFiles } from './wfpr_modules/github.com/icgc-argo/data-processing-utility-tools/helper-functions@1.0.0/main'
include { payloadAddUniformIds as pAddIdT; payloadAddUniformIds as pAddIdN } from './wfpr_modules/github.com/icgc-argo/data-processing-utility-tools/payload-add-uniform-ids@0.1.1/main'


workflow SangerWxs {
    take:
        study_id
        tumour_aln_analysis_id
        normal_aln_analysis_id
        tumour_aln_metadata
        tumour_aln_cram
        tumour_extra_info
        normal_aln_metadata
        normal_aln_cram
        normal_extra_info

    main:
        local_mode = false

        tumour_aln_seq = Channel.from()
        tumour_aln_seq_idx = Channel.from()
        normal_aln_seq = Channel.from()
        normal_aln_seq_idx = Channel.from()

        if (tumour_aln_analysis_id && normal_aln_analysis_id) {
            log.info "Run the workflow using input metadata / sequencing data from SONG/SCORE servers"

            // download tumour aligned seq and metadata from song/score (analysis type: sequencing_alignment)
            dnldT(study_id, tumour_aln_analysis_id)
            tumour_aln_seq = dnldT.out.files.flatten().first()
            tumour_aln_seq_idx = dnldT.out.files.flatten().last()
            tumour_aln_meta = dnldT.out.analysis_json

            // download normal aligned seq and metadata from song/score (analysis type: sequencing_alignment)
            dnldN(study_id, normal_aln_analysis_id)
            normal_aln_seq = dnldN.out.files.flatten().first()
            normal_aln_seq_idx = dnldN.out.files.flatten().last()
            normal_aln_meta = dnldN.out.analysis_json
        } else if (
            !tumour_aln_metadata.startsWith('NO_FILE') && \
            !tumour_aln_cram.startsWith('NO_FILE') && \
            !tumour_extra_info.startsWith('NO_FILE') && \
            !normal_aln_metadata.startsWith('NO_FILE') && \
            !normal_aln_cram.startsWith('NO_FILE') && \
            !normal_extra_info.startsWith('NO_FILE')
        ) {
            if (!params.publish_dir) {
                exit 1, "When use local inputs, params.publish_dir must be specified."
            } else {
                log.info "Use local inputs, outputs will be in: ${params.publish_dir}"
            }

            local_mode = true

            tumour_aln_seq = file(tumour_aln_cram)
            tumour_aln_seq_idx = Channel.fromPath(getSecondaryFiles(tumour_aln_cram, ['crai', 'bai']))
            pAddIdT(file(tumour_aln_metadata), file(tumour_extra_info))
            tumour_aln_meta = pAddIdT.out.payload

            normal_aln_seq = file(normal_aln_cram)
            normal_aln_seq_idx = Channel.fromPath(getSecondaryFiles(normal_aln_cram, ['crai', 'bai']))
            pAddIdN(file(normal_aln_metadata), file(normal_extra_info))
            normal_aln_meta = pAddIdN.out.payload
        } else {
            exit 1, "To download input aligned seq files from SONG/SCORE, please provide `params.tumour_aln_analysis_id` and `params.normal_aln_analysis_id`.\n" +
                "Or please provide `params.tumour_aln_metadata`, `params.tumour_aln_cram`, `params.tumour_extra_info`, `params.normal_aln_metadata`, `params.normal_aln_cram` and `params.normal_extra_info` to use local files as input."
        }

        // generate Bas for tumour
        basT(
            'tumour', tumour_aln_seq, tumour_aln_seq_idx.collect(),
            file(params.generateBas.ref_genome_fa),
            Channel.fromPath(getSecondaryFiles(params.generateBas.ref_genome_fa, ['fai']), checkIfExists: true).collect())

        // generate Bas for normal
        basN(
            'normal', normal_aln_seq, normal_aln_seq_idx.collect(),
            file(params.generateBas.ref_genome_fa),
            Channel.fromPath(getSecondaryFiles(params.generateBas.ref_genome_fa, ['fai']), checkIfExists: true).collect())

        // run Sanger WXS
        sangerWxs(
            file(params.sangerWxsVariantCaller.ref_genome_tar),
            file(params.sangerWxsVariantCaller.vagrent_annot),
            file(params.sangerWxsVariantCaller.ref_snv_indel_tar),
            tumour_aln_seq,  // aln seq
            tumour_aln_seq_idx.collect(),  // idx
            basT.out.bas_file,  // bas
            normal_aln_seq,  // aln seq
            normal_aln_seq_idx.collect(),  // idx
            basN.out.bas_file  // bas
        )

        // repack results
        repack(sangerWxs.out.result_archive, 'WXS')

        // fix caveman VCF
        cavemanFix(repack.out.caveman)

        // extract variant calls
        extractVarSnv(cavemanFix.out.fixed_tar, 'flagged.muts')
        extractVarIndel(repack.out.pindel, 'flagged')

        // prepare variant call supplements
        prepSupp(cavemanFix.out.fixed_tar.concat(repack.out.pindel, sangerWxs.out.timings).collect())

        pGenVarSnv(
            normal_aln_meta, tumour_aln_meta,
            extractVarSnv.out.extracted_files,
            name, short_name, version
        )
        pGenVarIndel(
            normal_aln_meta, tumour_aln_meta,
            extractVarIndel.out.extracted_files,
            name, short_name, version
        )
        pGenVarSupp(
            normal_aln_meta, tumour_aln_meta,
            prepSupp.out.supplement_tar.collect(),
            name, short_name, version
        )

        prepQc(basN.out.bas_file_with_tn.concat(basT.out.bas_file_with_tn).collect())
        pGenQc(normal_aln_meta, tumour_aln_meta,
                 prepQc.out.qc_metrics_tar,
                 name, short_name, version)

        // upload variant results in paralllel
        if (!local_mode) {
            upSnv(study_id, pGenVarSnv.out.payload, pGenVarSnv.out.files_to_upload, '')
            upIndel(study_id, pGenVarIndel.out.payload, pGenVarIndel.out.files_to_upload, '')
            upSupp(study_id, pGenVarSupp.out.payload, pGenVarSupp.out.files_to_upload, '')
            upQc(study_id, pGenQc.out.payload, pGenQc.out.files_to_upload, '')
        }

        if (params.cleanup) {
            if (local_mode) {
                cleanup(
                    basT.out.bas_file.concat(
                        basN.out, sangerWxs.out, pGenVarSnv.out, pGenVarIndel.out,
                        pGenVarSupp.out, pGenQc.out, repack.out.pindel, cavemanFix.out.fixed_tar,
                        extractVarSnv.out.extracted_files, extractVarIndel.out.extracted_files,
                        prepSupp.out, prepQc.out).collect(),
                    true
                )
            } else {
                cleanup(
                    dnldT.out.files.concat(
                        dnldN.out, basT.out, basN.out, sangerWxs.out, pGenVarSnv.out, pGenVarIndel.out,
                        pGenVarSupp.out, pGenQc.out, repack.out.pindel, cavemanFix.out.fixed_tar,
                        extractVarSnv.out.extracted_files, extractVarIndel.out.extracted_files,
                        prepSupp.out, prepQc.out).collect(),
                    upSnv.out.analysis_id.concat(
                        upIndel.out.analysis_id, upSupp.out.analysis_id, upQc.out.analysis_id).collect()
                )
            }
        }

}


workflow {
    SangerWxs(
        params.study_id,
        params.tumour_aln_analysis_id,
        params.normal_aln_analysis_id,
        params.tumour_aln_metadata,
        params.tumour_aln_cram,
        params.tumour_extra_info,
        params.normal_aln_metadata,
        params.normal_aln_cram,
        params.normal_extra_info
    )
}
