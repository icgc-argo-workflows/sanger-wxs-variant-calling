#!/usr/bin/env nextflow
nextflow.preview.dsl = 2
name = 'sanger-wxs-variant-calling'
short_name = 'sanger-wxs'
version = '0.1.0-dev'


/*
========================================================================================
                    ICGC-ARGO Sanger WGS Variant Calling Pipeline
========================================================================================
#### Homepage / Documentation
https://github.com/icgc-argo/sanger-wxs-variant-calling
#### Authors
Junjun Zhang @junjun-zhang <junjun.zhang@oicr.on.ca>
Linda Xiang @lindaxiang <linda.xiang@oicr.on.ca>
----------------------------------------------------------------------------------------

Required Parameters (no default):
--study_id                              SONG study ID
--tumour_aln_analysis_id                Tumour WGS sequencing_alignment SONG analysis ID
--normal_aln_analysis_id                Normal WGS sequencing_alignment SONG analysis ID
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
    library_strategy                    library strategy of input reads, default: WGS
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
params.api_token = ""
params.song_url = ""
params.score_url = ""
params.cleanup = true

params.cpus = 1
params.mem = 1

download_params = [
    'song_container_version': '4.0.0',
    'score_container_version': '3.0.1',
    'song_url': params.song_url,
    'score_url': params.score_url,
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
    'exclude': 'chrUn%,HLA%,%_alt,%_random,chrM,chrEBV',
    'ref_genome_tar': '',
    'vagrent_annot': '',
    'ref_snv_indel_tar': '',
    *:(params.sangerWxsVariantCaller ?: [:])
]

repackSangerResults_params = [
    'library_strategy': 'WXS',
    *:(params.repackSangerResults ?: [:])
]

extractSangerCall_params = [
    *:(params.extractSangerCall ?: [:])
]

payloadGenVariantCall_params = [
    *:(params.payloadGenVariantCall ?: [:])
]

payloadGenQcMetrics_params = [
    *:(params.payloadGenVariantCall ?: [:])
]

upload_params = [
    'song_container_version': '4.0.0',
    'score_container_version': '3.0.1',
    'song_url': params.song_url,
    'score_url': params.score_url,
    'api_token': params.api_token,
    *:(params.uploadVariant ?: [:])
]


// Include all modules and pass params
include { songScoreDownload as dnldT; songScoreDownload as dnldN } from './song-score-utils/song-score-download' params(download_params)
include { generateBas as basT; generateBas as basN; } from './modules/raw.githubusercontent.com/icgc-argo/variant-calling-tools/generate-bas.0.2.0.0/tools/generate-bas/generate-bas' params(generateBas_params)
include sangerWxsVariantCaller as sangerWxs from './modules/raw.githubusercontent.com/icgc-argo/variant-calling-tools/sanger-wxs-variant-caller.3.1.6-2/tools/sanger-wxs-variant-caller/sanger-wxs-variant-caller.nf' params(sangerWxsVariantCaller_params)
include repackSangerResults as repack from './modules/raw.githubusercontent.com/icgc-argo/variant-calling-tools/repack-sanger-results.0.2.0.0/tools/repack-sanger-results/repack-sanger-results' params(repackSangerResults_params)
include { extractFilesFromTarball as extractVarSnv; extractFilesFromTarball as extractVarIndel; extractFilesFromTarball as extractQC } from './modules/raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/extract-files-from-tarball.0.2.0.0/tools/extract-files-from-tarball/extract-files-from-tarball' params(extractSangerCall_params)
include { payloadGenVariantCalling as pGenVar } from "./modules/raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/payload-gen-variant-calling.0.1.0.0/tools/payload-gen-variant-calling/payload-gen-variant-calling" params(payloadGenVariantCall_params)
//include { payloadGenSangerQC as pGenQC } from "./modules/raw.githubusercontent.com/icgc-argo/data-processing-utility-tools/payload-gen-sanager-qc.0.1.0.0/tools/payload-gen-sanager-qc/payload-gen-sanager-qc" params(payloadGenSangerQc_params)
include { songScoreUpload as upVar; songScoreUpload as upQC} from './song-score-utils/song-score-upload' params(upload_params)
include cleanupWorkdir as cleanup from './modules/raw.githubusercontent.com/icgc-argo/nextflow-data-processing-utility-tools/b45093d3ecc3cb98407549158c5315991802526b/process/cleanup-workdir'


def getSecondaryFiles(main_file, exts){  //this is kind of like CWL's secondary files
  def all_files = []
  for (ext in exts) {
    all_files.add(main_file + ext)
  }
  return all_files
}


workflow SangerWxs {
    take:
        study_id
        tumour_aln_analysis_id
        normal_aln_analysis_id

    main:
        // download tumour aligned seq and metadata from song/score (analysis type: sequencing_alignment)
        dnldT(study_id, tumour_aln_analysis_id)

        // download normal aligned seq and metadata from song/score (analysis type: sequencing_alignment)
        dnldN(study_id, normal_aln_analysis_id)

        // generate Bas for tumour
        basT(
            dnldT.out.files.flatten().first(), dnldT.out.files.flatten().last(),
            file(params.generateBas.ref_genome_fa),
            Channel.fromPath(getSecondaryFiles(params.generateBas.ref_genome_fa, ['.fai']), checkIfExists: true).collect())

        // generate Bas for normal
        basN(
            dnldN.out.files.flatten().first(), dnldN.out.files.flatten().last(),
            file(params.generateBas.ref_genome_fa),
            Channel.fromPath(getSecondaryFiles(params.generateBas.ref_genome_fa, ['.fai']), checkIfExists: true).collect())

        // run Sanger WXS
        sangerWxs(
            file(params.sangerWxsVariantCaller.ref_genome_tar),
            file(params.sangerWxsVariantCaller.vagrent_annot),
            file(params.sangerWxsVariantCaller.ref_snv_indel_tar),
            dnldT.out.files.flatten().first(),  // aln seq
            dnldT.out.files.flatten().last(),  // idx
            basT.out.bas_file,  // bas
            dnldN.out.files.flatten().first(),  // aln seq
            dnldN.out.files.flatten().last(),  // idx
            basN.out.bas_file  // bas
        )

        // repack results
        repack(sangerWxs.out.result_archive, 'WXS')

        // extract variant calls
        extractVarSnv(repack.out.caveman, 'flagged.muts.vcf.gz')
        extractVarIndel(repack.out.pindel, 'flagged.vcf.gz')

        pGenVar(
            dnldT.out.song_analysis, dnldN.out.song_analysis,
            extractVarSnv.out[0].concat(extractVarIndel.out[0]).flatten(),
            name, short_name, version
        )

        // upload variant results in paralllel
        upVar(study_id, pGenVar.out[0], extractVarSnv.out[0])   // to be fixed here

        /*  // more to flesh out
        qc_result_patterns = Channel.from(
            '???', '???')
        extractQC(repack.out.collect(), qc_result_patterns.flatten())

        payloadGenSangerQC(dnldT.out.song_analysis, extractQC.out)

        upQC(study_id, payloadGenSangerQC.out.payload, extractQC.out)
        */


        if (params.cleanup) {
            cleanup(
                dnldT.out.files.concat(dnldN.out, basT.out, basN.out, sangerWxs.out,
                    repack.out, extractVarSnv.out, extractVarIndel.out).collect(),
                upVar.out.analysis_id.collect())
        }

}


workflow {
    SangerWxs(
        params.study_id,
        params.tumour_aln_analysis_id,
        params.normal_aln_analysis_id
    )
}