#!/bin/bash nextflow

/*
 * Copyright (c) 2019-2020, Ontario Institute for Cancer Research (OICR).
 *                                                                                                               
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program. If not, see <https://www.gnu.org/licenses/>.
 */

/*
 * author Junjun Zhang <junjun.zhang@oicr.on.ca>
 */

nextflow.enable.dsl=2

params.tumour_aln_analysis_id = ""
params.normal_aln_analysis_id = ""

// the following four if provided local files will be used
params.tumour_aln_metadata = "NO_FILE1"
params.tumour_aln_cram = "NO_FILE2"
params.normal_aln_metadata = "NO_FILE3"
params.normal_aln_cram = "NO_FILE4"

params.publish_dir = ""  // dir for outputs, must be set when running in local mode


params.wf_name = ""
params.wf_short_name = ""
params.wf_version = ""
params.api_token = ""
params.song_url = ""
params.score_url = ""
params.cpus = 1
params.mem = 1  // GB
params.cleanup = true

include { SangerWxs } from "../main" params(params)


workflow {
  SangerWxs(
    params.study_id,
    params.tumour_aln_analysis_id,
    params.normal_aln_analysis_id,
    params.tumour_aln_metadata,
    params.tumour_aln_cram,
    params.normal_aln_metadata,
    params.normal_aln_cram,
  )
}
