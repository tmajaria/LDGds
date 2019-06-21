task calculate_LD {
	File script # /LDGds/calculate_LD.R
	File gds_file
	File? sample_ids_file
	String? ref_var
	String? interval
	Int? half_interval
	Int? min_mac
	Int? max_mac
	Float? min_maf
	Float? max_maf
	String? ld_method
	String out_pref
	Int? memory
	Int? disk

	Int default_disk = ceil(size(gds_file, "GB")) + 20

	command {
		R --vanilla --args ${gds_file} ${default="NA" sample_ids_file} ${default="NA" ref_var} ${default="NA" interval} ${default="50000" half_interval} ${default="0" min_mac} ${default="10000000" max_mac} ${default="0" min_maf} ${default="1" max_maf} ${default="r" ld_method} ${out_pref} < ${script}
	}

	runtime {
		docker: "analysiscommon/genesis_wdl:v0.1"
		disks: "local-disk " + select_first([disk,default_disk]) + " HDD"
		memory: select_first([memory,"5"]) + " GB"
	}

	output {
		File out_file = select_first(glob("${out_pref}*.csv"))
	}
}

workflow LD_wf {
	File this_gds_file
	File? this_sample_ids_file
	String? this_ref_var
	String? this_interval
	Int? this_half_interval
	Int? this_min_mac
	Int? this_max_mac
	Float? this_min_maf
	Float? this_max_maf
	String? this_ld_method
	String this_out_pref

	Int? this_memory
	Int? this_disk

	call calculate_LD {
		input: 
			gds_file = this_gds_file,
			sample_ids_file = this_sample_ids_file,
			ref_var = this_ref_var,
			interval = this_interval,
			half_interval = this_half_interval,
			min_mac = this_min_mac,
			max_mac = this_max_mac,
			min_maf = this_min_maf,
			max_maf = this_max_maf,
			ld_method = this_ld_method,
			out_pref = this_out_pref,
			memory = this_memory,
			disk = this_disk
	}
	

	output {
        File ld_file = calculate_LD.out_file
    }
}
