version 1.0

workflow NextstrainPhylogeny {
    input {
        # Sequence inputs - user provides ONE of these options
        File? multi_fasta_file
        Array[File]? individual_fasta_files
        
        # Required inputs
        File metadata_file
        File reference_file
        
        # Optional filtering parameters
        String? min_date
        String? max_date
        Int? sequences_per_group
        String? group_by
        File? exclude_strains_file
        
        # Runtime parameters
        String docker_image = "nextstrain/base:latest"
        Int cpu = 4
        String memory = "8 GB"
        Int disk_size = 50
    }
    
    # Validate inputs - exactly one sequence input must be provided
    if (!defined(multi_fasta_file) && !defined(individual_fasta_files)) {
        call ErrorTask as NoInputError { 
            input: 
                error_message = "ERROR: You must provide either 'multi_fasta_file' OR 'individual_fasta_files'"
        }
    }
    
    if (defined(multi_fasta_file) && defined(individual_fasta_files)) {
        call ErrorTask as BothInputsError { 
            input: 
                error_message = "ERROR: You cannot provide both 'multi_fasta_file' AND 'individual_fasta_files'. Choose one."
        }
    }
    
    if (defined(individual_fasta_files)) {
        call CombineFastas {
            input:
                fasta_files = select_first([individual_fasta_files]),
                docker_image = docker_image
        }
    }
    
    # Use either the provided multi-fasta or the newly created one
    File final_sequences = select_first([multi_fasta_file, CombineFastas.combined_fasta])
    
    call NextstrainPipeline {
        input:
            sequences_file = final_sequences,
            metadata_file = metadata_file,
            reference_file = reference_file,
            min_date = min_date,
            max_date = max_date,
            sequences_per_group = sequences_per_group,
            group_by = group_by,
            exclude_strains_file = exclude_strains_file,
            docker_image = docker_image,
            cpu = cpu,
            memory = memory,
            disk_size = disk_size
    }
    
    output {
        File nextstrain_aligned_fasta      = NextstrainPipeline.aligned_fasta
        File nextstrain_phylogenetic_tree  = NextstrainPipeline.phylogenetic_tree
        File nextstrain_ancestral_json     = NextstrainPipeline.ancestral_json
        File nextstrain_traits_json        = NextstrainPipeline.traits_json
        File nextstrain_branch_lengths     = NextstrainPipeline.branch_lengths
        File nextstrain_auspice_file       = NextstrainPipeline.auspice_json
    }
}

task ErrorTask {
    input {
        String error_message
    }
    
    command <<<
        echo "~{error_message}" >&2
        exit 1
    >>>
    
    runtime {
        docker: "ubuntu:20.04"
        memory: "1 GB"
        cpu: 1
    }
}

task CombineFastas {
    input {
        Array[File] fasta_files
        String docker_image
    }
    
    command <<<
        
        echo "Combining ~{length(fasta_files)} FASTA files into a single multi-FASTA"
        
        # Validate that we have at least one file
        if [ ~{length(fasta_files)} -eq 0 ]; then
            echo "ERROR: No FASTA files provided" >&2
            exit 1
        fi
        
        # Combine all FASTA files
        cat ~{sep=' ' fasta_files} > combined_sequences.fasta
        
        # Validate the combined file
        python3 -c "
        from Bio import SeqIO
        import sys
        
        try:
            sequences = list(SeqIO.parse('combined_sequences.fasta', 'fasta'))
            print(f'Successfully combined {len(sequences)} sequences')
            
            if len(sequences) == 0:
                print('ERROR: No sequences found in combined file', file=sys.stderr)
                sys.exit(1)
                
            # Check for duplicate IDs
            seq_ids = [seq.id for seq in sequences]
            duplicates = set([x for x in seq_ids if seq_ids.count(x) > 1])
            if duplicates:
                print(f'WARNING: Found duplicate sequence IDs: {duplicates}')
                
        except Exception as e:
            print(f'ERROR: Failed to parse combined FASTA: {e}', file=sys.stderr)
            sys.exit(1)
        "
        
        echo "FASTA combination completed successfully"
    >>>
    
    output {
        File combined_fasta = "combined_sequences.fasta"
    }
    
    runtime {
        docker: docker_image
        memory: "4 GB"
        cpu: 2
        disks: "local-disk 20 SSD"
    }
}

task NextstrainPipeline {
    input {
        File sequences_file
        File metadata_file
        File reference_file
        String? min_date
        String? max_date
        Int? sequences_per_group
        String? group_by
        File? exclude_strains_file
        String docker_image
        Int cpu
        String memory
        Int disk_size
    }
    
    command <<<   
        output_dir="results"     
        mkdir -p ${output_dir}
    
     
        # Copy exclude file if provided
        ~{if defined(exclude_strains_file) then "cp " + exclude_strains_file + " config/dropped_strains.txt" else ""}
        
        echo "Starting Nextstrain pipeline..."
        echo "Input validation and setup completed" > validation.log
        
        # Index sequences
        echo "Creating sequence index..."
        augur index \
            --sequences ~{sequences_file} \
            --output ${output_dir}/sequence_index.tsv
        
        # Build filter command dynamically
        filtered_data="${output_dir}/filtered.fasta"
        FILTER_CMD="augur filter --sequences ~{sequences_file} --metadata ~{metadata_file} --sequence-index ${output_dir}/sequence_index.tsv --output-sequences ${filtered_data} "
        
        ~{if defined(min_date) then 'FILTER_CMD="$FILTER_CMD --min-date ' + min_date + '"' else ""}
        ~{if defined(max_date) then 'FILTER_CMD="$FILTER_CMD --max-date ' + max_date + '"' else ""}
        ~{if defined(sequences_per_group) then 'FILTER_CMD="$FILTER_CMD --sequences-per-group ' + sequences_per_group + '"' else ""}
        ~{if defined(group_by) then 'FILTER_CMD="$FILTER_CMD --group-by ' + group_by + '"' else ""}
        ~{if defined(exclude_strains_file) then 'FILTER_CMD="$FILTER_CMD --exclude config/dropped_strains.txt"' else ""}
        
        echo "Filtering sequences with command: $FILTER_CMD"
        eval $FILTER_CMD
        
        
        # Align sequences
        echo "Aligning sequences to reference..."
        augur align \
            --sequences ${output_dir}/filtered.fasta \
            --reference-sequence ~{reference_file} \
            --output ${output_dir}/aligned.fasta \
            --fill-gaps
        
        # Build phylogeny
        echo "Building raw phylogenetic tree..."
        augur tree \
            --alignment ${output_dir}/aligned.fasta \
            --output ${output_dir}/tree_raw.nwk
        
        # Time-resolve tree
        echo "Refining tree with TreeTime..."
        augur refine \
            --tree ${output_dir}/tree_raw.nwk \
            --alignment ${output_dir}/aligned.fasta \
            --metadata ~{metadata_file} \
            --output-tree ${output_dir}/tree.nwk \
            --output-node-data ${output_dir}/branch_lengths.json \
            --timetree \
            --coalescent opt \
            --date-confidence \
            --clock-filter-iqd 4
        
        # Infer traits
        echo "Inferring ancestral traits..."
        augur traits \
            --tree ${output_dir}/tree.nwk \
            --metadata ~{metadata_file} \
            --columns country \
            --output-node-data ${output_dir}/traits.json
        
        # Reconstruct ancestral sequences
        echo "Reconstructing ancestral sequences..."
        augur ancestral \
            --alignment ${output_dir}/aligned.fasta \
            --tree ${output_dir}/tree.nwk \
            --output-node-data ${output_dir}/ancestral.json
        
        # Export for visualization
        echo "Exporting results for Auspice..."
        augur export v2 \
            --tree ${output_dir}/tree.nwk \
            --node-data ${output_dir}/branch_lengths.json \
            --node-data ${output_dir}/traits.json \
            --node-data ${output_dir}/ancestral.json \
            --metadata ~{metadata_file} \
            --output ${output_dir}/auspice.json
        
        echo "Pipeline completed successfully!" >> validation.log
        echo "🎉 Nextstrain pipeline completed successfully!"
    >>>
    
    output {
        File aligned_fasta = "results/aligned.fasta"
        File phylogenetic_tree = "results/tree.nwk"
        File traits_json = "results/traits.json"
        File ancestral_json = "results/ancestral.json"
        File auspice_json = "results/auspice.json"
        File branch_lengths = "results/branch_lengths.json"
    }
    
    runtime {
        docker: docker_image
        memory: memory
        cpu: cpu
        disks: "local-disk ${disk_size} SSD"
    }
}
