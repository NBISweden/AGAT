// Nextflow pipeline: chunk GFF by sizes and benchmark AGAT 
// Usage example:
//   nextflow run t/performance.nf --gff debug_cases/ensembl/Homo_sapiens.GRCh38.114.chr.1000000.gff3 --sizes 1000,5000,10000 --cpus 1,2,4

nextflow.enable.dsl=2

// Parameters
params.gff   = params.gff   ?: 'Homo_sapiens.GRCh38.114.chr.4171206.gff3'
params.sizes = params.sizes ?: '100000,500000,1000000,2000000,4171206'
params.cpus  = params.cpus  ?: '0,1,2,4,8'
params.help  = false

// Help message
def helpMessage() {
    log.info """
    ========================================
    AGAT Performance Benchmark Pipeline
    ========================================
    
    This pipeline benchmarks agat_convert_sp_gxf2gxf.pl performance by:
      1. Creating GFF subsets of different sizes
      2. Running conversions with different CPU counts
      3. Capturing timing metrics (user, sys, real, memory)
      4. Generating a summary table
    
    Usage:
      nextflow run t/performance.nf [options]
    
    Required parameters:
      --gff <file>         Input GFF/GFF3 file to benchmark
                           Default: ${params.gff}
    
    Optional parameters:
      --sizes <list>       Comma-separated list of line counts for subsets
                           Default: ${params.sizes}
                           Example: --sizes 500,1000,5000
      
      --cpus <list>        Comma-separated list of CPU counts to test
                           Default: ${params.cpus}
                           Example: --cpus 1,2,4,8
      
      --help               Show this help message and exit
    
    Output:
      t/perf_results/performance_summary.tsv  - Summary table with all metrics
      t/perf_results/metrics_*_cpu*.tsv       - Individual run metrics
    
    Examples:
      # Basic run with 2 sizes and 1 CPU
      nextflow run t/performance.nf --gff myfile.gff3 --sizes 1000,5000 --cpus 1
      
      # Full benchmark with multiple CPUs
      nextflow run t/performance.nf --gff myfile.gff3 --sizes 1000,5000,10000 --cpus 1,2,4,8
      
      # Resume a previous run
      nextflow run t/performance.nf --gff myfile.gff3 --sizes 1000,5000 --cpus 1,2 -resume
    
    ========================================
    """.stripIndent()
}

// Show help message if requested
if (params.help) {
    helpMessage()
    exit 0
}

// Normalize params into lists
def sizes_list = (params.sizes instanceof List) ? params.sizes.collect{ it as Integer } : params.sizes.toString().split(/[ ,]+/).findAll{ it }.collect{ it as Integer }
def cpus_list  = (params.cpus  instanceof List) ? params.cpus.collect{ it as Integer }  : params.cpus.toString().split(/[ ,]+/).findAll{ it }.collect{ it as Integer }

// Channels will be created within the workflow block for DSL2 correctness

// Create subsets of the input GFF retaining headers and first N non-header lines
process makeSubset {
	tag { "size:${size}" }
	input:
		val size
		path gff
	output:
		tuple val(size), path("chunk_${size}.gff3")
		script:
			"""
			set -eo pipefail
			head -n ${size} ${gff} > chunk_${size}.gff3
			"""
}

// Cross product: each chunk with each CPU value
// Combos will be handled via process inputs in DSL2 take/emit style

// Run conversion and capture detailed timing + AGAT memory footprint
process convertRun {
	tag { "size:${size}-cpu:${cpu}" }
	publishDir "t/perf_results", mode: 'copy', overwrite: true
	input:
		tuple val(size), path(chunk)
		each cpu
	output:
		path "metrics_*.tsv"
	script:
		"""
		set -euo pipefail

		out_gff=out_${size}_cpu${cpu}.gff3

		{
		time -p \
			agat_convert_sp_gxf2gxf.pl \
			-g $chunk -o \${out_gff} --cpu ${cpu} \
			1> agat.out 2> agat.err
		} 2> time.txt

		file_size=\$(du -h $chunk | cut -f1)
		user=\$(grep '^ *user' time.txt | cut -d' ' -f2)
		sys=\$(grep '^ *sys' time.txt | cut -d' ' -f2)
		real=\$(grep '^ *real' time.txt | cut -d' ' -f2)
		agat_mem=\$(grep -E 'Total memory:' agat.out | tail -n1 | sed -E 's/.*Total memory:\s([^\s]+)\s.*/\\1/')
		seq_nb=\$(awk '!/^#/ {print \$1}' \${out_gff} | sort -u | wc -l)

		printf "%d\t%s\t%d\t%s\t%s\t%s\t%s\t%s\n" "${size}" "\${file_size}" "${cpu}" "\${user}" "\${sys}" "\${real}" "\${agat_mem}" "\${seq_nb}" > "metrics_${size}_cpu${cpu}.tsv"
		"""
}

// Aggregate all metrics into a single summary table
process aggregate {
	tag { "summary" }
	publishDir "t/perf_results", mode: 'copy', overwrite: true
	input:
		path(ms)
	output:
		path "performance_summary.tsv"
	script:
		"""
		set -euo pipefail
		echo -e "size\tfile_size_mo\tcpu\tuser_s\tsys_s\telapsed_s\ttotal_memory_mo\tseq_nb" > performance_summary.tsv
		cat ${ms} >> performance_summary.tsv
		"""
}

workflow {
	sizes_ch = Channel.fromList(sizes_list)
	cpus_ch  = Channel.fromList(cpus_list)
	gff_ch   = Channel.fromPath(params.gff)

	chunks_ch  = makeSubset(sizes_ch, gff_ch.collect())
	metrics_ch = convertRun(chunks_ch, cpus_ch)
	summary_ch = aggregate(metrics_ch.collect())
}

