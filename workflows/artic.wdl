version 1.0

workflow run_artic {
	input {
		String accession
		Array[File] fastqs
		String primer_version
		Int min_length
		Int max_length
		Int min_reads

		String container_registry
	}

	call artic {
		input:
			accession = accession,
			fastqs = fastqs,
			primer_version = primer_version,
			min_length = min_length,
			max_length = max_length,
			min_reads = min_reads,
			container_registry = container_registry
	}

	output {
		File consensus_fa = artic.consensus_fa
		File vcf = artic.vcf
		File vcf_index = artic.vcf_index
		File bam = artic.bam
		File summary = artic.summary
	}

	meta {
		author: "Heather Ward"
		email: "heather@dnastack.com"
	}
}

task artic {
	input {
		String accession
		Array[File] fastqs
		String primer_version
		Int min_length
		Int max_length
		Int min_reads

		String container_registry
	}

	String primer_version_base = sub(primer_version, "ARTIC ", "")

	Int threads = 2

	command <<<
		set -e

		cp ~{sep=' ' fastqs} .

		# articGuppyPlex
		artic guppyplex \
			--min-length ~{min_length} \
			--max-length ~{max_length} \
			--output "~{accession}.guppyplex.fastq" \
			--directory .

		# articMinionMedaka
		num_guppyplex_reads=$(($(wc -l < "~{accession}.guppyplex.fastq") / 4))
		if [ "$num_guppyplex_reads" -lt ~{min_reads} ]; then
			echo "[ERROR] Num guppyplex reads [$num_guppyplex_reads] does not pass minimum read filter [~{min_reads}]"
			exit 1
		fi

		artic minion \
			--medaka \
			--minimap2 \
			--threads ~{threads} \
			--scheme-directory "$SCHEME_REPO/primer_schemes" \
			--read-file "~{accession}.guppyplex.fastq" \
			nCoV-2019/~{primer_version_base} \
			~{accession}

		# articRemoveUnmappedReads
		samtools view \
			-F4 \
			-o ~{accession}.mapped.sorted.bam \
			~{accession}.sorted.bam

		# makeQCCSV
		qc.py \
			--nanopore \
			--outfile ~{accession}.qc.tmp.csv \
			--sample ~{accession} \
			--ref "$SCHEME_REPO/primer_schemes/nCoV-2019/~{primer_version_base}/nCoV-2019.reference.fasta" \
			--bam ~{accession}.primertrimmed.rg.sorted.bam \
			--fasta ~{accession}.consensus.fasta

		# remove meta chars from the qc output
		tr -d $'\r' < ~{accession}.qc.tmp.csv > ~{accession}.qc.csv

		qc_pass=$(tail -1 "~{accession}.qc.csv" | awk -F ',' '{print $NF}')
		if [ ! "$qc_pass" = "TRUE" ]; then
			echo "[ERROR] Sample did not pass QC [$qc_pass]"
			cat ~{accession}.qc.csv
			exit 1
		fi

		# make sample summary zip, including some useful plots and metrics
		mkdir ~{accession}_summary
		mv ~{accession}-barplot.png ~{accession}-boxplot.png ~{accession}.alignreport.txt ~{accession}.depth.png ~{accession}.qc.csv ~{accession}_summary
		zip -r ~{accession}_summary.zip ~{accession}_summary

		# change VCF sample name from SAMPLE to accession
		# add processing pipeline, version to header
		chrom_line=$(bcftools view -h ~{accession}.pass.vcf.gz | tail -1)
		echo "SAMPLE ~{accession}" > new_samplename.txt
		bcftools view --no-version -h ~{accession}.pass.vcf.gz | sed '$ d' > header.txt
		echo -e "##processing_pipeline=https://github.com/connor-lab/ncov2019-artic-nf/tree/v$ARTIC_VERSION\n$chrom_line" >> header.txt
		bcftools reheader \
			-h header.txt \
			-s new_samplename.txt \
			~{accession}.pass.vcf.gz \
			-o ~{accession}.malformed.vcf.gz
		tabix ~{accession}.malformed.vcf.gz

		# Create a formatted VCF that can be imported into bigquery
		# The longranger output VCF has a few problems:
		#    - it includes a semicolon at the end of the info line
		#    - the PH field does not match its schema (sometimes outputs Floats but claims to be Integer); also has 4 (identical?) entries per allele
		# To import into bigquery in a way that can be flattened, we remove the semicolon from the info line; delete the PH field
		bcftools annotate --no-version \
			-x INFO/PH,FORMAT/UG,FORMAT/UQ \
			~{accession}.malformed.vcf.gz \
		> ~{accession}.vcf

		bgzip ~{accession}.vcf
		tabix ~{accession}.vcf.gz

		# validate VCF index
		count=0
		while [ "$count" -lt 3 ]; do
			validate_index.sh \
				-v "~{accession}.vcf.gz" \
				-i "~{accession}.vcf.gz.tbi"
			retc=$?
			if [ "$retc" -ne 0 ]; then
				tabix -f "~{accession}.vcf.gz"
			else
				break
			fi
			count=$((count + 1))
		done
		if [ "$retc" -ne 0 ]; then
			echo "Failed to validate index"
			echo "Index dump:"
			bgzip -d < "~{accession}.vcf.gz.tbi" | xxd
			echo; echo "VCF dump:"
			xxd "~{accession}.vcf.gz"
			exit 1
		fi
	>>>

	output {
		File consensus_fa = "~{accession}.consensus.fasta"
		File vcf = "~{accession}.vcf.gz"
		File vcf_index = "~{accession}.vcf.gz.tbi"
		File bam = "~{accession}.primertrimmed.rg.sorted.bam"
		File summary = "~{accession}_summary.zip"
	}

	runtime {
		docker: "~{container_registry}/artic:1.3.0"
		cpu: threads
		memory: "7.5 GB"
		disks: "local disk 200 HDD"
	}
}
