# vOTU-related rules from "Snakefile" for acme-map4


rule votu_kma:
    # Map a given sample to the vOTU kma-index
    input:
        fq1 = get_fq1,
        fq2 = get_fq2
    output:
        votuKma = "20.votu/kma/kma.{sample}.tsv",
        votuReads = "20.votu/kma/votu.reads.{sample}.tsv.gz"
    log:
        summary = "logs/kma/{sample}.premapping",
        kma = "logs/kma/{sample}.log"
    threads: 20
    params:
        votudb = reference_paths["virusDB"],
        outfile = "20.votu/kma/kma.{sample}"
    resources:
        nodes = 1,
        ppn = 20,
        mem_gb = kma_mem_gb,
        walltime = 18000
    envmodules:
        "tools",
        "kma/1.4.7"
    shell:
        """
        echo $(date -Is) Starting KMA > {log.summary}
        kma -ipe "{input.fq1}"  "{input.fq2}" \
            -o {params.outfile} \
            -t_db {params.votudb} \
            -t {threads} \
            -tsv \
            -nc \
            -na \
            -mrs 0.01 \
            -apm p \
            &> {log.kma}

        # Save list of read-ids mapping to vOTU for downstream rarefying
        # NOTE hard-coded filter to remove HMR0 genes !!!
        # Also, the zgrep construct is to allow empty files to be processed without crashing
        echo $(date -Is) Creating reads file >> {log.summary}
        {{ zgrep -v "HMR0" 20.votu/kma/kma.{wildcards.sample}.frag.gz || test $? = 1; }} \
            | cut -f 7 | cut -d "/" -f 1 |  awk '!seen[$0]++' \
            | gzip -c > {output.votuReads}.tmp
        mv {output.votuReads}.tmp {output.votuReads}

        # Delete unnecessary KMA output files
        rm 20.votu/kma/kma.{wildcards.sample}.frag.gz
        rm 20.votu/kma/kma.{wildcards.sample}.res
        echo $(date -Is) Finished >> {log.summary}

        # kma arguments explanation 
        # -mrs: minimum read score normalized to alignment length 
        # -apm: p - reward for pairing mapped reads 
        # -nc : No consensus file
        # -na : No aln file
        """


rule votu_summary:
    input:
        *votu_summary_input
    output:
        '04.rds/%s.vOTU.counts.rds' % config["outputPrefix"],
        '04.rds/%s.vOTU.effectiveLength.rds' % config["outputPrefix"],
        '04.rds/%s.vOTU.comp.rds' % config["outputPrefix"],
        '04.rds/%s.vOTU.accounting.rds' % config["outputPrefix"],
        '20.votu/%s.vOTU.accounting.tsv' % config["outputPrefix"]
    log:
        "logs/votu_summary.log"
    threads: 1
    params:
        virusCatalogName = reference_paths["virusCatalogName"],
        virusObjectFile = reference_paths["virusObjectFile"],
        outputPrefix = config["outputPrefix"],
        virusAnnotation = reference_paths["virusAnnotation"],
        votuDirectory = "20.votu"
    resources:
        nodes = 1,
        ppn = 1,
        mem_gb = 4,
        walltime = 10800
    envmodules:
        "tools",
        "gcc/14.2.0",
        "intel/basekit/INITIALIZE/2023.0.0",
        "intel/basekit/mkl/2023.0.0",
        "R/4.5.0",
        "/home/projects/cu_10051/tools/Modules/modulefiles-shared/kmaprofiling/20230919"
    localrule: True
    shell:
        """
        mkdir -p 04.rds # If not already existing

        # Parse KMA files and write abundance matrices
        kmaprofiling --matrices \
            --votudir {params.votuDirectory} \
            --annodir {params.virusAnnotation} &> {log}

        # Convert files to RDS objects
        Rscript acme-map4/votu/votu_summary.R   \
            -v {params.virusObjectFile}     \
            -o {params.outputPrefix}             \
            -c {params.virusCatalogName}       \
            &>> {log}
        """

rule votu_kma_rarefy:
    # Map only rarefied vOTU reads to the vOTU kma-index
    input:
        fq1 = get_fq1,
        fq2 = get_fq2,
        votuReads = "20.votu/kma/votu.reads.{sample}.tsv.gz",
        accountingtable = '20.votu/%s.vOTU.accounting.tsv' % config["outputPrefix"]
    output:
        "21.votu.rarefied/kma/kma.{sample}.tsv"
    log:
        summary = "logs/kma/{sample}.rarefy.premapping",
        kma = "logs/kma/{sample}.rarefy.log"
    threads: 20
    params:
        votudb = reference_paths["virusDB"],
        outfile = "21.votu.rarefied/kma/kma.{sample}",
        downsizeAlgo = config["vOTU"]["downsizeAlgorithm"],
        downsizeTarget = config["vOTU"]["downsizeTarget"],
        readSuffixes = config["vOTU"]["readSuffixes"]
    resources:
        nodes = 1,
        ppn = 20,
        mem_gb = kma_mem_gb,
        walltime = 18000
    envmodules:
        "tools",
        "kma/1.4.7",
        "seqtk/1.3",
        "python/3.13.2"
    shell:
        """
        clean_up () {{
            tempSize=$( du -hs "$TMPDIR/{wildcards.sample}" | cut -f 1 )
            rm -Rf "$TMPDIR/{wildcards.sample}"
            echo Cleaned up temporary files ${{tempSize}} >> {log.summary}
        }}
        trap clean_up EXIT
        mkdir -p  "$TMPDIR/{wildcards.sample}"

        echo $(date -Is) Creating vOTU-specific fastq files > {log.summary}
        # Subset to vOTU reads from 1st round of kma-mapping
        gzip -dc {input.votuReads} > $TMPDIR/{wildcards.sample}/tmp.votu.reads

        if [ "{params.readSuffixes}" = "true" ]; then 
            awk 'NF{{print $0 "/1"}}' $TMPDIR/{wildcards.sample}/tmp.votu.reads \
                > $TMPDIR/{wildcards.sample}/votu.reads1
            awk 'NF{{print $0 "/2"}}' $TMPDIR/{wildcards.sample}/tmp.votu.reads \
                > $TMPDIR/{wildcards.sample}/votu.reads2
            seqtk subseq "{input.fq1}" \
                $TMPDIR/{wildcards.sample}/votu.reads1 > $TMPDIR/{wildcards.sample}/votu.reads.fq1 &
            seqtk subseq "{input.fq2}" \
                $TMPDIR/{wildcards.sample}/votu.reads2 > $TMPDIR/{wildcards.sample}/votu.reads.fq2 &
            wait
        else
            seqtk subseq "{input.fq1}" \
                $TMPDIR/{wildcards.sample}/tmp.votu.reads > $TMPDIR/{wildcards.sample}/votu.reads.fq1 &
            seqtk subseq "{input.fq2}" \
                $TMPDIR/{wildcards.sample}/tmp.votu.reads > $TMPDIR/{wildcards.sample}/votu.reads.fq2 &
            wait 
        fi

        # Rarefy - select the column with reads mapped to vOTU - sort and pick lowest number
        if [ "{params.downsizeTarget}" = "0" ]; then
            downsizeTarget=`python3 acme-map4/votu/votu_downsizeTarget.py \
                -a {params.downsizeAlgo} -i {input.accountingtable} | cut -f1`
        else
            downsizeTarget="{params.downsizeTarget}"
        fi

        echo $(date -Is) Sampling ${{downsizeTarget}} read pairs >> {log.summary}
        # Sample with a fixed seed ("-s100") to ensure matched pairs are retained
        seqtk sample -s100 $TMPDIR/{wildcards.sample}/votu.reads.fq1 ${{downsizeTarget}} \
            > $TMPDIR/{wildcards.sample}/votu.rarefied.reads.fq1 &
        seqtk sample -s100 $TMPDIR/{wildcards.sample}/votu.reads.fq2 ${{downsizeTarget}} \
            > $TMPDIR/{wildcards.sample}/votu.rarefied.reads.fq2 &
        wait

        if [ ! -s "$TMPDIR/{wildcards.sample}/votu.rarefied.reads.fq1" ] \
        || [ ! -s "$TMPDIR/{wildcards.sample}/votu.rarefied.reads.fq2" ]; then
            echo "No rarefied reads available; creating empty KMA output file." >> {log.summary}
            : > {output}
            echo $(date -Is) Finished >> {log.summary}
            exit 0
        fi

        echo $(date -Is) Starting KMA >> {log.summary}
        kma -ipe "$TMPDIR/{wildcards.sample}/votu.rarefied.reads.fq1" \
                 "$TMPDIR/{wildcards.sample}/votu.rarefied.reads.fq2" \
            -o {params.outfile} \
            -t_db {params.votudb} \
            -t {threads} \
            -tsv \
            -nc \
            -na \
            -nf \
            -mrs 0.01 \
            -apm p &> {log.kma}
        # -nf : No frag file

        # Delete unnecessary KMA output files
        rm -f {params.outfile}.res

        echo $(date -Is) Finished >> {log.summary}
        """


rule votu_summary_rarefy:
    input:
        expand("21.votu.rarefied/kma/kma.{sample}.tsv", sample=SAMPLES),
    output:
        '04.rds/%s.vOTU.dsCounts.rds' % config["outputPrefix"],
        '04.rds/%s.vOTU.dsComp.rds' % config["outputPrefix"],
        # ! '04.rds/%s.vOTU.accounting.rds' % config["outputPrefix"],
    log:
        "logs/votu_ds_summary.log"
    threads: 1
    params:
        downsizeTarget = config["vOTU"]["downsizeTarget"],
        downsizeAlgorithm = config["vOTU"]["downsizeAlgorithm"],
        virusCatalogName = reference_paths["virusCatalogName"],
        outputPrefix = config["outputPrefix"],
        virusObjectFile = reference_paths["virusObjectFile"],
        votuDirectory = "21.votu.rarefied",
        virusAnnotation = reference_paths["virusAnnotation"]
    resources:
        nodes = 1,
        ppn = 1,
        mem_gb = 4,
        walltime = 10800
    envmodules:
        "tools",
        "gcc/14.2.0",
        "intel/basekit/INITIALIZE/2023.0.0",
        "intel/basekit/mkl/2023.0.0",
        "R/4.5.0",
        "/home/projects/cu_10051/tools/Modules/modulefiles-shared/kmaprofiling/20230919"
    localrule: True
    shell:
        """
        # Parse KMA files and write abundance matrices
        kmaprofiling --matrices \
            --votudir {params.votuDirectory} \
            --annodir {params.virusAnnotation} &> {log}

        # Convert files to RDS objects
        Rscript acme-map4/votu/votu_summary_rarefy.R   \
            -t {params.downsizeTarget}                 \
            -a {params.downsizeAlgorithm}              \
            -o {params.outputPrefix}                   \
            -c {params.virusCatalogName}               \
            -v {params.virusObjectFile}                \
            &>> {log}
        """
