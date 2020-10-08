FASTQ_DIR = fastq
CUT_DIR = cut
FASTQC_DIR = fastqc
REFERENCE_DIR = reference
REFERENCE_PREFIX = genome
BAM_DIR = bam
COUNT_DIR = count
GENES = genes/genes.gtf

ACCESSIONS = $(shell cat accessions.list)
FASTQ = $(addprefix $(FASTQ_DIR)/, $(addsuffix .fastq.gz, $(ACCESSIONS)))
CUT = $(addprefix $(CUT_DIR)/, $(addsuffix .fastq.gz, $(ACCESSIONS)))
FASTQC = $(addprefix $(FASTQC_DIR)/, $(addsuffix _fastqc.html, $(ACCESSIONS)))
BAM = $(addprefix $(BAM_DIR)/, $(addsuffix .bam, $(ACCESSIONS)))
INDEX = $(addprefix $(BAM_DIR)/, $(addsuffix .bam.bai, $(ACCESSIONS)))
COUNT = $(addprefix $(COUNT_DIR)/, $(addsuffix .tsv, $(ACCESSIONS)))

N_THREADS = 12

download: $(FASTQ)
fastqc: $(FASTQC)
align: $(BAM)
count: $(COUNT)

$(FASTQ): $(FASTQ_DIR)/%.fastq.gz:
	@mkdir -p $(FASTQ_DIR)
	@echo "Downloading fastq.gz file for $@."
	fastq-dump --gzip -O $(FASTQ_DIR) $*  

$(CUT): $(CUT_DIR)/%.fastq.gz: $(FASTQ_DIR)/%.fastq.gz
	@mkdir -p $(CUT_DIR)
	@echo "Cutting adapters for $@."
	cutadapt -j $(N_THREADS) -a "A{100}" -o $@ $<

$(FASTQC): $(FASTQC_DIR)/%_fastqc.html: $(CUT_DIR)/%.fastq.gz
	@mkdir -p $(FASTQC_DIR)
	fastqc --outdir $(FASTQC_DIR) -t $(N_THREADS) $<

%.fastq_view:
	gzip --stdout -d $*.fastq.gz

$(BAM): $(BAM_DIR)/%.bam: $(CUT_DIR)/%.fastq.gz
	@mkdir -p $(BAM_DIR)
	bowtie2 --local --very-sensitive-local -p $(N_THREADS) -x $(REFERENCE_DIR)/$(REFERENCE_PREFIX) $(word 2,$^) > $@

$(INDEX): $(BAM_DIR)/%.bam.bai: $(BAM_DIR)/%.bam
	samtools sort $< -o $<
	samtools index $<

$(COUNT): $(COUNT_DIR)/%.tsv: $(BAM_DIR)/%.bam.bai
	@mkdir -p $(COUNT_DIR)
	htseq-count -f bam -s no $(BAM_DIR)/$*.bam $(GENES) > $@