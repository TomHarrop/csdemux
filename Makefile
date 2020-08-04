graph: graph.svg

graph.svg: csdemux/Snakefile
	snakemake \
	-n \
	-s csdemux/Snakefile \
	--cores 8 \
	--dag \
	--forceall \
	--config \
	outdir=out \
	threads=8 \
	mem_gb=8 \
	r1_file=data/test_r1.fastq.gz \
	r2_file=data/test_r2.fastq.gz \
	samples_csv=data/example.samples.csv \
	| grep -v "^[[:space:]+]0" | grep -v "\->[[:space:]]0" \
	| dot -Tsvg \
	> graph.svg

readme: README.rst

README.rst: README.md
	pandoc -f markdown -t rst README.md > README.rst
