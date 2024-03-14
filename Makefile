DAG:
	snakemake --rulegraph | dot -Tsvg > dag.svg
	convert dag.svg dag.pdf