

```bash
mamba activate snake
cd /data/gpfs-1/users/martinl_c/work/hcpython/New
snakemake --latency-wait 300 --use-singularity --rerun-incomplete --singularity-args "--bind $HOME:$HOME --bind /data:/data --bind /fast:/fast"
```