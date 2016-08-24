## **bamR**: Utilities for work wih **bam** files in **R**


**bamR** is an R package for generating two-dimensional heat maps useful for different genome-wide analyses.

## Download
To download this package you can use the Github user interface and just click the **Download** button. If you want to download the package from the terminal, then you need to install first a `git` client of your choice, using the package manager that is available for your system. For example, in Ubuntu or other Debian-based distributions of Linux, you can use `apt-get`:
```
$ sudo apt-get install git
```
while in Fedora, CentOS, or Red Hat Linux, you can use `yum`:
```
$ sudo yum install git
```
Installers for OSX (Mac) and Windows (PC) are available for download at the Git websites: `http://git-scm.com/download/mac` and `http://git-scm.com/download/win`.
After `git` has been installed, run the following command from the folder where you want to download the **plot2DO** package:
```
$ git clone https://github.com/rchereji/plot2DO.git
```


## Dependencies
**bamR** uses the following R packages: `caTools, colorRamps, GenomicAlignments, GenomicRanges, optparse`. To install these packages, open R and execute the following commands:
```{r}
> # Install caTools, colorRamps and optparse packages from CRAN:
> install.packages(c("caTools", "colorRamps", "optparse"))
>                  
> # Get the latest version of Bioconductor:
> source("https://bioconductor.org/biocLite.R")
> biocLite()
> 
> # Install the remaining Bioconductor packages:
> biocLite(c("GenomicAlignments", "GenomicRanges"))
```

## Workflow
After the R packages have been installed, the utilities included with **bamR** can be executed from `bash`. Below is a workflow example.

First you need to add the bamR folder to PATH. For example, you can copy the bamR folder to `/data/USER/scripts/bamR` and run in the terminal
```
$ export PATH=$PATH:/data/USER/scripts/bamR
```

Go to the folder where the bam files are located.
```
$ cd /data/USER/test/bam_files
$ ls *.bam
AGH01_WTI24.bam  AGH0220_1.bam  AGH0220_5.bam  AGH25_3.bam       Rpb3_AGH03_4.bam
AGH01_WTI44.bam  AGH0220_2.bam  AGH0220_6.bam  Rpb3_AGH03_1.bam  Rpb3_AGH03_5.bam
AGH01_WTU14.bam  AGH0220_3.bam  AGH08_5.bam    Rpb3_AGH03_2.bam  Rpb3_AGH03_6.bam
AGH01_WTU25.bam  AGH0220_4.bam  AGH09_5.bam    Rpb3_AGH03_3.bam
```

Check the length histogram of the sequenced reads, using `compute_length_histogram.R`. To see all the possible options run
```
$ compute_length_histogram.R --help
Usage: compute_length_histogram.R [options]


Options:
	-f FILE, --file=FILE
		Dataset file name (BAM format)

	-o OUTPUTS, --outputs=OUTPUTS
		Types of outputs to be generated, separated by commas only [options: pdf, csv, RData; default = pdf,csv,RData]

	-h, --help
		Show this help message and exit
```

Examples:
```
$ compute_length_histogram.R --file=AGH01_WTI44.bam
$ compute_length_histogram.R --file=AGH01_WTI24.bam
$ compute_length_histogram.R --file=AGH01_WTU25.bam
$ compute_length_histogram.R --file=AGH01_WTU14.bam
```

The results are put in the folder `Length_histograms`. These histograms are useful for selecting a reasonable size range for the “good reads” that are kept for further analyses. For example, for MNase-seq experiments, a reasonable interval of mono-nucleosomal sizes is [120, 160] bp.

To check the distribution of reads (occupancy/coverage or dyads/centers) at the gene ends, use `align_and_average.R`. To see all the possible options run
```
$ align_and_average.R --help
Usage: align_and_average.R [options]


Options:
	-f FILES, --files=FILES
		Dataset file names; replicates separated by commas [e.g.: -f rep_1.bam,rep_2.bam,rep_3.bam]

	--sampleLabel=SAMPLELABEL
		Sample label [default = Sample?]

	-t TYPE, --type=TYPE
		Types of averages to be computed, separated by commas only [options: occ, dyads; default = occ]

	-r REFERENCE, --reference=REFERENCE
		Reference points to align [options: TSS (default), TTS, Plus1]

	-l MINLENGTH, --minLength=MINLENGTH
		The smallest DNA fragment to be considered [default = 50]

	-L MAXLENGTH, --maxLength=MAXLENGTH
		The largest DNA fragment to be considered [default = 300]

	-u UPSTREAMPLOTWINDOW, --upstreamPlotWindow=UPSTREAMPLOTWINDOW
		Length of the upstream region that will be plotted [default = 1000]

	-d DOWNSTREAMPLOTWINDOW, --downstreamPlotWindow=DOWNSTREAMPLOTWINDOW
		Length of the downstream region that will be plotted [default = 1000]

	-p PRESORTEDLIST, --presortedList=PRESORTEDLIST
		Presorted list of genes (csv file)

	-h, --help
		Show this help message and exit
```

Examples:
```
$ align_and_average.R --files=AGH01_WTU14.bam,AGH01_WTU25.bam --sampleLabel=WTU --type=dyads --reference=Plus1 --minLength=120 --maxLength=160
$ align_and_average.R --files=AGH01_WTU14.bam,AGH01_WTU25.bam --sampleLabel=WTU --type=occ --reference=Plus1 --minLength=120 --maxLength=160
$ align_and_average.R --files=AGH01_WTI24.bam,AGH01_WTI44.bam --sampleLabel=WTI --type=dyads --reference=Plus1 --minLength=120 --maxLength=160
$ align_and_average.R --files=AGH01_WTI24.bam,AGH01_WTI44.bam --sampleLabel=WTI --type=occ --reference=Plus1 --minLength=120 --maxLength=160
```
The results are put in the folder `Avg_Dyads` or `Avg_Occ`, depending on the selected type of the plot (options: `--type=dyads` or `--type=occ`).

To compare different alignments use `plot_multiple_alignments.R`. To see all the possible options run
```
$ plot_multiple_alignments.R --help
Usage: plot_multiple_alignments.R [options]


Options:
	-f FILES, --files=FILES
		Alignment file names, separated by commas only [e.g. -f Avg_Occ_1.RData,Avg_Occ_2.Rdata]

	-l LABELS, --labels=LABELS
		Labels for the plot, separated by commas only [e.g. -l Exp_1,Exp_2]

	-t TYPE, --type=TYPE
		Type of average plot to be generated [options: occ, dyads; default = occ]

	-r REFERENCE, --reference=REFERENCE
		Reference points to align [options: TSS (default), TTS, Plus1]

	-h, --help
		Show this help message and exit
```

Examples:
```
$ cd Avg_Occ
$ ls *.RData
Avg_Occ_Plus1.WTI.120_160.RData  Avg_Occ_Plus1.WTU.120_160.Rdata

$ plot_multiple_alignments.R --files=Avg_Occ_Plus1.WTI.120_160.RData,Avg_Occ_Plus1.WTU.120_160.RData --labels=WTI,WTU --type=occ –reference=Plus1
$ cd ..
```

Using the aligned dyads, one can estimate the average nucleosome spacing, using `regression_analysis_dyads.R`.
Examples:
```
$ cd Avg_Dyads
$ ls *.RData
Avg_Dyads_Plus1.WTI.120_160.RData  Avg_Dyads_Plus1.WTU.120_160.Rdata
$ regression_analysis_dyads.R --file=Avg_Dyads_Plus1.WTU.120_160.RData
$ cd ..
```

To compare 2 samples by generating heat maps which are sorted by the difference between the 2 samples use `compare_heat_maps_treatment_vs_control.R`. To see all the possible options run
```
$ compare_heat_maps_treatment_vs_control.R --help
Usage: compare_heat_maps_treatment_vs_control.R [options]


Options:
	-c CONTROLFILES, --controlFiles=CONTROLFILES
		Control file names, separated by comma [e.g.: -c control_1.bam,control_2.bam,control_3.bam]

	-t TREATMENTFILES, --treatmentFiles=TREATMENTFILES
		Treatment file names, separated by comma [e.g.: -c treatment_1.bam,treatment_2.bam,treatment_3.bam]

	--controlLabel=CONTROLLABEL
		Control label [default = Control]

	--treatmentLabel=TREATMENTLABEL
		Treatment label [default = Treatment]

	-r REFERENCE, --reference=REFERENCE
		Reference points to align [options: TSS (default), TTS, Plus1]

	-p PRESORTEDLIST, --presortedList=PRESORTEDLIST
		Presorted list of genes (csv file)

	-l MINLENGTH, --minLength=MINLENGTH
		The smallest DNA fragment to be considered [default = 50]

	-L MAXLENGTH, --maxLength=MAXLENGTH
		The largest DNA fragment to be considered [default = 300]

	-u UPSTREAMPLOTWINDOW, --upstreamPlotWindow=UPSTREAMPLOTWINDOW
		Length of the upstream region that will be plotted [default = 1000]

	-d DOWNSTREAMPLOTWINDOW, --downstreamPlotWindow=DOWNSTREAMPLOTWINDOW
		Length of the downstream region that will be plotted [default = 1000]

	-a UPSTREAMSORTWINDOW, --upstreamSortWindow=UPSTREAMSORTWINDOW
		Length of the upstream region that will be used for sorting the differences [default = 500]

	-z DOWNSTREAMSORTWINDOW, --downstreamSortWindow=DOWNSTREAMSORTWINDOW
		Length of the downstream region that will be used for sorting the differences [default = 100]

	--colorScale=COLORSCALE
		Maximum color scale in single condition heat maps [default = 2]

	--colorScaleDiff=COLORSCALEDIFF
		Maximum color scale in difference heat map [default = 0.5]

	-h, --help
		Show this help message and exit
```

Example:
```
$ compare_heat_maps_treatment_vs_control.R --controlFiles=AGH0220_1.bam,AGH0220_2.bam,AGH0220_3.bam --controlLabel=WTU --treatmentFiles=AGH0220_4.bam,AGH0220_5.bam,AGH0220_6.bam --treatmentLabel=WTI --reference=TSS --minLength=50 --maxLength=300 --upstreamSortWindow=500 --downstreamSortWindow=100 --colorScale=2 --colorScaleDiff=0.5
```

In general, the parameters that have the default values can be omitted (this is true for all other R functions). The same result is obtained by executing the following command
```
$ compare_heat_maps_treatment_vs_control.R --controlFiles=AGH0220_1.bam,AGH0220_2.bam,AGH0220_3.bam --controlLabel=WTU --treatmentFiles=AGH0220_4.bam,AGH0220_5.bam,AGH0220_6.bam --treatmentLabel=WTI
```

`WTI_vs_WTU.50_300.csv` contains the average occupancies in the regions defined by `reference`, `upstreamSortWindow`, and `downstreamSortWindow`.

To get the 1000 genes that have the biggest change in occupancy, run
```
$ cd Heatmap_Occ
$ tail -n +4 WTI_vs_WTU.50_300.csv | sort -t ',' -g -k 4,4 | cut -d ',' -f 1 | head -n 1000 > WTI_vs_WTU.50_300.top1000.csv
```

Copy the list of top 1000 genes back to the parent folder
```
$ cp WTI_vs_WTU.50_300.top1000.csv ../
$ cd ..
```
Compare the top 1000 genes, using the extra argument `–presortedList=WTI_vs_WTU.50_300.top1000.csv`.
Examples:
1. Comparison WTI vs. WTU:
```
$ compare_heat_maps_treatment_vs_control.R --controlFiles=AGH0220_1.bam,AGH0220_2.bam,AGH0220_3.bam --controlLabel=WTU --treatmentFiles=AGH0220_4.bam,AGH0220_5.bam,AGH0220_6.bam --treatmentLabel=WTI --presortedList=WTI_vs_WTU.50_300.top1000.csv
```
2. Comparison gcn5 vs. WT:
```
$ compare_heat_maps_treatment_vs_control.R --controlFiles=AGH0220_4.bam,AGH0220_5.bam,AGH0220_6.bam --controlLabel=WT --treatmentFiles=AGH08_5.bam,AGH09_5.bam,AGH25_3.bam --treatmentLabel=gcn5 –presortedList=WTI_vs_WTU.50_300.top1000.csv
```

The results are put in the folder `Heatmap_Occ`.

To compute the average Rpb3 over the gene bodies of the same 1000 genes, use `compute_average_occ_over_gene_bodies.R`. To see all the possible options run
```
$ compute_average_occ_over_gene_bodies.R --help
Usage: compute_average_occ_over_gene_bodies.R [options]


Options:
	-f FILES, --files=FILES
		Dataset file names; replicates separated by commas [e.g.: -f rep_1.bam,rep_2.bam,rep_3.bam]

	--sampleLabel=SAMPLELABEL
		Sample label [default = Sample?]

	-l MINLENGTH, --minLength=MINLENGTH
		The smallest DNA fragment to be considered [default = 50]

	-L MAXLENGTH, --maxLength=MAXLENGTH
		The largest DNA fragment to be considered [default = 300]

	-p PRESORTEDLIST, --presortedList=PRESORTEDLIST
		Presorted list of genes (csv file)

	--colorScale=COLORSCALE
		Maximum color scale in single condition heat maps [default = 5]

	-h, --help
		Show this help message and exit
```

Examples:
```
$ compute_average_occ_over_gene_bodies.R --files=Rpb3_AGH03_1.bam,Rpb3_AGH03_2.bam,Rpb3_AGH03_3.bam --sampleLabel=WT --presortedList=WTI_vs_WTU.50_300.top1000.csv
$ Rscript compute_average_occ_over_gene_bodies.R --files=Rpb3_AGH03_4.bam,Rpb3_AGH03_5.bam,Rpb3_AGH03_6.bam --sampleLabel=WTI --presortedList=WTI_vs_WTU.50_300.top1000.csv
```

The results are put in the folder `Avg_Occ_over_gene_bodies`.

To plot the heat map for Rpb3, use `plot_single_heatmap.R. To see all the possible options run
```
$ plot_single_heat_map.R --help
Usage: plot_single_heat_map.R [options]


Options:
	-f FILES, --files=FILES
		Data file names; replicates separated by comma [e.g.: -c data_1.bam,data_2.bam,data_3.bam]

	--sampleLabel=SAMPLELABEL
		Sample label [default = Control]

	-t TYPE, --type=TYPE
		Type of heat map; multiple types separated by commas only [options: occ, dyads; default = occ]

	-r REFERENCE, --reference=REFERENCE
		Reference points to align [options: TSS (default), TTS, Plus1]

	-p PRESORTEDLIST, --presortedList=PRESORTEDLIST
		Presorted list of genes (csv file)

	-l MINLENGTH, --minLength=MINLENGTH
		The smallest DNA fragment to be considered [default = 50]

	-L MAXLENGTH, --maxLength=MAXLENGTH
		The largest DNA fragment to be considered [default = 300]

	-u UPSTREAMPLOTWINDOW, --upstreamPlotWindow=UPSTREAMPLOTWINDOW
		Length of the upstream region that will be plotted [default = 1000]

	-d DOWNSTREAMPLOTWINDOW, --downstreamPlotWindow=DOWNSTREAMPLOTWINDOW
		Length of the downstream region that will be plotted [default = 1000]

	--colorScale=COLORSCALE
		Maximum color scale in single condition heat maps [default = 2]

	-h, --help
		Show this help message and exit
```

Examples:
```
$ plot_single_heat_map.R --files=Rpb3_AGH03_1.bam,Rpb3_AGH03_2.bam,Rpb3_AGH03_3.bam --sampleLabel=Rpb3_WTU --presortedList=WTI_vs_WTU.50_300.top1000.csv --colorScale=5 --type=occ
$ plot_single_heat_map.R --files=Rpb3_AGH03_4.bam,Rpb3_AGH03_5.bam,Rpb3_AGH03_6.bam --sampleLabel=Rpb3_WTI --presortedList=WTI_vs_WTU.50_300.top1000.csv --colorScale=5 --type=occ
```

The same function can be used to plot a heat map of dyad distributions. Examples:
```
$ plot_single_heat_map.R --files=AGH01_WTU14.bam,AGH01_WTU25.bam --reference=Plus1 --sampleLabel=Nuc_WTU --minLength=120 --maxLength=160 --colorScale=2 --type=dyads 
$ plot_single_heat_map.R --files=AGH01_WTI24.bam,AGH01_WTI44.bam --reference=Plus1 --sampleLabel=Nuc_WTI --minLength=120 --maxLength=160 --colorScale=2 --type=dyads 
```

## License
**bamR** is freely available under the MIT License.
