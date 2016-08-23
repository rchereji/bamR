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

## Usage
After the R packages have been installed, the utilities included with **bamR** can be executed from `bash`. Below is a workflow example.

## License
**bamR** is freely available under the MIT License.
