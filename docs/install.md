# Installation

## Using Docker

      
First you must have [Docker](https://docs.docker.com/get-docker/) installed and running.  
Secondly have look at the availabe AGAT biocontainers at [quay.io](https://quay.io/repository/biocontainers/agat?tab=tags).

Then:
```
# get the chosen AGAT container version
docker pull quay.io/biocontainers/agat:0.8.0--pl5262hdfd78af_0
# use an AGAT's tool e.g. agat_convert_sp_gxf2gxf.pl
docker run quay.io/biocontainers/agat:0.8.0--pl5262hdfd78af_0 agat_convert_sp_gxf2gxf.pl --help
```

## Using Singularity
      
First you must have [Singularity](https://sylabs.io/guides/3.5/user-guide/quick_start.html) installed and running.  
Secondly have look at the availabe AGAT biocontainers at [quay.io](https://quay.io/repository/biocontainers/agat?tab=tags).

Then:
```
# get the chosen AGAT container version
singularity pull docker://quay.io/biocontainers/agat:1.0.0--pl5321hdfd78af_0
# run the container
singularity run agat_1.0.0--pl5321hdfd78af_0.sif
```

You are now in the container. You can use an AGAT's tool e.g. agat_convert_sp_gxf2gxf.pl doing
```
agat_convert_sp_gxf2gxf.pl --help
```
   </details>

## Using Bioconda
      
### Install AGAT

  ```
  conda install -c bioconda agat
  ```

or in a fresh environment:

  ```
  conda create -c bioconda -n agat agat
  ```

### Update AGAT

  ```
  conda update agat
  ```

### Uninstall AGAT
  ```
  conda uninstall agat  
  ```

   
## Old school - Manually

You will have to install all prerequisites and AGAT manually.

### Install prerequisites
  * R (optional)  
    You can install it by conda (`conda install r-base`), through [CRAN](https://cran.r-project.org) ([See here for a nice tutorial](https://www.datacamp.com/community/tutorials/installing-R-windows-mac-ubuntu)) or using your package management tool (e.g apt for Debian, Ubuntu, and related Linux distributions). R is optional and can be used to perform some plots. You will need to install the perl depency Statistics::R

  * Perl >= 5.8  
    It should already be available on your computer. If you are unlucky [perl.org](https://www.perl.org/get.html) is the place to go.

  * Perl modules  
    They can be installed in different ways:

    * using cpan or cpanm

    ```
    cpanm install bioperl Clone Graph::Directed LWP::UserAgent Carp Sort::Naturally File::Share File::ShareDir::Install Moose YAML LWP::Protocol::https Term::ProgressBar
    ```

    * using conda

      * using the provided yaml file

      ```
      conda env create -f conda_environment_AGAT.yml
      conda activate agat
      ```

      * manually  

      ```
      conda install perl-bioperl perl-clone perl-graph perl-lwp-simple perl-carp perl-sort-naturally perl-file-share perl-file-sharedir-install perl-moose perl-yaml perl-lwp-protocol-https perl-term-progressbar
      ```

    * using your package management tool (e.g apt for Debian, Ubuntu, and related Linux distributions)

    ```
    apt install libbio-perl-perl libclone-perl libgraph-perl liblwp-useragent-determined-perl libstatistics-r-perl libcarp-clan-perl libsort-naturally-perl libfile-share-perl libfile-sharedir libfile-sharedir-install-perl libyaml-perl liblwp-protocol-https-perl libterm-progressbar-perl
    ```

  * Optional
    Some scripts offer the possibility to perform plots. You will need R and Statistics::R which are not included by default.

    * R   
      You can install it by conda (`conda install r-base`), through [CRAN](https://cran.r-project.org) ([See here for a nice tutorial](https://www.datacamp.com/community/tutorials/installing-R-windows-mac-ubuntu)) or using your package management tool (e.g apt for Debian, Ubuntu, and related Linux distributions).

    * Statistics::R
        You can install it through conda  (`conda install perl-statistics-r`), using cpan/cpanm (`cpanm install Statistics::R`), or your package management tool  (`apt install libstatistics-r-perl`)



### Install AGAT

  ```
  git clone https://github.com/NBISweden/AGAT.git # Clone AGAT
  cd AGAT                                         # move into AGAT folder
  perl Makefile.PL                                # Check all the dependencies*
  make                                            # Compile
  make test                                       # Test
  make install                                    # Install
  ```

<sup>*</sup>If dependencies are missing you will be warn. Please refer to the [Install prerequisites](#install-prerequisites) section.

**Remark**: On MS Windows, instead of make you'd probably have to use dmake or nmake depending the toolchain you have.

### Update AGAT
From the folder where the repository is located.

  ```
  git pull                                        # Update to last AGAT
  perl Makefile.PL                                # Check all the dependencies*
  make                                            # Compile
  make test                                       # Test
  make install                                    # Install
  ```
<sup>*</sup>If dependencies are missing you will be warn. Please refer to the [Install prerequisites](#install-prerequisites) section.

### Change to a specific version
From the folder where the repository is located.  

  ```
  git pull                                        # Update the code
  git checkout v0.1                               # use version v0.1 (See releases tab for a list of available versions)
  perl Makefile.PL                                # Check all the dependencies*
  make                                            # Compile
  make test                                       # Test
  make install                                    # Install
  ```
<sup>*</sup>If dependencies are missing you will be warn. Please refer to the [Install prerequisites](#install-prerequisites) section.

### Uninstall AGAT

  ```
  perl uninstall_AGAT
  ```
   </details>
