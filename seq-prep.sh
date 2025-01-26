#!bin/bash
# run with sudo!

# minimum linux system requirements: bash curl wget git make sed runtest
# FOR x86-64 architecture only

set -e # exit on error
set -o pipefall # fail on piped commands

### Bioinformatics setup script to initialize naive linux system ###


# Temp folder set up
TEMP_DIR="seqtemp"
if [ ! -d $TEMP_DIR ]; then
    mkdir -p $TEMP_DIR
fi
cd $TEMP_DIR

# LANGUAGE SUPPORT
# verify python installation
PYTHON3_VERSION=3.12.8

if command -v python3 &>/dev/null; then
    PYTHON_VERSION=$(python3 -c 'import sys; print(".".join(map(str, sys.version_info[:3])))')
    if [[ $(echo -e "$PYTHON_VERSION\n3.10" | sort -V | head -n1) == "3.10" ]]; then
        echo "Python 3.10 or higher is already installed."
    else
    	# uninstall and rebuild python version
    	echo "Python version is lower than 3.10. Upgrading Python"

    	# remove
    	rm -rf /usr/local/bin/python3
    	rm -rf /usr/local/bin/python3*

    	#build
    	mkdir -p /usr/src
    	cd /usr/src
    	wget https://www.python.org/ftp/python/$PYTHON3_VERSION/Python-$PYTHON3_VERSION.tgz
    	tar xzf Python-$PYTHON3_VERSION.tgz
    	cd Python-$PYTHON3_VERSION
    	./configure --enable-optimizations
    	make -j$(nproc)
    	make altinstall
    fi
else
    # build python
    echo "Python is not installed. Installing Python $PYTHON3_VERSION"
    mkdir -p /usr/src
    cd /usr/src
    wget https://www.python.org/ftp/python/$PYTHON3_VERSION/Python-$PYTHON3_VERSION.tgz
    tar xzf Python-$PYTHON3_VERSION.tgz
    cd Python-$PYTHON3_VERSION
    ./configure --enable-optimizations
    make -j$(nproc)
    make altinstall
fi

# build anaconda3
ANACONDA_VERSION=latest

echo "Installing Miniconda"
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-${ANACONDA_VERSION}-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh 
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3
rm ~/miniconda3/miniconda.sh
source ~/miniconda3/bin/activate
conda init --all

# build R 

# build Java
JAVA_VERSION=21

# BIOINFORMATICS TOOLS

# build samtools
SAMTOOLS_VERSION=1.21

echo "Installing Samtools"
wget https://github.com/samtools/samtools/releases/download/${SAMTOOLS_VERSION}/samtools-${SAMTOOLS_VERSION}.tar.bz2
tar -xjf samtools-${SAMTOOLS_VERSION}.tar.bz2
cd samtools-${SAMTOOLS_VERSION}
./configure
make
make install
cd ..

# build bcftool
BCFTOOLS_VERSION=1.21

echo "Installing BCF tools"
wget https://github.com/samtools/bcftools/releases/download/${BCFTOOLS_VERSION}/bcftools-${BCFTOOLS_VERSION}.tar.bz2
tar -xjf bcftools-${BCFTOOLS_VERSION}.tar.bz2
cd bcftools-${BCFTOOLS_VERSION}
./configure
make
make install
cd ..

# build gatk
GATK_VERSION=4.6.1.0

echo "Installing GATK"
wget https://github.com/broadinstitute/gatk/releases/download/4.6.1.0/gatk-4.6.1.0.zip
unzip gatk-4.6.1.0.zip -d ~/gatk

# build bwa
echo "Installing BWA"
git clone https://github.com/lh3/bwa.git
cd bwa
make
cp bwa /usr/local/bin
cd ..

# build fastqc
FASTQC_VERSION=0.12.1

echo "Installing FASTQC"
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_{$FASTQC_VERSION}.zip
unzip fastqc_{$FASTQC_VERSION}.zip
chmod +x FastQC/fastqc
cp FastQC/fastqc /usr/local/bin


# build bowtie2
BOWTIE2_VERSION=2.5.4

echo "Installing Bowtie2"
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/{$BOWTIE2_VERSION}/bowtie2-{$BOWTIE2_VERSION}-linux-x86_64.zip/download
unzip bowtie2.zip -d ~/bowtie2
cp ~/bowtie2/bowtie2* /usr/local/bin


# Remove Temp
cd ..
rm -rf $TEMP_DIR
