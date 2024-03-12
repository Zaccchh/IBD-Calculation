set -eu

if ! [ -f "sample/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz" ]; then
    echo "Downloading sample VCF file..."
    mkdir -p sample
    wget https://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/phase3/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz -O sample/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
fi

if ! [ -x "$(command -v bcftools)" ]; then
    echo "Downloading bcftools..."
    wget https://github.com/samtools/bcftools/releases/download/1.19/bcftools-1.19.tar.bz2
    tar -xjvf bcftools-1.19.tar.bz2
    new_dir="$PWD/bcftools"
    cd bcftools-1.19
    ./configure --prefix=$new_dir
    make
    make install
    cd ..
    rm -rf bcftools-1.19
    
    export PATH=$new_dir:$PATH
fi

echo "Indexing VCF file..."
bcftools index -t "sample/ALL.chr15.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"

echo "Installing requirements..."
pip install --user -r requirements.txt
