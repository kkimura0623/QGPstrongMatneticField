 #!/bin/bash -ex

export PATH=$PATH:/usr/bin:/bin
export HOME=/mnt/disk1/kimura/QGPstrongMatneticField/photon_prop_v1.0_magstrng
workdir=/mnt/disk1/kimura/QGPstrongMatneticField/photon_prop_v1.0_magstrng/work/job$1
outdir=/mnt/disk1/kimura/QGPstrongMatneticField/photon_prop_v1.0_magstrng/out
mkdir -p $workdir

cp -r /mnt/disk1/kimura/QGPstrongMatneticField/photon_prop_v1.0_magstrng $workdir
cd $workdir
make
echo $1
./epem_pair_prod_test.F90 $1 1
./epem_pair_prod_test.F90 $1 2
cp *.dat $outdir
rm -rf $workdir
