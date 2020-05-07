#!/bin/bash

curl -O https://bbp.epfl.ch/nmc-portal/documents/10184/7288948/hoc_combos_syn.1_0_10.allzips.tar/cef96ba5-35ea-45a9-a2ef-2f27ec12ae6b
tar -xvf cef96ba5-35ea-45a9-a2ef-2f27ec12ae6b
rm -rf cef96ba5-35ea-45a9-a2ef-2f27ec12ae6b
mv hoc_combos_syn.1_0_10.allzips AllLayers

cd AllLayers
mkdir L1 L23 L4 L5 L6
 
mv L1* L1/
for f in L1/*
do
	echo "unzip L1 - $f"
	unzip $f -d L1/
done
rm -rf L1/*zip


mv L23* L23/ 
for f in L23/*
do
	echo "unzip L23 - $f"
	unzip $f -d L23/
done
rm -rf L23/*zip



mv L4* L4/ 
for f in L4/*
do
	echo "unzip L4 - $f"
	unzip $f -d L4/
done
rm -rf L4/*zip


mv L5* L5/ 
for f in L5/*
do
	echo "unzip L5 - $f"
	unzip $f -d L5/
done
rm -rf L5/*zip


mv L6* L6/ 
for f in L6/*
do
	echo "unzip L6 - $f"
	unzip $f -d L6/
done
rm -rf L6/*zip
cd ../