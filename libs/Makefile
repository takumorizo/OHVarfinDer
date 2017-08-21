default: boost_1_45_0_subset/boost samtools-0.1.19/samtools

boost_1_45_0_subset/boost:
	echo "preparing required boost headers ..."
	tar zxvf boost_1_45_0_subset.tar.gz

samtools-0.1.19/samtools:
	echo "preparing samtools v0.1.19 ..."
	tar zxvf samtools-0.1.19.tar.gz
	cd samtools-0.1.19; make
	cd ..
