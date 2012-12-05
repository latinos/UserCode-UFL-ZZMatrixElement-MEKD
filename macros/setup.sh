export MEKD_COMPILE_WITH_ROOT=No

if [ -f ${ROOTSYS}/bin/root-config ]; then
	export MEKD_COMPILE_WITH_ROOT=Yes
fi

make clean
make
. setLocalLibrary.sh
