TBB_INSTALL_DIR := $(shell dirname -- "$0";)/lib/tbb
TBB_INCLUDE := ${TBB_INSTALL_DIR}/include
TBB_LIBRARY_RELEASE := ${TBB_INSTALL_DIR}/build/linux_intel64_gcc_cc7.5.0_libc2.27_kernel4.15.0_release
TBB_LIBRARY_DEBUG := ${TBB_INSTALL_DIR}/build/linux_intel64_gcc_cc7.5.0_libc2.27_kernel4.15.0_debug

all: consult_classify

clean:
	@echo "Succesfully cleaned."
	rm -f ./consult_classify

consult_classify:
	g++ consult_classify.cpp -o consult_classify -I${TBB_INCLUDE} -Wl,-rpath,${TBB_LIBRARY_RELEASE} -L${TBB_LIBRARY_RELEASE} -ltbb -Wfatal-errors
