#!/usr/bin/make -rRf
my_mkfile_path := $(abspath $(lastword $(MAKEFILE_LIST)))
my_current_dir :=  $(patsubst %/,%,$(dir $(my_mkfile_path)))/

CXX?=g++
CMAKE?=cmake
cc_files?= calib_cons.cc ../clustering/fastq_reader.cc
cc_flags?=
cmake_args?=-DCMAKE_BUILD_TYPE=Release -DCMAKE_POLICY_VERSION=3.25 -DCMAKE_POLICY_VERSION_MINIMUM=3.5
ZSTD_DIR?=../vendor/zstd
ZSTD_LIB?=$(ZSTD_DIR)/lib/libzstd.a
cc_args?= $(cc_flags)-std=c++11 -O3 -pthread -Ispoa_v1.1.3/include/ -I$(ZSTD_DIR)/lib -Lspoa_v1.1.3/build/lib/ -L$(ZSTD_DIR)/lib
cmake_major:=$(shell $(CMAKE) --version | head -n1 | awk '{print $$3}' | cut -d. -f1)
cmake_ok:=$(shell test "$(cmake_major)" -ge 3 >/dev/null 2>&1 && echo yes)
$(info Detected CMake major version $(cmake_major))
calib_cons: $(cc_files) $(ZSTD_LIB) spoa_v1.1.3/build/lib/libspoa.a 
	$(CXX) $(cc_files) $(cc_args) -lspoa $(LDLIBS) $(ZLIB_LIB) -o calib_cons

$(ZSTD_LIB):
	$(MAKE) -C .. zstd

spoa_v1.1.3/build/lib/libspoa.a:
ifeq ($(cmake_ok),yes)
	@echo Using $(CMAKE)
	@if [ ! -d spoa_v1.1.3 ]; then \
		git clone --recursive https://github.com/rvaser/spoa spoa_v1.1.3; \
	else \
		echo "spoa_v1.1.3 already present, skipping remote fetch"; \
		cd spoa_v1.1.3 && git submodule update --init --recursive || true; \
	fi
	cd spoa_v1.1.3; git checkout 1.1.3; git submodule update --init --recursive; mkdir -p build; cd build; $(CMAKE) $(cmake_args) ..;
	make -C spoa_v1.1.3/build
else
	$(error "cmake needs to be version 3 or newer")
endif
