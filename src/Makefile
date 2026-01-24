CC = gcc
CXX = g++
JPEGLI_DIR = ../third_party/jpegli
JPEGLI_BUILD_DIR = $(JPEGLI_DIR)/build
TURBO_DIR = ../third_party/jpegli/third_party/libjpeg-turbo
TURBO_BUILD_DIR = $(TURBO_DIR)/build
TURBO_LIB = $(TURBO_BUILD_DIR)/libjpeg.a
MOZJPEG_DIR = ../third_party/mozjpeg
MOZJPEG_BUILD_DIR = $(MOZJPEG_DIR)/build
MOZJPEG_LIB = $(MOZJPEG_BUILD_DIR)/libjpeg.a
BOOST_DIR = ../third_party/boost
BOOST_DIR_ABS = $(abspath $(BOOST_DIR))
BOOST_HEADERS_DIR = $(BOOST_DIR_ABS)/boost
BOOST_PO_SRC_DIR = $(BOOST_DIR_ABS)/libs/program_options/src
BOOST_PO_BUILD_DIR = .build/boost_program_options
BOOST_PO_SRCS = $(wildcard $(BOOST_PO_SRC_DIR)/*.cpp)
BOOST_PO_OBJS = $(patsubst $(BOOST_PO_SRC_DIR)/%.cpp,$(BOOST_PO_BUILD_DIR)/%.o,$(BOOST_PO_SRCS))
BOOST_PO_LIB = $(BOOST_PO_BUILD_DIR)/libboost_program_options.a
BOOST_INCLUDE_FLAGS = -I$(BOOST_DIR_ABS)
LDFLAGS = -L$(JPEGLI_BUILD_DIR)/lib -L$(JPEGLI_BUILD_DIR)/third_party/highway -L$(TURBO_BUILD_DIR)
CFLAGS = -I$(JPEGLI_DIR) -I$(JPEGLI_BUILD_DIR)/lib/include -I$(JPEGLI_BUILD_DIR)/lib/include/jpegli $(BOOST_INCLUDE_FLAGS)
CXXFLAGS = -I$(JPEGLI_DIR) -I$(JPEGLI_BUILD_DIR)/lib/include -I$(JPEGLI_BUILD_DIR)/lib/include/jpegli $(BOOST_INCLUDE_FLAGS)
TURBO_CXXFLAGS = -I$(TURBO_DIR) -I$(TURBO_BUILD_DIR)
MOZJPEG_CXXFLAGS = -I$(MOZJPEG_DIR) -I$(MOZJPEG_BUILD_DIR)
UCFLAGS = -Ofast -std=gnu11 -fsigned-char $(CFLAGS)
UCXXFLAGS = -pthread -std=gnu++20 -fsigned-char $(CXXFLAGS)
PREFIX ?= /usr/local
BINPREFIX ?= ..
BINDIR ?= $(PREFIX)/bin
INSTALL ?= install
MAKEDIR ?= install -d
CMAKE = cmake
NPROC = $(shell nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 4)

ifeq ($(OS),Windows_NT)
	CXXFLAGS += -mno-ms-bitfields
	CMAKE += -G "MSYS Makefiles"
endif
#OBJECTS = transupp.o
CXXSRC = jpegultrascan.cpp turbo_encode.cpp
OBJS = .build/obj/jpegultrascan.o .build/obj/turbo_encode.o
OBJS_MOZ = .build/obj_moz/jpegultrascan.o .build/obj_moz/mozjpeg_encode.o

.PHONY: deps bin all clean jpegli boost_program_options boost_headers turbo mozjpeg mozjpeg-bin
all: deps bin

jpegli:
	@if [ ! -f "$(JPEGLI_BUILD_DIR)/lib/include/jpegli/jpeglib.h" ] || [ ! -f "$(JPEGLI_BUILD_DIR)/lib/libjpegli-static.a" ]; then \
		echo "Building jpegli..."; \
		mkdir -p $(JPEGLI_BUILD_DIR); \
		cd $(JPEGLI_BUILD_DIR) && \
		$(CMAKE) -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF -DJPEGXL_ENABLE_TOOLS=OFF -DJPEGXL_ENABLE_EXAMPLES=OFF -DJPEGXL_ENABLE_SKCMS=OFF -DCMAKE_POLICY_DEFAULT_CMP0007=NEW -DCMAKE_POLICY_VERSION_MINIMUM=3.5 .. || exit 1; \
		cd - > /dev/null; \
		if [ ! -f "$(JPEGLI_BUILD_DIR)/lib/include/jpegli/jpeglib.h" ]; then \
			echo "Error: jpeglib.h was not generated during CMake configuration."; \
			echo "Expected at: $(JPEGLI_BUILD_DIR)/lib/include/jpegli/jpeglib.h"; \
			exit 1; \
		fi; \
		cd $(JPEGLI_BUILD_DIR) && \
		$(CMAKE) --build . --target jpegli-static -- -j$(NPROC) || exit 1; \
		$(CMAKE) --build . --target hwy -- -j$(NPROC) || exit 1; \
		cd - > /dev/null; \
		if [ ! -f "$(JPEGLI_BUILD_DIR)/lib/libjpegli-static.a" ]; then \
			echo "Error: libjpegli-static.a was not built."; \
			exit 1; \
		fi; \
	fi

boost_program_options:
	@if [ ! -f "$(BOOST_PO_LIB)" ]; then \
		echo "Building boost program_options..."; \
		mkdir -p $(BOOST_PO_BUILD_DIR); \
		for src in $(BOOST_PO_SRCS); do \
			obj=$(BOOST_PO_BUILD_DIR)/$$(basename $$src .cpp).o; \
			$(CXX) $(UCXXFLAGS) $(BOOST_INCLUDE_FLAGS) -c $$src -o $$obj || exit 1; \
		done; \
		ar rcs $(BOOST_PO_LIB) $(BOOST_PO_OBJS) || exit 1; \
	fi

deps: jpegli boost_program_options turbo

turbo:
	@if [ ! -f "$(TURBO_LIB)" ]; then \
		echo "Building libjpeg-turbo..."; \
		mkdir -p $(TURBO_BUILD_DIR); \
		cd $(TURBO_BUILD_DIR) && \
		$(CMAKE) -DCMAKE_BUILD_TYPE=Release -DENABLE_SHARED=OFF -DENABLE_STATIC=ON -DWITH_TURBOJPEG=OFF -DWITH_SIMD=OFF -DCMAKE_POLICY_VERSION_MINIMUM=3.5 .. || exit 1; \
		cd - > /dev/null; \
		cd $(TURBO_BUILD_DIR) && \
		$(CMAKE) --build . --target jpeg-static -- -j$(NPROC) || exit 1; \
		cd - > /dev/null; \
		if [ ! -f "$(TURBO_LIB)" ]; then \
			echo "Error: libjpeg.a was not built."; \
			exit 1; \
		fi; \
	fi

mozjpeg:
	@if [ ! -f "$(MOZJPEG_LIB)" ]; then \
		echo "Building mozjpeg..."; \
		mkdir -p $(MOZJPEG_BUILD_DIR); \
		cd $(MOZJPEG_BUILD_DIR) && \
		$(CMAKE) -DCMAKE_BUILD_TYPE=Release -DBUILD_SHARED_LIBS=OFF -DWITH_TURBOJPEG=OFF -DWITH_SIMD=OFF -DCMAKE_POLICY_VERSION_MINIMUM=3.5 .. || exit 1; \
		cd - > /dev/null; \
		cd $(MOZJPEG_BUILD_DIR) && \
		$(CMAKE) --build . --target jpeg-static -- -j$(NPROC) || exit 1; \
		cd - > /dev/null; \
		if [ ! -f "$(MOZJPEG_LIB)" ]; then \
			echo "Error: mozjpeg libjpeg.a was not built."; \
			exit 1; \
		fi; \
	fi

bin: deps $(OBJS)
	@if [ ! -f "$(JPEGLI_BUILD_DIR)/lib/include/jpegli/jpeglib.h" ]; then \
		echo "Error: jpeglib.h not found. Run 'make clean-all' and 'make' again."; \
		exit 1; \
	fi
	@if [ -f "$(JPEGLI_BUILD_DIR)/lib/libjpegli-static.a" ]; then \
		$(CXX) $(UCXXFLAGS) $(OBJECTS) $(OBJS) -o ${BINPREFIX}/jpegultrascan $(BOOST_PO_LIB) $(LDFLAGS) -ljpegli-static -ljpeg -lhwy -lpthread -g; \
	elif [ -f "$(JPEGLI_BUILD_DIR)/lib/libjpegli.a" ]; then \
		$(CXX) $(UCXXFLAGS) $(OBJECTS) $(OBJS) -o ${BINPREFIX}/jpegultrascan $(BOOST_PO_LIB) $(LDFLAGS) -ljpegli -ljpeg -lhwy -lpthread -g; \
	else \
		echo "Error: jpegli library not found. Run 'make jpegli' first."; \
		exit 1; \
	fi

.build/obj:
	mkdir -p .build/obj

.build/obj_moz:
	mkdir -p .build/obj_moz

.build/obj/jpegultrascan.o: jpegultrascan.cpp | .build/obj
	$(CXX) $(UCXXFLAGS) -DENABLE_TURBO=1 -DENABLE_MOZJPEG=0 -c $< -o $@

.build/obj/turbo_encode.o: turbo_encode.cpp | .build/obj
	$(CXX) $(UCXXFLAGS) $(TURBO_CXXFLAGS) -DENABLE_TURBO=1 -DENABLE_MOZJPEG=0 -c $< -o $@

.build/obj_moz/jpegultrascan.o: jpegultrascan.cpp | .build/obj_moz
	$(CXX) $(MOZJPEG_CXXFLAGS) $(UCXXFLAGS) -DENABLE_TURBO=0 -DENABLE_MOZJPEG=1 -DUSE_MOZJPEG_API=1 -c $< -o $@

.build/obj_moz/mozjpeg_encode.o: mozjpeg_encode.cpp | .build/obj_moz
	$(CXX) $(MOZJPEG_CXXFLAGS) $(UCXXFLAGS) -DENABLE_TURBO=0 -DENABLE_MOZJPEG=1 -DUSE_MOZJPEG_API=1 -c $< -o $@

mozjpeg-bin: mozjpeg boost_program_options $(OBJS_MOZ)
	@if [ -f "$(JPEGLI_BUILD_DIR)/lib/libjpegli-static.a" ]; then \
		:; \
	fi
	$(CXX) $(UCXXFLAGS) $(OBJECTS) $(OBJS_MOZ) -o ${BINPREFIX}/jpegultrascan_mozjpeg $(BOOST_PO_LIB) $(MOZJPEG_LIB) -lpthread -g

clean:
	rm -f *.o *.a
	rm -rf .build

clean-all: clean
	rm -rf $(JPEGLI_BUILD_DIR)
