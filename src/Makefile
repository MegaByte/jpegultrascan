CC = gcc
CXX = g++
UCFLAGS = -Ofast -std=gnu11 -fsigned-char $(CFLAGS)
UCXXFLAGS = -pthread -std=gnu++20 -fsigned-char $(CXXFLAGS)
PREFIX ?= /usr/local
BINPREFIX ?= ..
BINDIR ?= $(PREFIX)/bin
INSTALL ?= install
MAKEDIR ?= install -d
CMAKE = cmake

ifeq ($(OS),Windows_NT)
	CXXFLAGS += -mno-ms-bitfields
	CMAKE += -G "MSYS Makefiles"
endif
OBJECTS = transupp.o
CXXSRC = jpegultrascan.cpp

.PHONY: mozjpeg deps bin all
all: deps bin

bin: deps
	$(CC) -c $(UCFLAGS) mozjpeg/transupp.c
	$(CXX) $(UCXXFLAGS) $(OBJECTS) $(CXXSRC) mozjpeg/libjpeg.a -o ${BINPREFIX}/jpegultrascan $(LDFLAGS) -lboost_program_options -g
clean:
	rm -f *.o *.a
	make -C mozjpeg clean
deps: mozjpeg
mozjpeg:
	cd mozjpeg/; \
	export CC="$(CC)"; \
	$(CMAKE) -DPNG_SUPPORTED=0 -DWITH_ARITH_DEC=1 -DWITH_ARITH_ENC=1 -DWITH_12BIT=1 -DENABLE_SHARED=0 .; \
	make jpegtran