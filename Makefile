# ceigen_lite — builds on Linux, macOS, Windows (MSYS2/MinGW), x86_64/arm64
# Run natively on each target platform (e.g. a CI matrix), not cross-compiled.

UNAME_S := $(shell uname -s)
UNAME_M := $(shell uname -m)

# --- OS detection ---
ifneq (,$(findstring MINGW,$(UNAME_S)))
    OS := windows
else ifneq (,$(findstring MSYS,$(UNAME_S)))
    OS := windows
else ifneq (,$(findstring CYGWIN,$(UNAME_S)))
    OS := windows
else ifeq ($(UNAME_S),Darwin)
    OS := mac
else ifeq ($(UNAME_S),Linux)
    OS := linux
else
    $(error Unsupported OS: $(UNAME_S))
endif

# --- Arch normalization (uname spells arm64 differently per OS) ---
ifneq (,$(filter x86_64 amd64 AMD64,$(UNAME_M)))
    ARCH := x86_64
else ifneq (,$(filter aarch64 arm64 ARM64,$(UNAME_M)))
    ARCH := arm64
else
    $(error Unsupported architecture: $(UNAME_M))
endif

# --- Output naming ---
ifeq ($(OS),windows)
    EXT := dll
else ifeq ($(OS),mac)
    EXT := dylib
else
    EXT := so
endif

TARGET := libceigen_lite-$(ARCH)-$(OS).$(EXT)

# On Mac, g++ is clang
CXX      ?= g++
CXXFLAGS := --std=c++11 -O3 -fpic -I./ -Wno-enum-compare
LDFLAGS  := -shared

# --- Vectorization: AVX2 is x86-only, use it there; NEON baseline on arm64 ---
ifeq ($(ARCH),x86_64)
    CXXFLAGS += -mavx2
else
    CXXFLAGS += -march=armv8-a
endif

# --- OpenMP: Apple clang has no built-in OpenMP; needs Homebrew's libomp ---
ifeq ($(OS),mac)
    OMP_PREFIX := $(shell brew --prefix libomp 2>/dev/null)
    ifeq ($(OMP_PREFIX),)
        $(warning libomp not found via brew — building WITHOUT OpenMP. Run: brew install libomp)
        OMPFLAGS :=
    else
        OMPFLAGS := -Xpreprocessor -fopenmp -I$(OMP_PREFIX)/include -L$(OMP_PREFIX)/lib -lomp
    endif
else
    OMPFLAGS := -fopenmp
endif

.PHONY: all clean
all: $(TARGET)

$(TARGET): ceigen_lite.cpp
	$(CXX) $(CXXFLAGS) $(OMPFLAGS) $(LDFLAGS) -o $@ $<

clean:
	rm -f libceigen_lite-*.so libceigen_lite-*.dylib libceigen_lite-*.dll
