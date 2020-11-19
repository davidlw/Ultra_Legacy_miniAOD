ROOT=`root-config --cflags --glibs`
CXX=g++
CXXFLAGS=-Wall -O2 -Wextra -std=c++11
ifeq "$(GCCVERSION)" "1"
  CXXFLAGS += -Wno-error=misleading-indentation
endif

MKDIR_BIN=mkdir -p $(PWD)/bin

all: mkdirBin bin/MyClass.exe bin/genOnlyAnalyzer.exe bin/HBT.exe bin/Cumulants.exe

mkdirBin:
	$(MKDIR_BIN)

bin/genOnlyAnalyzer.exe: src/genOnlyAnalyzer.C
	$(CXX) $(CXXFLAGS) -Wno-error=maybe-uninitialized $(ROOT) -I $(PWD) -o bin/genOnlyAnalyzer.exe src/genOnlyAnalyzer.C
bin/HBT.exe: src/HBT.C
	$(CXX) $(CXXFLAGS) -Wno-error=maybe-uninitialized $(ROOT) -I $(PWD) -o bin/HBT.exe src/HBT.C
bin/Cumulants.exe: src/Cumulants.C
	$(CXX) $(CXXFLAGS) -Wno-error=maybe-uninitialized $(ROOT) -I $(PWD) -o bin/Cumulants.exe src/Cumulants.C
bin/MyClass.exe: src/MyClass.C
	$(CXX) $(CXXFLAGS) -Wno-error=maybe-uninitialized $(ROOT) -I $(PWD) -o bin/MyClass.exe src/MyClass.C

clean:
	rm -f $(PWD)/include/*~
	rm -f $(PWD)/src/*~
	rm -f $(PWD)/src/*.so
	rm -f $(PWD)/src/*.d
	rm -f $(PWD)/src/*.pcm
	rm -f $(PWD)/bin/*.exe
	rmdir bin
