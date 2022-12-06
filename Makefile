GCC:=g++
BUILD_DIR:=build
LIBS:=-Lincludes

CFLAGS:=-c -g -Wall 
CFLAGS+=$(shell root-config --cflags)\
		$(shell root-config --auxcflags)

LDFLAGS:=$(shell root-config --ldflags)
LIBS:=$(LIBS) $(shell root-config --libs)

SRC:=$(wildcard *.cc)
OBJ:=$(patsubst %.cc, $(BUILD_DIR)/%.o, $(SRC))

EXE:=$(patsubst %.cc, %, $(SRC))

MKDIR=[ -d $(@D) ] || mkdir -p $(@D)
all: $(EXE)
obj: $(OBJ)

$(EXE) : $(OBJ)
	$(MKDIR)
	$(GCC) $(LDFLAGS) $(patsubst %, $(BUILD_DIR)/%.o, $@) -o $@ $(LIBS)

$(OBJ) : $(SRC)
	$(MKDIR)
	$(GCC) $(CFLAGS) $(patsubst $(BUILD_DIR)/%.o, %.cc ,$@) -o $@

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR) *~
	rm -f $(EXE)
