GCC:=g++
BUILD_DIR:=build
LIBS_DIR:=include

CFLAGS:=-c -g -Wall 
CFLAGS+=$(shell root-config --cflags)

LDFLAGS:=$(shell root-config --libs)

SRC:=$(wildcard *.cc)
OBJ:=$(patsubst %.cc, $(BUILD_DIR)/%.o, $(SRC))

EXE:=clusterise

MKDIR=[ -d $(@D) ] || mkdir -p $(@D)

all: $(EXE)

$(EXE) : $(OBJ) $(OBJ_LIBS)
	$(MKDIR)
	$(GCC) $(LDFLAGS) $< -o $@

$(OBJ) : $(SRC)
	$(MKDIR)
	$(GCC) $(CFLAGS) $< -o $@

.PHONY: clean
clean:
	rm -rf $(BUILD_DIR) *~
	rm -f $(EXE)
