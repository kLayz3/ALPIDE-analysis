BUILD_DIR:=build
LIBS_DIR:=include

CFLAGS:=-c -g -Wall 
CFLAGS+=$(shell root-config --cflags)

LDFLAGS:=$(shell root-config --libs)

SRC:=$(wildcard *.cpp)
OBJ:=$(patsubst %.cpp, $(BUILD_DIR)/%.o, $(SRC))

EXE:=analyse

MKDIR=[ -d $(@D) ] || mkdir -p $(@D)

all: $(EXE)

$(EXE) : $(OBJ) $(OBJ_LIBS)
	$(MKDIR)
	g++ $(LDFLAGS) $< -o $@

$(OBJ) : $(SRC)
	$(MKDIR)
	g++ $(CFLAGS) $< -o $@

clean:
	rm -rf $(BUILD_DIR) *~
	rm -f $(EXE)
