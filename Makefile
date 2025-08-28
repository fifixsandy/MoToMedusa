SRC_DIR:=src
OBJ_DIR:=obj
BIN_DIR:=.
LIB_DIR:=lib
LACE_DIR:=$(LIB_DIR)/sylvan/build/_deps
BUDDY_DIR:=$(LIB_DIR)/buddy

SRCS:=$(wildcard $(SRC_DIR)/*.c)
EXEC:=$(BIN_DIR)/MEDUSA
OBJS:=$(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRCS))

INTERFACE_DIR := $(SRC_DIR)/interface_implementations
INTERFACE_OBJ_default := $(OBJ_DIR)/interface_sylvan_gmp.o

INTERFACE_OBJ_sylvan_gmp := $(OBJ_DIR)/interface_sylvan_gmp.o
INTERFACE_OBJ_sylvan_doubles := $(OBJ_DIR)/interface_sylvan_doubles.o
INTERFACE_OBJ_buddy_gmp := $(OBJ_DIR)/interface_buddy_gmp.o
INTERFACE_OBJ_buddy_doubles := $(OBJ_DIR)/interface_buddy_doubles.o

CC:=gcc
CFLAGS:=-g -O2
CLIBS=-lgmp -lpthread -lm
INC_DIRS:=-I $(LIB_DIR)/sylvan/src/ -I $(LACE_DIR)/lace-src/src/ -I $(LACE_DIR)/lace-build/
INC_DIRS += -I $(BUDDY_DIR)/src/

N_JOBS=4

OF_TYPE=pdf
F_OUT_NAME=res
LONG_NUMS_OUT_FILE=res-vars.txt
BSCRIPT_PATH=benchmark-utils/scripts

.DEFAULT : all
.PHONY : clean clean-all clean-artifacts clean-deps clean-benchmark plot benchmarks \
           init make-sylvan download-sylvan make-sliqsim sylvan_gmp buddy_gmp \
		   make-motobuddy download-motobuddy init-motobuddy
all: $(OBJS) $(INTERFACE_OBJ_default) $(LIB_DIR)/sylvan/build/src/lib/libsylvan.a $(LACE_DIR)/lace-build/lib/liblace.a | $(BIN_DIR)
	$(CC) $(INC_DIRS) $(CFLAGS) -o $(EXEC) $^ $(CLIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(INC_DIRS) $(CFLAGS) -c $< -o $@

$(BIN_DIR) $(OBJ_DIR):
	mkdir -p $@

-include $(OBJ:.o=.d)

plot:
	@dot -T$(OF_TYPE) $(F_OUT_NAME).dot -o $(F_OUT_NAME).$(OF_TYPE)

benchmarks:
	@bash ./$(BSCRIPT_PATH)/run-benchmarks.sh

# BENCHMARK INIT:
make-sliqsim:
	cd .. &&\
	git clone https://github.com/NTU-ALComLab/SliQSim.git || true &&\
	cd SliQSim/cudd &&\
	./configure --enable-dddmp --enable-obj --enable-shared --enable-static &&\
	cd .. &&\
	make

# INIT:
init: make-sylvan
	mkdir -p $(LIB_DIR) && mv sylvan $(LIB_DIR)

init-motobuddy: make-motobuddy
	mkdir -p $(LIB_DIR) && mv motobuddy $(LIB_DIR)

make-sylvan: download-sylvan
	cd sylvan;			\
	mkdir build;		\
	cd build;			\
	cmake ..;			\
	make -j $(N_JOBS);

make-motobuddy: download-motobuddy
	cd MoToBuddy; \
	mkdir build; \
	cd build; \
	cmake ..;\
	make;

download-sylvan:
	@git clone https://github.com/trolando/sylvan.git || true

download-motobuddy:
	@git clone https://github.com/VeriFIT/MoToBuddy.git || true

# CLEAN:
clean: clean-artifacts

clean-all: clean-artifacts clean-deps clean-benchmark

clean-artifacts:
	rm -rf $(EXEC) $(F_OUT_NAME).dot $(F_OUT_NAME).$(OF_TYPE) $(LONG_NUMS_OUT_FILE) $(OBJ_DIR)

clean-deps:
	rm -rf $(LIB_DIR)

clean-benchmark:
	cd .. && rm -rf SliQSim

$(OBJ_DIR)/%.o: $(INTERFACE_DIR)/%.c | $(OBJ_DIR)
	$(CC) $(INC_DIRS) $(CFLAGS) -c $< -o $@

# sylvan_gmp target
sylvan_gmp: $(OBJS) $(INTERFACE_OBJ_sylvan_gmp) $(LIB_DIR)/sylvan/build/src/lib/libsylvan.a $(LACE_DIR)/lace-build/lib/liblace.a | $(BIN_DIR)
	$(CC) $(INC_DIRS) $(CFLAGS) -o $(BIN_DIR)/MEDUSA_sylvan_gmp $^ $(CLIBS)

sylvan_doubles:  $(OBJS) $(INTERFACE_OBJ_sylvan_doubles) $(LIB_DIR)/sylvan/build/src/lib/libsylvan.a $(LACE_DIR)/lace-build/lib/liblace.a | $(BIN_DIR)
	$(CC) $(INC_DIRS) $(CFLAGS) -o $(BIN_DIR)/MEDUSA_sylvan_doubles $^ $(CLIBS)

# todo check if libsylvan and liblace needed
buddy_gmp: $(OBJS) $(INTERFACE_OBJ_buddy_gmp) $(LIB_DIR)/sylvan/build/src/lib/libsylvan.a $(LACE_DIR)/lace-build/lib/liblace.a $(LIB_DIR)/buddy/src/libbuddy.a | $(BIN_DIR)
	$(CC) $(INC_DIRS) $(CFLAGS) -o $(BIN_DIR)/MEDUSA_buddy_gmp $^ $(CLIBS)

buddy_doubles: $(OBJS) $(INTERFACE_OBJ_buddy_doubles) $(LIB_DIR)/sylvan/build/src/lib/libsylvan.a $(LACE_DIR)/lace-build/lib/liblace.a $(LIB_DIR)/buddy/src/libbuddy.a | $(BIN_DIR)
	$(CC) $(INC_DIRS) $(CFLAGS) -o $(BIN_DIR)/MEDUSA_buddy_doubles $^ $(CLIBS)