CFLAGS = -Wall -Wextra -mkl -O3
CC = mpiicc
FUNC_DIR = func
funcObj = $(FUNC_DIR)/calcBlock.o $(FUNC_DIR)/matrix.o $(FUNC_DIR)/matrixIO.o $(FUNC_DIR)/pivot.o $(FUNC_DIR)/luFunc.o $(FUNC_DIR)/randomGen.o
tsluObj = $(FUNC_DIR)/tslu.o
caluObj = $(FUNC_DIR)/calu.o
all: calu


$(FUNC_DIR)/tslu.o: $(objects) $(FUNC_DIR)/tslu.c
	$(CC) $(CFLAGS) -c $^
	mv tslu.o func
calu: $(tsluObj) $(funcObj) $(caluObj) main.c
	$(CC) $(CFLAGS) -o $@ $^

random: $(FUNC_DIR)/randomGen.o $(FUNC_DIR)/matrixIO.o random.c
	$(CC) $(CFLAGS) -o $@ $^

test: $(tsluObj) $(funcObj) test.c
	$(CC) $(CFLAGS) -o $@ $^

serial: $(tsluObj) $(funcObj) serial.c
	$(CC) $(CFLAGS) -o $@ $^

.PHONY: all clean

clean:
	rm -f *.o pivot random $(FUNC_DIR)/*.o

