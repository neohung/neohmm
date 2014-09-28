CC = gcc
INCLUDE = .
CFLAGS =
OBJ = main.o libneohmm.a
TARGET = main

all: $(OBJ)
	gcc -o $(TARGET) main.o -L. -lneohmm -lm

%.o: %.c
	echo *** [GCC] $@ : $<
	$(CC) -c -o $@ $< $(CFLAGS) -I $(INCLUDE) 

%.a: %.o
	echo *** [AR] $@ : $<
	ar rcs $@ $<

.PHONY: clean
clean:
	rm -f $(TARGET)
	rm -f $(OBJ)
run:
	./$(TARGET)

