CFILES = $(wildcard s21*.c)
OBJ = $(patsubst %.c, %.o, $(CFILES))
TARGET = s21_math
CHECK:=$(shell pkg-config --cflags --libs check)
FLAG_COV = --coverage
CC=gcc
CFLAGS=-Wall -Werror -Wextra -std=c11
S21_MATH=-L. s21_math.a
DEBUG =-g
CTESTFLAGS=-lcheck -lm -lpthread 
TEST = tests_math.c
LIBRARY = s21_math.a


all: clean $(LIBRARY) test gcov_report

$(TARGET): s21_math.o
	$(CC) s21_math.o -o $(TARGET)

s21_math.o : s21_math.c
	$(CC) $(DEBUG) -c s21_math.c

$(LIBRARY): s21_math.h $(CFILES)
	$(CC) $(CFLAGS) -c $(CFILES)
	ar rc $(LIBRARY) $(OBJ)
	ranlib $(LIBRARY)
	rm -f $(OBJ)

clean:
	rm -rf *.o *.a test *.gcda *.gcno *.info coverage

style:
	clang-format -style=google -i *.c
	clang-format -style=google -i *.h
	clang-format -style=google -n *.c
	clang-format -style=google -n *.h

checkstyle:
	clang-format -style=google -n *.c
	clang-format -style=google -n *.h

de: $(TARGET)
		
gcov_report: $(LIBRARY)
	$(CC) $(FLAG_COV) $(CFLAGS) $(CTESTFLAGS) $(TEST) -o test $(CHECK) $(S21_MATH)
	./test
	lcov -t "test" -o test.info -c -d .
	genhtml -o coverage test.info
	open coverage/index.html 

test: $(LIBRARY)
#	checkmk s21_math.check > tests_math.c
#	checkmk clean_mode=1 s21_math.check > tests_math.c
	$(CC) $(CTESTFLAGS) $(TEST) $(LIBRARY) -o test
	./test

rebuild:
	make clean
	make
