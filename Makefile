CC=gcc
#CFLAGS= -c -O -march=native -Wall
#CFLAGS= -c -O3 -march=broadwell -mavx2  -Wall
CFLAGS= -c -O3 -march=broadwell -mavx2  -Rpass=loop-vectorize -Rpass-analysis=loop-vectorize -Wall
SOURCE=strassen.c verify.c 
LDFLAGS=
OBJECTS= $(SOURCE:.c=.o)

EXECUTABLE=verify


all:  $(SOURCE) $(EXECUTABLE)


$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.c.o:
	$(CC) $(CFLAGS) $< -o $@


clean: 
	rm -f *.o $(EXECUTABLE)
