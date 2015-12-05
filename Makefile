CFLAGS  = -L/usr/5lib
LFLAGS  = -lm
CC      = gcc 


pol3scan:	pol3scan.c
		$(CC) $(CFLAGS)  pol3scan.c $(LFLAGS)  -o pol3scan

