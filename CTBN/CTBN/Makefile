PROG = CTBN
CC = g++
CPPFLAGS = -O3 -std=c++14 -Wall -Wextra -lpthread -g
LDFLAGS = -lpthread
OBJS = CTBN.o parameters.o

$(PROG) : $(OBJS)
	$(CC) $(LDFLAGS)  -o $(PROG) $(OBJS) -lpthread
CTBN.o:
	$(CC) $(CPPFLAGS) -c CTBN.cpp
parameters.o: utils/parameters.h
	$(CC) $(CPPFLAGS) -c utils/parameters.cpp 
clean:
	rm -f core $(PROG) $(OBJS)