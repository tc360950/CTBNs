PROG = CTBN
CC = g++
CPPFLAGS = -O3 -std=c++14 -Wall -Wextra -lpthread -g
LDFLAGS = -lpthread
OBJS = CTBN.o 

$(PROG) : $(OBJS)
	$(CC) $(LDFLAGS)  -o $(PROG) $(OBJS) -lpthread
CNVTree.o:
	$(CC) $(CPPFLAGS) -c CNVTree.cpp
clean:
	rm -f core $(PROG) $(OBJS)
