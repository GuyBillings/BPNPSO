CC     = mpicc
CFLAGS = -o3 -n -c
LINK   = -lmpfr -lgmp -lgsl -lmpi
 
all: bpnpso
 
bpnpso: bpnpso.o bpnpsofuncs.o swarm.o
	$(CC) $(LINK) bpnpso.o bpnpsofuncs.o swarm.o -o bpnpso

swarm.o: swarm.c
	$(CC) $(CFLAGS) $(LINK) swarm.c -o swarm.o

bpnpsofuncs.o: bpnpsofuncs.c
	$(CC) $(CFLAGS) $(LINK) bpnpsofuncs.c -o bpnpsofuncs.o
 
bpnpso.o: bpnpso.c
	$(CC) $(CFLAGS) $(LINK) bpnpso.c -o bpnpso.o
 
clean:
	rm -f bpnpso bpnpso.o bpnpsofuncs.o swarm.o
 
  
