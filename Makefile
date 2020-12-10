
CC=gcc
CILKCC=/usr/local/OpenCilk-9.0.1-Linux/bin/clang
CFLAGS=-w -O3 -fopenmp -pthread
CILKFLAGS=-w -O3 -fcilkplus
BINFILES = v3_2 v3_3 v3_3merge v3_3mergeOmpOutDynamic v4 v4pthreads v4ompInStatic v4ompOutDynamic v3_3mergeOmpInStatic v3_3mergeCilk v4cilk
EXTRALIBS = mmio.o coo2csc.o

all: $(BINFILES)

%.o: %.c 
    $(CC) $(CFLAGS) -c $^

v3_2: v3_2.o $(EXTRALIBS)
    $(CC) $(CFLAGS) -o $@ $^

v3_3: v3_3.o $(EXTRALIBS)
    $(CC) $(CFLAGS) -o $@ $^

v3_3merge: v3_3merge.o $(EXTRALIBS)
    $(CC) $(CFLAGS) -o $@ $^

v3_3mergeOmpOutDynamic: v3_3mergeOmpOutDynamic.o $(EXTRALIBS)
    $(CC) $(CFLAGS) -o $@ $^

v4: v4.o $(EXTRALIBS)
    $(CC) $(CFLAGS) -o $@ $^

v4pthreads: v4pthreads.o $(EXTRALIBS)
    $(CC) $(CFLAGS) -o $@ $^

v3_3mergeOmpInStatic: v3_3mergeOmpInStatic.o $(EXTRALIBS)
    $(CC) $(CFLAGS) -o $@ $^

v4ompInStatic: v4ompInStatic.o $(EXTRALIBS)
    $(CC) $(CFLAGS) -o $@ $^

v4ompOutDynamic: v4ompOutDynamic.o $(EXTRALIBS)
    $(CC) $(CFLAGS) -o $@ $^

v3_3mergeCilk:
    $(CILKCC) $(CILKFLAGS) -o v3_3mergeCilk v3_3mergeCilk.c mmio.c coo2csc.c

v4cilk:
    $(CILKCC) $(CILKFLAGS) -o v4cilk v4cilk.c mmio.c coo2csc.c

clean:
    rm -f $(BINFILES) *.o
