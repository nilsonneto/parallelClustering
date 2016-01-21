CC  := gcc

CFLAGS  := -O3 -fopenmp
LDFLAGS := -pthread -lm

TARGET = kmeansp kmeans

SOURCE_FILES = main.c em.c kmeans.c kmeans-gen.h
SOURCE_K = kmeans.c
SOURCE_KP = kmeans-parallel.c kmeans-gen.h

ALL: $(TARGET)

parallelClustering: $(SOURCE_FILES)
	$(CC) $(CFLAGS) $(SOURCE_FILES) -o parallelClustering $(LDFLAGS)

kmeans: $(SOURCE_K)
	$(CC) $(CFLAGS) $(SOURCE_K) -o kmeans $(LDFLAGS)

kmeansp: $(SOURCE_KP)
	$(CC) $(CFLAGS) $(SOURCE_KP) -o kmeansp $(LDFLAGS)
	
clean:
	rm -f $(TARGET) *.o 
