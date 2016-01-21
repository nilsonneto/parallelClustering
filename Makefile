CC  := gcc

CFLAGS  := -O3
LDFLAGS := -pthread -lm

TARGET = parallelClustering kmeans

SOURCE_FILES = main.c em.c kmeans.c
SOURCE_K = kmeans.c

ALL: $(TARGET)

parallelClustering: $(SOURCE_FILES)
	$(CC) $(CFLAGS) $(SOURCE_FILES) -o parallelClustering $(LDFLAGS)

kmeans: $(SOURCE_K)
	$(CC) $(CFLAGS) $(SOURCE_K) -o kmeans $(LDFLAGS)
	
clean:
	rm -f $(TARGET) *.o 
