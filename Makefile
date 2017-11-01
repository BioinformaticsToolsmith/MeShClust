all: bin/Red.o bin/meshclust

bin/Red.o:
	mkdir -p bin
	mkdir -p bin/exception
	mkdir -p bin/nonltr
	mkdir -p bin/utility
	$(MAKE) -C src
bin/meshclust: bin/Red.o
	$(MAKE) -C src/cluster
	cp src/cluster/meshclust bin

clean:
	$(MAKE) clean -C src
	$(MAKE) clean -C src/cluster
	$(RM) -r bin

rebuild: clean all
.PHONY: all clean
