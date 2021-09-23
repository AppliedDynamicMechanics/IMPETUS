CXX=mpiCC
INCLUDE=-Ilib/include
LDLIBS=-Llib -lvamde-vamde-p

bin/vamde-run.exe: object/main.o
	$(CXX) $(INCLUDE) $< -o bin/impetus-run.exe $(LDLIBS)

object/main.o: src/main.cpp
	$(CXX) -c -w $(INCLUDE) $< -o $@

clean:
	rm -f bin/* object/*
