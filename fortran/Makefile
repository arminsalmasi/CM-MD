FC = gfortran
FFLAGS = -O2 -Wall
TEST_DIR = tests

OBJS = md_global.o md_utilities.o

all: md.out

md.out: md_global.o md_utilities.o md_monitor.o main.o
	$(FC) $(FFLAGS) -o $@ $^

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

test_calcKin: $(OBJS) $(TEST_DIR)/test_do_calcKin.f90
	$(FC) $(FFLAGS) -o $(TEST_DIR)/test_do_calcKin $(TEST_DIR)/test_do_calcKin.f90 $(OBJS)
	./$(TEST_DIR)/test_do_calcKin

test: test_calcKin

clean:
	rm -f *.o *.mod md.out $(TEST_DIR)/test_do_calcKin
