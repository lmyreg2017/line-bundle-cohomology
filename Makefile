headers := cohomology.h chomp.h
objects := cohomology.o scanbox.o chomp.o
program := cech_cohomology

$(program): $(objects)
	gcc $(objects) -o $@

%.o: %.c $(headers)
	gcc -c $< -o $@ -Wall -Wextra -O2 -g

clean:
	rm -f $(program) $(objects)

all: cech_cohomology