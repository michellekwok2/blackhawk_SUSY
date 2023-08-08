.KEEP_STATE:

#
VERSION = v2.0
YEAR = 2021

# Choose your compilers here (usually gcc on Linux systems):
CC = gcc
CFLAGS= -O3 -pipe -fomit-frame-pointer
#CFLAGS_MP= -fopenmp

#CC = clang
#CFLAGS= -O3 -pipe -fomit-frame-pointer
#CFLAGS_MP= -fopenmp

#CC = icc
#CFLAGS = -O3
#CFLAGS_MP= -qopenmp

MAKE = make
#MAKE_MP = -j$(nproc)

AR = ar
ARFLAGS = rcs
#ARFLAGS = rcsU


.SUFFIXES:	.o .c .h
.PRECIOUS:	.c .h libblackhawk.a

CINCLUDE= -I./src -L./src

all: libblackhawk.a
	@case `uname` in \
	   Linux) RANL=;;\
	   OSF1) CFLAGS="$(CFLAGS) $(CFLAGS_MP) -ieee";;\
	   *) RANL="ranlib libnr.a";;\
	   esac;\
	echo ' ';\
   	echo 'Please run "make name" to compile "name.c" into "name.x"';\
	echo ' '

%.c:: %.c libblackhawk.a
	$(CC) -c $(CFLAGS) $@;
	$(CC) -o $*.x $(CFLAGS) $(CFLAGS_MP) $(CINCLUDE) $*.o -lblackhawk -lm;
	@rm -f $*.o;
	@touch $*.x

%:: %.c libblackhawk.a
	$(CC) -c $(CFLAGS) $*.c;
	$(CC) -o $*.x $(CFLAGS) $(CFLAGS_MP) $(CINCLUDE) $*.o -lblackhawk -lm;
	@rm -f $*.o;
	@touch $*.x

BlackHawk_tot: BlackHawk_tot.c libblackhawk.a

BlackHawk_inst: BlackHawk_inst.c libblackhawk.a

clean:
	rm -rf tmp *.x;
	@echo > src/FlagsForMake;
	$(MAKE) -C src/ clean
	
distclean: 
	rm -rf tmp *.x;
	@echo > src/FlagsForMake;
	$(MAKE) -C src/ distclean
	
libblackhawk.a: 
	@echo;
	@echo BlackHawk $(VERSION) - A. Arbey \& J. Auffinger $(YEAR);
	@echo;
	@echo CC = $(CC) > src/FlagsForMake;\
	echo CFLAGS = $(CFLAGS) >> src/FlagsForMake;\
	echo CFLAGS_MP = $(CFLAGS_MP) >> src/FlagsForMake;\
	echo MAKE = $(MAKE) >> src/FlagsForMake;\
	echo AR = $(AR) >> src/FlagsForMake;\
	echo ARFLAGS = $(ARFLAGS) >> src/FlagsForMake;
	$(MAKE) $(MAKE_mp) -C src/ libblackhawk.a

save: 
	rm -f blackhawk_$(VERSION).tar.bz2;\
	mkdir blackhawk_$(VERSION);\
	cp -p README.txt blackhawk_$(VERSION)/;\
	cp -p BlackHawk_tot.c blackhawk_$(VERSION)/;\
	cp -p BlackHawk_inst.c blackhawk_$(VERSION)/;\
	cp -p Makefile blackhawk_$(VERSION)/;\
	cp -p parameters.txt blackhawk_$(VERSION)/;\
	mkdir blackhawk_$(VERSION)/src;\
	cp -p src/*.h blackhawk_$(VERSION)/src/;\
	cp -p src/*.c blackhawk_$(VERSION)/src/;\
	cp -p src/Makefile blackhawk_$(VERSION)/src/;\
	cp -rp src/tables blackhawk_$(VERSION)/src/;\
	mkdir blackhawk_$(VERSION)/manual;\
	cp -p manual/*.pdf blackhawk_$(VERSION)/manual/;\
	cp -rp scripts blackhawk_$(VERSION)/;\
	mkdir blackhawk_$(VERSION)/results;\
	tar cjvf blackhawk_$(VERSION).tar.bz2 blackhawk_$(VERSION);\
	rm -rf blackhawk_$(VERSION)
